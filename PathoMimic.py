"""
PathoMimic: batch Foldseek host structural-mimicry screen.

Run directly in Spyder or with Python. The program:

1. Creates Effector_Structure and Results beside this script, then finds every
   .pdb, .cif, and .mmcif file in Effector_Structure.
2. Resolves one or more host species to NCBI taxonomy IDs.
3. Runs a separate Foldseek search for every effector-host pair against the
   current AlphaFold/Proteome database and, optionally, PDB100.
4. Creates one result folder per effector, with host-specific and combined
   CSV files, raw JSON, candidate target FASTA sequences, and a summary.

Only the third-party ``requests`` package is required. It is normally included
with Anaconda/Spyder. If it is missing, run this in the Anaconda Prompt:

    conda install requests

Scientific interpretation
-------------------------
This is a hypothesis-generating structural resemblance screen. A Foldseek hit
is not, by itself, evidence of molecular mimicry. High-ranking candidates need
follow-up assessment of fold ubiquity, domain architecture, surface chemistry,
binding-site similarity, biological localisation, and experimental function.
The score produced here is an explicitly heuristic ranking score, not a
calibrated probability of mimicry.
"""

from __future__ import annotations

import csv
import hashlib
import json
import math
import re
import sys
import time
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Dict, Iterable, List, Optional, Sequence, Tuple
from urllib.parse import quote

try:
    import requests
except ImportError as exc:
    raise SystemExit(
        "PathoMimic requires the 'requests' package. Open Anaconda Prompt and "
        "run: conda install requests"
    ) from exc


# =============================================================================
# USER SETTINGS -- edit this section in Spyder
# =============================================================================

BASE_DIR = Path(__file__).resolve().parent
EFFECTOR_STRUCTURE_DIR = BASE_DIR / "Effector_Structure"
RESULTS_DIR = BASE_DIR / "Results"

# Enter scientific names, NCBI TaxIDs, or both. If left empty, the script asks
# for a comma-separated list in the Spyder console when it starts.
# Example:
# HOST_SPECIES = ["Arabidopsis thaliana", "Solanum tuberosum", "3702"]
HOST_SPECIES: List[str] = []

# Optional exact overrides. These bypass online name resolution and are useful
# if NCBI returns an ambiguous name. Keys must match entries in HOST_SPECIES.
# Example: HOST_TAXID_OVERRIDES = {"Arabidopsis thaliana": 3702}
HOST_TAXID_OVERRIDES: Dict[str, int] = {}

# The AlphaFold host-proteome search is always performed. PDB100 adds a search
# for experimentally determined host structures, although many hosts have few
# or no structures in PDB.
INCLUDE_PDB100 = True

# 3diaa is the recommended discovery mode. It is substantially faster than
# exhaustive TM-align mode and combines Foldseek 3Di and amino-acid information.
# Other server-supported values include "3di" and "tmalign".
SEARCH_MODE = "3diaa"

# Candidate classification thresholds. These are deliberately editable and
# should be treated as screening criteria, not biological proof thresholds.
MIN_FOLDSEEK_PROBABILITY = 0.80
MIN_QUERY_COVERAGE = 0.50
MIN_TARGET_COVERAGE_WHOLE_PROTEIN = 0.50
MAX_SEQUENCE_IDENTITY_PERCENT = 30.0
MAX_CANDIDATE_EVALUE = 1.0

# Resume/retry behaviour.
SKIP_IDENTICAL_COMPLETED_SEARCHES = True
POLL_SECONDS = 10
JOB_TIMEOUT_MINUTES = 45
SUBMISSION_PAUSE_SECONDS = 3
HTTP_TIMEOUT_SECONDS = 90
MAX_HTTP_RETRIES = 5
MAX_SUBMISSION_RETRIES = 6

# Optional notification/contact email passed to the Foldseek server.
CONTACT_EMAIL = ""

# Set True to validate inputs, resolve hosts, and show selected databases
# without submitting Foldseek jobs.
DRY_RUN = False


# =============================================================================
# PROGRAM CONSTANTS
# =============================================================================

SCRIPT_VERSION = "1.1.0"
FOLDSEEK_BASE_URL = "https://search.foldseek.com"
FOLDSEEK_API_URL = FOLDSEEK_BASE_URL + "/api"
NCBI_TAXON_SUGGEST_URL = (
    "https://api.ncbi.nlm.nih.gov/datasets/v2alpha/taxonomy/taxon_suggest/{}"
)

INPUT_SUFFIXES = {".pdb", ".cif", ".mmcif"}
CANDIDATE_CLASSES = {
    "candidate_whole_protein_mimic",
    "candidate_domain_mimic",
}

CSV_FIELDS = [
    "effector",
    "input_structure",
    "host_requested",
    "host_taxid",
    "database",
    "query_chain",
    "target",
    "uniprot_accession",
    "target_description",
    "target_taxid",
    "target_taxname",
    "sequence_identity_percent",
    "alignment_length",
    "query_length",
    "target_length",
    "query_coverage",
    "target_coverage",
    "foldseek_probability",
    "evalue",
    "bitscore",
    "query_start",
    "query_end",
    "target_start",
    "target_end",
    "global_mimicry_score",
    "domain_mimicry_score",
    "candidate_class",
    "ticket_id",
    "foldseek_results_url",
]


class PathoMimicError(RuntimeError):
    """Expected error that should be reported cleanly to the user."""


def utc_now() -> str:
    return datetime.now(timezone.utc).isoformat(timespec="seconds")


def safe_name(value: str, fallback: str = "unnamed") -> str:
    """Return a Windows-safe file/folder component."""
    cleaned = re.sub(r'[<>:"/\\|?*\x00-\x1f]+', "_", value.strip())
    cleaned = re.sub(r"\s+", "_", cleaned).strip(" ._")
    return cleaned[:150] or fallback


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def write_json(path: Path, data: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    with temporary.open("w", encoding="utf-8", newline="\n") as handle:
        json.dump(data, handle, indent=2, ensure_ascii=False)
        handle.write("\n")
    temporary.replace(path)


def write_text(path: Path, text: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    with temporary.open("w", encoding="utf-8", newline="\n") as handle:
        handle.write(text)
    temporary.replace(path)


def write_csv(path: Path, rows: Sequence[Dict[str, Any]], fields: Sequence[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    with temporary.open("w", encoding="utf-8-sig", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(fields), extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)
    temporary.replace(path)


def read_csv(path: Path) -> List[Dict[str, str]]:
    if not path.exists():
        return []
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle))


def validate_structure_file(path: Path) -> None:
    if not path.is_file() or path.stat().st_size == 0:
        raise PathoMimicError(f"Structure file is empty or unavailable: {path}")

    sample = path.read_text(encoding="utf-8", errors="replace")[:2_000_000]
    suffix = path.suffix.lower()
    if suffix == ".pdb" and not re.search(r"^(ATOM  |HETATM)", sample, re.MULTILINE):
        raise PathoMimicError(f"No ATOM/HETATM records found in PDB: {path.name}")
    if suffix in {".cif", ".mmcif"} and "_atom_site." not in sample:
        raise PathoMimicError(f"No _atom_site records found in mmCIF: {path.name}")


def request_json(
    session: requests.Session,
    method: str,
    url: str,
    *,
    retry_count: int = MAX_HTTP_RETRIES,
    **kwargs: Any,
) -> Dict[str, Any]:
    """Make a JSON HTTP request with bounded exponential backoff."""
    retry_statuses = {429, 500, 502, 503, 504}
    last_error: Optional[Exception] = None

    for attempt in range(retry_count + 1):
        try:
            response = session.request(
                method,
                url,
                timeout=HTTP_TIMEOUT_SECONDS,
                **kwargs,
            )
            if response.status_code in retry_statuses and attempt < retry_count:
                retry_after = response.headers.get("Retry-After", "")
                try:
                    wait_seconds = min(float(retry_after), 60.0)
                except (TypeError, ValueError):
                    wait_seconds = min(2.0 ** attempt * 5.0, 60.0)
                print(
                    f"  Server returned HTTP {response.status_code}; retrying in "
                    f"{wait_seconds:.0f} s..."
                )
                time.sleep(wait_seconds)
                continue

            response.raise_for_status()
            try:
                payload = response.json()
            except ValueError as exc:
                snippet = response.text[:500].replace("\n", " ")
                raise PathoMimicError(
                    f"Server returned non-JSON content from {url}: {snippet}"
                ) from exc
            if not isinstance(payload, dict):
                raise PathoMimicError(f"Unexpected JSON response from {url}")
            return payload
        except (requests.RequestException, PathoMimicError) as exc:
            last_error = exc
            if attempt >= retry_count:
                break
            wait_seconds = min(2.0 ** attempt * 5.0, 60.0)
            print(f"  Network/API error; retrying in {wait_seconds:.0f} s: {exc}")
            time.sleep(wait_seconds)

    raise PathoMimicError(f"Request failed after retries: {last_error}")


def get_user_hosts() -> List[str]:
    hosts = [str(item).strip() for item in HOST_SPECIES if str(item).strip()]
    if hosts:
        return hosts

    entered = input(
        "Enter host scientific names or NCBI TaxIDs, separated by commas:\n> "
    ).strip()
    hosts = [item.strip() for item in entered.split(",") if item.strip()]
    if not hosts:
        raise PathoMimicError("No host species were supplied.")
    return hosts


def resolve_host_taxid(
    session: requests.Session, host: str
) -> Tuple[str, int, str]:
    """Return requested label, TaxID, and resolved scientific name."""
    if host in HOST_TAXID_OVERRIDES:
        taxid = int(HOST_TAXID_OVERRIDES[host])
        if taxid <= 0:
            raise PathoMimicError(f"Invalid TaxID override for {host!r}: {taxid}")
        return host, taxid, host

    if host.isdigit():
        taxid = int(host)
        if taxid <= 0:
            raise PathoMimicError(f"Invalid NCBI TaxID: {host}")
        return host, taxid, f"TaxID_{taxid}"

    url = NCBI_TAXON_SUGGEST_URL.format(quote(host, safe=""))
    payload = request_json(
        session,
        "GET",
        url,
        params={"tax_rank_filter": "higher_taxon"},
    )
    suggestions = payload.get("sci_name_and_ids", [])
    if not isinstance(suggestions, list):
        suggestions = []

    wanted = " ".join(host.casefold().split())
    exact = []
    for item in suggestions:
        if not isinstance(item, dict):
            continue
        sci_name = " ".join(str(item.get("sci_name", "")).casefold().split())
        matched = " ".join(str(item.get("matched_term", "")).casefold().split())
        if sci_name == wanted or matched == wanted:
            exact.append(item)

    if not exact:
        alternatives = ", ".join(
            str(item.get("sci_name", ""))
            for item in suggestions[:5]
            if isinstance(item, dict)
        )
        extra = f" Suggestions: {alternatives}." if alternatives else ""
        raise PathoMimicError(
            f"NCBI did not return an exact taxonomy match for {host!r}.{extra} "
            "Use the accepted scientific name or add HOST_TAXID_OVERRIDES."
        )

    # Prefer an exact species-rank match, then the first exact textual match.
    chosen = next(
        (item for item in exact if str(item.get("rank", "")).upper() == "SPECIES"),
        exact[0],
    )
    try:
        taxid = int(chosen["tax_id"])
    except (KeyError, TypeError, ValueError) as exc:
        raise PathoMimicError(f"NCBI returned no usable TaxID for {host!r}") from exc
    return host, taxid, str(chosen.get("sci_name", host))


def resolve_all_hosts(
    session: requests.Session, host_entries: Sequence[str]
) -> List[Dict[str, Any]]:
    resolved: List[Dict[str, Any]] = []
    seen_taxids = set()
    for host in host_entries:
        requested, taxid, scientific_name = resolve_host_taxid(session, host)
        if taxid in seen_taxids:
            print(f"Skipping duplicate host TaxID {taxid} ({requested}).")
            continue
        seen_taxids.add(taxid)
        resolved.append(
            {
                "requested": requested,
                "taxid": taxid,
                "scientific_name": scientific_name,
            }
        )
        print(f"Host: {requested} -> {scientific_name} (NCBI TaxID {taxid})")
    return resolved


def version_key(value: str) -> Tuple[int, ...]:
    numbers = re.findall(r"\d+", value)
    return tuple(int(number) for number in numbers) or (0,)


def choose_database(
    databases: Sequence[Dict[str, Any]], name: str, required: bool
) -> Optional[Dict[str, Any]]:
    candidates = [
        db
        for db in databases
        if db.get("name") == name
        and db.get("status") == "COMPLETE"
        and bool(db.get("taxonomy"))
        and not bool(db.get("motif"))
        and not bool(db.get("rna"))
        and not bool(db.get("interface"))
    ]
    if not candidates:
        if required:
            raise PathoMimicError(
                f"Foldseek currently exposes no usable taxonomy-enabled {name} database."
            )
        print(f"Warning: optional Foldseek database {name} is unavailable; skipping it.")
        return None
    return max(candidates, key=lambda item: version_key(str(item.get("version", ""))))


def discover_databases(session: requests.Session) -> List[Dict[str, Any]]:
    payload = request_json(session, "GET", FOLDSEEK_API_URL + "/databases")
    available = payload.get("databases", [])
    if not isinstance(available, list):
        raise PathoMimicError("Foldseek returned an invalid database catalogue.")

    chosen = [choose_database(available, "AlphaFold/Proteome", required=True)]
    if INCLUDE_PDB100:
        chosen.append(choose_database(available, "PDB100", required=False))
    resolved = [item for item in chosen if item is not None]
    for item in resolved:
        print(
            f"Database: {item['name']} {item.get('version', '')} "
            f"(server id: {item['path']})"
        )
    return resolved


def submit_foldseek_job(
    session: requests.Session,
    structure_path: Path,
    database_paths: Sequence[str],
    host_taxid: int,
) -> Dict[str, Any]:
    data: List[Tuple[str, str]] = [
        ("mode", SEARCH_MODE),
        ("email", CONTACT_EMAIL),
        ("taxfilter", str(host_taxid)),
    ]
    data.extend(("database[]", path) for path in database_paths)

    for attempt in range(MAX_SUBMISSION_RETRIES + 1):
        try:
            # Reopen the structure on every retry. A requests file handle is at
            # EOF after an attempted upload and must not be silently reused.
            with structure_path.open("rb") as handle:
                files = {"q": (structure_path.name, handle, "application/octet-stream")}
                payload = request_json(
                    session,
                    "POST",
                    FOLDSEEK_API_URL + "/ticket",
                    files=files,
                    data=data,
                    retry_count=0,
                )
        except PathoMimicError:
            if attempt >= MAX_SUBMISSION_RETRIES:
                raise
            wait_seconds = min(10 * (attempt + 1), 60)
            print(f"  Submission failed; retrying in {wait_seconds} s...")
            time.sleep(wait_seconds)
            continue

        status = str(payload.get("status", "UNKNOWN")).upper()
        if status in {"PENDING", "RUNNING", "COMPLETE"} and payload.get("id"):
            return payload
        if status in {"RATELIMIT", "MAINTENANCE"} and attempt < MAX_SUBMISSION_RETRIES:
            wait_seconds = min(15 * (attempt + 1), 60)
            print(f"  Foldseek status {status}; retrying in {wait_seconds} s...")
            time.sleep(wait_seconds)
            continue
        raise PathoMimicError(f"Foldseek rejected the job: {payload}")

    raise PathoMimicError("Foldseek submission retry limit was reached.")


def wait_for_ticket(session: requests.Session, ticket: str) -> Dict[str, Any]:
    deadline = time.monotonic() + JOB_TIMEOUT_MINUTES * 60
    last_status = ""
    last_report = 0.0

    while time.monotonic() < deadline:
        payload = request_json(
            session,
            "GET",
            FOLDSEEK_API_URL + f"/ticket/{quote(ticket, safe='')}",
        )
        status = str(payload.get("status", "UNKNOWN")).upper()
        now = time.monotonic()
        if status != last_status or now - last_report >= 60:
            print(f"  Ticket {ticket}: {status}")
            last_status = status
            last_report = now

        if status == "COMPLETE":
            return payload
        if status in {"ERROR", "UNKNOWN"}:
            raise PathoMimicError(f"Foldseek ticket ended with status {status}: {payload}")
        if status in {"RATELIMIT", "MAINTENANCE"}:
            time.sleep(min(max(POLL_SECONDS, 15), 60))
        else:
            time.sleep(min(max(POLL_SECONDS, 2), 60))

    raise PathoMimicError(
        f"Foldseek ticket {ticket} did not finish within {JOB_TIMEOUT_MINUTES} minutes."
    )


def fetch_results(session: requests.Session, ticket: str) -> Dict[str, Any]:
    return request_json(
        session,
        "GET",
        FOLDSEEK_API_URL + f"/result/{quote(ticket, safe='')}/0",
        retry_count=MAX_HTTP_RETRIES + 2,
    )


def as_float(value: Any, default: float = 0.0) -> float:
    try:
        return float(value)
    except (TypeError, ValueError):
        return default


def as_int(value: Any, default: int = 0) -> int:
    try:
        return int(value)
    except (TypeError, ValueError):
        return default


def parse_target_header(target: str) -> Tuple[str, str]:
    """Extract an AlphaFold UniProt accession and description when possible."""
    match = re.match(
        r"^AF-([A-Za-z0-9]+)-F\d+-model_v\d+(?:\.pdb(?:\.gz)?)?\s*(.*)$",
        target.strip(),
    )
    if match:
        return match.group(1), match.group(2).strip()
    first, _, remainder = target.strip().partition(" ")
    return "", remainder.strip() if remainder else first


def classify_hit(
    probability: float,
    query_coverage: float,
    target_coverage: float,
    sequence_identity: float,
    evalue: float,
) -> str:
    if probability < MIN_FOLDSEEK_PROBABILITY or query_coverage < MIN_QUERY_COVERAGE:
        return "weak_or_partial_resemblance"
    if sequence_identity > MAX_SEQUENCE_IDENTITY_PERCENT:
        return "likely_conventional_homologue"
    if evalue > MAX_CANDIDATE_EVALUE:
        return "lower_confidence_structural_resemblance"
    if target_coverage >= MIN_TARGET_COVERAGE_WHOLE_PROTEIN:
        return "candidate_whole_protein_mimic"
    return "candidate_domain_mimic"


def flatten_results(
    raw: Dict[str, Any],
    *,
    effector_name: str,
    input_structure: Path,
    host: Dict[str, Any],
    ticket: str,
) -> List[Dict[str, Any]]:
    rows: List[Dict[str, Any]] = []
    results_url = f"{FOLDSEEK_BASE_URL}/result/{ticket}/0"

    for database_result in raw.get("results", []):
        if not isinstance(database_result, dict):
            continue
        database = str(database_result.get("db", ""))
        groups = database_result.get("alignments", [])
        if isinstance(groups, dict):
            groups = [[groups]]
        elif groups and isinstance(groups, list) and isinstance(groups[0], dict):
            groups = [groups]

        for group in groups or []:
            if isinstance(group, dict):
                group = [group]
            if not isinstance(group, list):
                continue
            for hit in group:
                if not isinstance(hit, dict):
                    continue
                alignment_length = as_int(hit.get("alnLength"))
                query_length = as_int(hit.get("qLen"))
                target_length = as_int(hit.get("dbLen"))
                query_coverage = (
                    alignment_length / query_length if query_length > 0 else 0.0
                )
                target_coverage = (
                    alignment_length / target_length if target_length > 0 else 0.0
                )
                probability = as_float(hit.get("prob"))
                sequence_identity = as_float(hit.get("seqId"))
                evalue = as_float(hit.get("eval"), default=float("inf"))
                novelty = max(0.0, 1.0 - sequence_identity / 100.0)
                global_score = (
                    100.0
                    * probability
                    * math.sqrt(max(0.0, query_coverage * target_coverage))
                    * novelty
                )
                domain_score = 100.0 * probability * query_coverage * novelty
                category = classify_hit(
                    probability,
                    query_coverage,
                    target_coverage,
                    sequence_identity,
                    evalue,
                )
                target = str(hit.get("target", ""))
                accession, description = parse_target_header(target)
                rows.append(
                    {
                        "effector": effector_name,
                        "input_structure": str(input_structure),
                        "host_requested": host["requested"],
                        "host_taxid": host["taxid"],
                        "database": database,
                        "query_chain": str(hit.get("query", "")),
                        "target": target,
                        "uniprot_accession": accession,
                        "target_description": description,
                        "target_taxid": hit.get("taxId", ""),
                        "target_taxname": hit.get("taxName", ""),
                        "sequence_identity_percent": round(sequence_identity, 3),
                        "alignment_length": alignment_length,
                        "query_length": query_length,
                        "target_length": target_length,
                        "query_coverage": round(query_coverage, 4),
                        "target_coverage": round(target_coverage, 4),
                        "foldseek_probability": round(probability, 6),
                        "evalue": evalue,
                        "bitscore": hit.get("score", ""),
                        "query_start": hit.get("qStartPos", ""),
                        "query_end": hit.get("qEndPos", ""),
                        "target_start": hit.get("dbStartPos", ""),
                        "target_end": hit.get("dbEndPos", ""),
                        "global_mimicry_score": round(global_score, 3),
                        "domain_mimicry_score": round(domain_score, 3),
                        "candidate_class": category,
                        "ticket_id": ticket,
                        "foldseek_results_url": results_url,
                        "_target_sequence": str(hit.get("tSeq", "")),
                    }
                )

    rows.sort(
        key=lambda row: (
            as_float(row.get("domain_mimicry_score")),
            as_float(row.get("global_mimicry_score")),
            as_float(row.get("foldseek_probability")),
        ),
        reverse=True,
    )
    return rows


def candidate_rows(rows: Sequence[Dict[str, Any]]) -> List[Dict[str, Any]]:
    return [row for row in rows if row.get("candidate_class") in CANDIDATE_CLASSES]


def write_candidate_fasta(path: Path, rows: Sequence[Dict[str, Any]]) -> None:
    seen = set()
    lines: List[str] = []
    for index, row in enumerate(rows, start=1):
        sequence = re.sub(r"[^A-Za-z]", "", str(row.get("_target_sequence", ""))).upper()
        target = str(row.get("target", ""))
        unique_key = (str(row.get("database", "")), target, sequence)
        if not sequence or unique_key in seen:
            continue
        seen.add(unique_key)
        accession = str(row.get("uniprot_accession", "")) or f"candidate_{index}"
        header = "|".join(
            [
                safe_name(accession),
                safe_name(str(row.get("host_requested", "host"))),
                safe_name(str(row.get("candidate_class", "candidate"))),
                f"score={row.get('domain_mimicry_score', '')}",
            ]
        )
        lines.append(">" + header)
        lines.extend(sequence[pos : pos + 80] for pos in range(0, len(sequence), 80))
    write_text(path, "\n".join(lines) + ("\n" if lines else ""))


def settings_fingerprint(
    structure_hash: str,
    host: Dict[str, Any],
    databases: Sequence[Dict[str, Any]],
) -> str:
    settings = {
        "script_version": SCRIPT_VERSION,
        "structure_sha256": structure_hash,
        "host_taxid": host["taxid"],
        "database_paths": [db["path"] for db in databases],
        "database_versions": [db.get("version", "") for db in databases],
        "search_mode": SEARCH_MODE,
        "thresholds": {
            "min_probability": MIN_FOLDSEEK_PROBABILITY,
            "min_query_coverage": MIN_QUERY_COVERAGE,
            "min_target_coverage_whole": MIN_TARGET_COVERAGE_WHOLE_PROTEIN,
            "max_sequence_identity_percent": MAX_SEQUENCE_IDENTITY_PERCENT,
            "max_candidate_evalue": MAX_CANDIDATE_EVALUE,
        },
    }
    encoded = json.dumps(settings, sort_keys=True, separators=(",", ":")).encode("utf-8")
    return hashlib.sha256(encoded).hexdigest()


def write_effector_summary(
    path: Path,
    structure: Path,
    host_stats: Sequence[Dict[str, Any]],
    failures: Sequence[Dict[str, Any]],
) -> None:
    lines = [
        f"# PathoMimic summary: {structure.stem}",
        "",
        f"- Input structure: `{structure}`",
        f"- Generated: {utc_now()}",
        f"- Search mode: `{SEARCH_MODE}`",
        "",
        "## Results by host",
        "",
        "| Host | TaxID | Total hits | Candidate mimics | Foldseek result |",
        "|---|---:|---:|---:|---|",
    ]
    for item in host_stats:
        result_link = (
            f"[open]({item['results_url']})" if item.get("results_url") else "-"
        )
        lines.append(
            f"| {item['host']} | {item['taxid']} | {item['total_hits']} | "
            f"{item['candidate_hits']} | {result_link} |"
        )

    if failures:
        lines.extend(["", "## Failures", ""])
        for failure in failures:
            lines.append(f"- **{failure['host']}**: {failure['error']}")

    lines.extend(
        [
            "",
            "## Interpretation",
            "",
            "`candidate_whole_protein_mimic` means the hit passed the configured "
            "Foldseek probability, query coverage, target coverage, sequence identity, "
            "and E-value screens.",
            "",
            "`candidate_domain_mimic` passed the same screen but covers less than the "
            "configured fraction of the host target, so it may represent a local/domain "
            "resemblance rather than whole-protein mimicry.",
            "",
            "The global and domain mimicry scores are heuristic ranking values. They are "
            "not probabilities and should not be used as evidence of mimicry without "
            "follow-up structural and biological validation.",
            "",
        ]
    )
    write_text(path, "\n".join(lines))


def make_session() -> requests.Session:
    session = requests.Session()
    contact = CONTACT_EMAIL.strip() or "not-provided"
    session.headers.update(
        {
            "User-Agent": f"PathoMimic/{SCRIPT_VERSION} (contact: {contact})",
            "Accept": "application/json",
        }
    )
    return session


def effector_output_names(structures: Sequence[Path]) -> Dict[Path, str]:
    counts: Dict[str, int] = {}
    for structure in structures:
        key = safe_name(structure.stem).casefold()
        counts[key] = counts.get(key, 0) + 1

    names: Dict[Path, str] = {}
    for structure in structures:
        stem = safe_name(structure.stem)
        if counts[stem.casefold()] > 1:
            stem += "_" + safe_name(structure.suffix.lstrip("."))
        names[structure] = stem + "_Foldseek_Mimicry"
    return names


def main() -> None:
    print(f"PathoMimic Foldseek screen v{SCRIPT_VERSION}")
    print("=" * 72)

    EFFECTOR_STRUCTURE_DIR.mkdir(parents=True, exist_ok=True)
    RESULTS_DIR.mkdir(parents=True, exist_ok=True)

    structures = sorted(
        (
            path
            for path in EFFECTOR_STRUCTURE_DIR.iterdir()
            if path.is_file() and path.suffix.lower() in INPUT_SUFFIXES
        ),
        key=lambda path: path.name.casefold(),
    )
    if not structures:
        raise PathoMimicError(
            f"No .pdb, .cif, or .mmcif files were found in {EFFECTOR_STRUCTURE_DIR}. "
            "Add structure files there and run PathoMimic again."
        )

    for structure in structures:
        validate_structure_file(structure)
    print(f"Found {len(structures)} structure file(s).")

    session = make_session()
    host_entries = get_user_hosts()
    hosts = resolve_all_hosts(session, host_entries)
    databases = discover_databases(session)
    database_paths = [str(db["path"]) for db in databases]

    if DRY_RUN:
        print("\nDRY_RUN is True. Inputs are valid; no Foldseek jobs were submitted.")
        return

    output_names = effector_output_names(structures)
    batch_rows: List[Dict[str, Any]] = []
    batch_failures: List[Dict[str, Any]] = []

    for structure_index, structure in enumerate(structures, start=1):
        effector_name = structure.stem
        effector_dir = RESULTS_DIR / output_names[structure]
        effector_dir.mkdir(parents=True, exist_ok=True)
        structure_hash = sha256_file(structure)
        print(
            f"\n[{structure_index}/{len(structures)}] Effector: {structure.name}\n"
            f"Results: {effector_dir}"
        )

        combined_rows: List[Dict[str, Any]] = []
        host_stats: List[Dict[str, Any]] = []
        effector_failures: List[Dict[str, Any]] = []

        for host_index, host in enumerate(hosts, start=1):
            host_label = host["scientific_name"]
            host_dir = effector_dir / "per_host" / (
                safe_name(host_label) + f"_taxid_{host['taxid']}"
            )
            host_dir.mkdir(parents=True, exist_ok=True)
            status_path = host_dir / "status.json"
            fingerprint = settings_fingerprint(structure_hash, host, databases)

            print(
                f"  [{host_index}/{len(hosts)}] Host: {host_label} "
                f"(TaxID {host['taxid']})"
            )

            if SKIP_IDENTICAL_COMPLETED_SEARCHES and status_path.exists():
                try:
                    prior_status = json.loads(status_path.read_text(encoding="utf-8"))
                except (OSError, json.JSONDecodeError):
                    prior_status = {}
                if (
                    prior_status.get("status") == "complete"
                    and prior_status.get("fingerprint") == fingerprint
                    and (host_dir / "all_hits.csv").exists()
                    and (host_dir / "raw_foldseek_result.json").exists()
                ):
                    # Re-flatten the raw result so private target sequences are
                    # retained when rebuilding the combined candidate FASTA.
                    try:
                        raw = json.loads(
                            (host_dir / "raw_foldseek_result.json").read_text(
                                encoding="utf-8"
                            )
                        )
                        prior_ticket = str(prior_status.get("ticket_id", ""))
                        rows = flatten_results(
                            raw,
                            effector_name=effector_name,
                            input_structure=structure,
                            host=host,
                            ticket=prior_ticket,
                        )
                    except (OSError, json.JSONDecodeError, TypeError, ValueError) as exc:
                        print(f"    Existing result is unreadable; rerunning: {exc}")
                    else:
                        candidates = candidate_rows(rows)
                        combined_rows.extend(rows)
                        results_url = str(prior_status.get("foldseek_results_url", ""))
                        host_stats.append(
                            {
                                "host": host_label,
                                "taxid": host["taxid"],
                                "total_hits": len(rows),
                                "candidate_hits": len(candidates),
                                "results_url": results_url,
                            }
                        )
                        print("    Reusing identical completed search.")
                        continue

            started_at = utc_now()
            try:
                ticket_payload = submit_foldseek_job(
                    session,
                    structure,
                    database_paths,
                    int(host["taxid"]),
                )
                ticket = str(ticket_payload["id"])
                print(f"    Submitted Foldseek ticket: {ticket}")
                if str(ticket_payload.get("status", "")).upper() != "COMPLETE":
                    wait_for_ticket(session, ticket)
                raw = fetch_results(session, ticket)

                write_json(host_dir / "raw_foldseek_result.json", raw)
                rows = flatten_results(
                    raw,
                    effector_name=effector_name,
                    input_structure=structure,
                    host=host,
                    ticket=ticket,
                )
                candidates = candidate_rows(rows)
                write_csv(host_dir / "all_hits.csv", rows, CSV_FIELDS)
                write_csv(host_dir / "candidate_mimics.csv", candidates, CSV_FIELDS)
                write_candidate_fasta(host_dir / "candidate_target_sequences.fasta", candidates)

                results_url = f"{FOLDSEEK_BASE_URL}/result/{ticket}/0"
                status = {
                    "status": "complete",
                    "script_version": SCRIPT_VERSION,
                    "started_at": started_at,
                    "completed_at": utc_now(),
                    "fingerprint": fingerprint,
                    "input_structure": str(structure),
                    "input_sha256": structure_hash,
                    "host": host,
                    "databases": databases,
                    "search_mode": SEARCH_MODE,
                    "ticket_id": ticket,
                    "foldseek_results_url": results_url,
                    "total_hits": len(rows),
                    "candidate_hits": len(candidates),
                }
                write_json(status_path, status)
                combined_rows.extend(rows)
                host_stats.append(
                    {
                        "host": host_label,
                        "taxid": host["taxid"],
                        "total_hits": len(rows),
                        "candidate_hits": len(candidates),
                        "results_url": results_url,
                    }
                )
                print(f"    Saved {len(rows)} hits; {len(candidates)} candidates.")
            except Exception as exc:  # isolate one failed host job from the batch
                error = f"{type(exc).__name__}: {exc}"
                failure = {
                    "effector": effector_name,
                    "input_structure": str(structure),
                    "host": host_label,
                    "host_taxid": host["taxid"],
                    "error": error,
                    "time": utc_now(),
                }
                effector_failures.append(failure)
                batch_failures.append(failure)
                write_json(
                    status_path,
                    {
                        "status": "failed",
                        "fingerprint": fingerprint,
                        "started_at": started_at,
                        "failed_at": utc_now(),
                        "failure": failure,
                    },
                )
                write_text(host_dir / "ERROR.txt", error + "\n")
                print(f"    FAILED: {error}")

            if SUBMISSION_PAUSE_SECONDS > 0:
                time.sleep(min(SUBMISSION_PAUSE_SECONDS, 60))

        combined_rows.sort(
            key=lambda row: (
                as_float(row.get("domain_mimicry_score")),
                as_float(row.get("global_mimicry_score")),
            ),
            reverse=True,
        )
        combined_candidates = candidate_rows(combined_rows)
        write_csv(effector_dir / "all_hits.csv", combined_rows, CSV_FIELDS)
        write_csv(
            effector_dir / "candidate_mimics.csv",
            combined_candidates,
            CSV_FIELDS,
        )
        write_candidate_fasta(
            effector_dir / "candidate_target_sequences.fasta",
            combined_candidates,
        )
        write_effector_summary(
            effector_dir / "summary.md",
            structure,
            host_stats,
            effector_failures,
        )
        write_json(
            effector_dir / "run_metadata.json",
            {
                "script_version": SCRIPT_VERSION,
                "generated_at": utc_now(),
                "input_structure": str(structure),
                "input_sha256": structure_hash,
                "hosts": hosts,
                "databases": databases,
                "search_mode": SEARCH_MODE,
                "thresholds": {
                    "min_foldseek_probability": MIN_FOLDSEEK_PROBABILITY,
                    "min_query_coverage": MIN_QUERY_COVERAGE,
                    "min_target_coverage_whole_protein": MIN_TARGET_COVERAGE_WHOLE_PROTEIN,
                    "max_sequence_identity_percent": MAX_SEQUENCE_IDENTITY_PERCENT,
                    "max_candidate_evalue": MAX_CANDIDATE_EVALUE,
                },
                "total_hits": len(combined_rows),
                "candidate_hits": len(combined_candidates),
                "failures": effector_failures,
            },
        )
        for stat in host_stats:
            batch_rows.append(
                {
                    "effector": effector_name,
                    "input_structure": str(structure),
                    "result_folder": str(effector_dir),
                    **stat,
                }
            )

    batch_fields = [
        "effector",
        "input_structure",
        "result_folder",
        "host",
        "taxid",
        "total_hits",
        "candidate_hits",
        "results_url",
    ]
    failure_fields = [
        "effector",
        "input_structure",
        "host",
        "host_taxid",
        "error",
        "time",
    ]
    write_csv(RESULTS_DIR / "PathoMimic_batch_summary.csv", batch_rows, batch_fields)
    write_csv(RESULTS_DIR / "PathoMimic_failures.csv", batch_failures, failure_fields)

    print("\n" + "=" * 72)
    print("PathoMimic batch finished.")
    print(f"Batch summary: {RESULTS_DIR / 'PathoMimic_batch_summary.csv'}")
    if batch_failures:
        print(
            f"{len(batch_failures)} host search(es) failed. Details: "
            f"{RESULTS_DIR / 'PathoMimic_failures.csv'}"
        )


if __name__ == "__main__":
    try:
        main()
    except KeyboardInterrupt:
        print("\nPathoMimic stopped by user. Completed host searches remain reusable.")
    except PathoMimicError as exc:
        print(f"\nERROR: {exc}")
        sys.exit(1)
