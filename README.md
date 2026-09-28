# PathoMimic

PathoMimic is a Python script for screening pathogen effector structures against host protein structures using the [Foldseek search server](https://search.foldseek.com/). It processes multiple effectors and multiple host species in one run, then ranks structural resemblances for further investigation.

**A structural match is a hypothesis, not evidence that an effector mimics a host protein.** The scores and candidate labels are screening aids; they are not probabilities or validated biological classifications.

## Requirements

- Python 3.9 or newer and an internet connection.
- The Python package `requests`.
- One or more effector structure files in PDB (`.pdb`) or mmCIF (`.cif`, `.mmcif`) format. Each file must contain atomic coordinates.

Install the dependency in the Python environment used to run the script:

```bash
python -m pip install requests
```

If you use Anaconda and Spyder, you can instead run `conda install requests` in the Anaconda Prompt for the environment used by Spyder.

## Quick start

1. Download `PathoMimic.py` into a folder where you have permission to create files.
2. Run it once with `python PathoMimic.py`, or open and run it in Spyder. It creates `Effector_Structure` and `Results` alongside the script. On this first run, a message explaining that no structure files were found is expected.
3. Put your `.pdb`, `.cif`, or `.mmcif` effector files directly inside `Effector_Structure`.
4. Run the script again. When prompted, enter one or more host scientific names or NCBI TaxIDs, separated by commas, for example `Arabidopsis thaliana, Solanum tuberosum`.
5. Open `Results/PathoMimic_batch_summary.csv` for the batch overview. Each effector has its own results folder.

The folders follow this layout:

```text
PathoMimic.py
Effector_Structure/
    effector_1.pdb
    effector_2.cif
Results/
    PathoMimic_batch_summary.csv
    PathoMimic_failures.csv
    effector_1_Foldseek_Mimicry/
        summary.md
        all_hits.csv
        candidate_mimics.csv
        candidate_target_sequences.fasta
        run_metadata.json
        per_host/
            Host_name_taxid_1234/
                all_hits.csv
                candidate_mimics.csv
                candidate_target_sequences.fasta
                raw_foldseek_result.json
                status.json
```

The script uses its own location for these folders, regardless of the current working directory. Input files must be directly within `Effector_Structure`, not in nested subfolders.

## Settings

Edit the **USER SETTINGS** section near the start of `PathoMimic.py` if you want to change the defaults:

| Setting | Default | Purpose |
|---|---|---|
| `HOST_SPECIES` | `[]` | Empty means prompt on each run. Set a list such as `["Arabidopsis thaliana", "Solanum tuberosum"]` to reuse hosts without prompting. TaxID strings also work. |
| `HOST_TAXID_OVERRIDES` | `{}` | Supply a TaxID for a scientific name if automatic NCBI resolution is ambiguous. |
| `INCLUDE_PDB100` | `True` | Include PDB100 when a suitable taxonomy enabled database is available. The AlphaFold/Proteome search is always included. |
| `SEARCH_MODE` | `"3diaa"` | Foldseek search mode. |
| `DRY_RUN` | `False` | Check structures, hosts, and available databases without submitting Foldseek searches. This still requires internet access. |
| `SKIP_IDENTICAL_COMPLETED_SEARCHES` | `True` | Reuse completed results when the structure, host, database versions, mode, and thresholds match. |
| `CONTACT_EMAIL` | `""` | Optional contact email sent to the Foldseek service. |

The probability, coverage, identity, and E-value thresholds are also editable in that section. Changes to those thresholds affect candidate classification and the cache fingerprint. PDB100 can have sparse representation for some hosts; an unavailable optional PDB100 database is skipped.

## Output and interpretation

Each effector folder contains `all_hits.csv` for every returned match, `candidate_mimics.csv` for hits passing the configured candidate screens, `candidate_target_sequences.fasta` for candidate sequences when the server supplies them, and a human readable `summary.md`. The `per_host` folders retain separate CSV files, raw Foldseek JSON, and search status. `run_metadata.json` records settings and counts. The batch CSV files collect results and failures across effectors.

The script calls a hit `candidate_whole_protein_mimic` when it passes the Foldseek probability, query coverage, target coverage, sequence identity, and E-value screens. A `candidate_domain_mimic` passes those screens except the whole protein target coverage criterion. The `global_mimicry_score` and `domain_mimicry_score` columns are heuristic sorting scores, **not** probabilities of mimicry. Other labels indicate weak or partial resemblance, likely conventional homology, or lower confidence structural resemblance.

Inspect the alignment, whether the fold is common, domain context, surface and binding site similarity, biological localisation, and experimental function before proposing mimicry. Foldseek and NCBI are external services: server availability, database releases, and returned hits can change over time.

## Troubleshooting

- **`No .pdb, .cif, or .mmcif files were found`**: Put coordinate files directly in the newly created `Effector_Structure` folder and rerun.
- **`No ATOM/HETATM records` or `No _atom_site records`**: Check that the file contains a structure in the indicated format, rather than a sequence or an empty placeholder.
- **`No host species were supplied`**: Enter hosts at the console prompt or populate `HOST_SPECIES` in the script. Spyder must have an interactive console available when the list is empty.
- **Taxonomy lookup fails**: Try an accepted scientific name or a numeric NCBI TaxID. You can specify an exact mapping in `HOST_TAXID_OVERRIDES`.
- **A host job fails**: See `Results/PathoMimic_failures.csv` and that host's `ERROR.txt`. Correct the issue and rerun; completed matching searches can be reused.
- **`requests` is missing**: Install it in the same Python environment that runs PathoMimic.

## Citation and service acknowledgement

If this tool contributes to a publication, cite the Foldseek method and any structure databases used, following their current citation guidance. PathoMimic sends uploaded structures and the selected host TaxID to the public Foldseek search service; host name resolution also contacts NCBI. Avoid submitting structures you cannot share with those services.
