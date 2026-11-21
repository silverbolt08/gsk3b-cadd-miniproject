# CADD mini project — GSK-3β virtual screen

Group "Breaking Bad" · Assignment 1 mechanism: the amyloid–tau cascade in
Alzheimer's disease · target: **GSK-3β (UniProt P49841)**, node 4 of the
cascade (tau hyperphosphorylation).

The report is `07_report/CADD_miniproject_GSK3B_report.docx`. This file
explains what is in the data folder and how to regenerate all of it.

## How to reproduce

Scripts are in `scripts/`, numbered in run order. From `scripts/`:

```
python 01_fetch_actives.py          # ChEMBL: GSK-3b activities + 11 positive controls
python 02a_modeller_alignment.py    # salign: 3 templates + target -> gsk3b-mult.ali
python 02b_modeller_build.py        # AutoModel: 5 models, DOPE / DOPE-HR / GA341
python 02c_dope_profile.py          # per-residue DOPE profile
python 02d_validate_model.py        # RMSD vs 1Q5K, Ramachandran, figures 2 and 3
python 02e_template_table.py        # template search summary table
python 03_build_library.py          # 80 CNS-like approved drugs + 48 natural products
python 04_prepare_ligands.py        # pH 7.4 -> 3D -> PDBQT for all 139 ligands
python 05a_align_model.py           # put the model in the 1Q5K frame
#   P2Rank (see below) on 02_structure/models/gsk3b_model_aligned.pdb
python 05b_define_site.py           # pocket + search box -> 03_pocket/box.json
python 05c_prepare_receptors.py     # PDB2PQR pH 7.4 + Meeko -> receptor PDBQTs
python 06_dock.py --receptor model     # the screen
python 06_dock.py --receptor crystal   # cross-check on 1Q5K chain A
python 07_redock_validation.py      # redock AR-A014418, RMSD vs the crystal pose
python 08_complexes_and_plip.py     # complexes + PLIP interaction profiling
python 09_analyse.py                # ranking, ROC/EF/BEDROC, figures 4-7
python 11_fingerprint_figure.py     # interaction fingerprint + hinge summary
python 11b_shortlist.py             # apply the four selection criteria
python 15_pymol_figures.py          # candidate binding-mode figures (PyMOL, ray-traced)
python 15b_pymol_overview.py        # model overview + redocking overlay (PyMOL)
python 05d_pocket_druggability.py   # cavity volume, burial, hydrophobicity
python write_candidates.py          # the three candidate write-ups
python 12_build_report.py           # assemble the Word report
python 13_package.py                # zip the data folder
```

P2Rank is run separately (it is a Java tool):

```
java -Xmx4g -cp p2rank_2.5/bin/p2rank.jar;p2rank_2.5/bin/lib/* \
     cz.siret.prank.program.Main predict \
     -f 02_structure/models/gsk3b_model_aligned.pdb -o 03_pocket/p2rank
```

MODELLER is imported through `scripts/modeller_env.py`, which puts the
installed MODELLER 10.7 (`C:\Program Files\Modeller10.7`) on `sys.path` for
the system Python 3.13.

## Folder map

| path | contents |
|---|---|
| `01_target/` | UniProt P49841 FASTA/JSON, target rationale with references |
| `02_structure/templates/` | 1Q5K, 1J1B, 4ACC, 1PYX as downloaded; `1Q5K_A.pdb` (chain A, protein only) and `1Q5K_A_wat.pdb` (with the 10 conserved waters) |
| `02_structure/alignment/` | BLAST/RCSB search output, `templates.ali/.pap`, `gsk3b-mult.ali` (the modelling alignment), `template_table.csv` |
| `02_structure/models/` | the 5 MODELLER models, `model_scores.csv`, `gsk3b_model_best.pdb`, `gsk3b_model_aligned.pdb` (1Q5K frame), `..._aligned_wat.pdb` |
| `02_structure/validation/` | DOPE profiles, `validation_summary.json/.md` |
| `03_pocket/` | P2Rank output, `box.json` (centre, size, pocket residues), `pocket_residues.csv`, `druggability.json/.md` (cavity volume, burial, hydrophobicity) |
| `04_ligands/` | `actives/known_actives.csv` + the full ChEMBL activity table, `library/screening_library.csv`, `sdf/`, `pdbqt/`, `ligand_manifest.csv` |
| `05_docking/receptor/` | protonated receptors (`*_H.pdb`, `.pqr`) and the rigid `*.pdbqt` |
| `05_docking/model/`, `05_docking/crystal/` | per-ligand `configs/`, `logs/`, `poses/` and `scores.csv` |
| `05_docking/` | `ref_ligand_TMU.pdb` (crystal pose of AR-A014418), `redock_validation.json/.md` |
| `06_analysis/complexes/` | final receptor–ligand complexes (PDB) for the controls and top hits |
| `06_analysis/plip/` | `interactions.csv`, `fingerprint.csv`, `interaction_geometry.json` |
| `06_analysis/figures/pymol/` | the .pml scripts behind every structural figure |
| `06_analysis/` | `ranking_model.csv`, `enrichment.json`, `hinge_summary.json`, `shortlist.csv` (the four selection criteria, pass/fail per compound), `structural_alerts.csv` (PAINS/Brenk on the top 40), `figures/` |
| `07_report/` | the report, `candidates.json` (the judgement calls that go into §6.3–8) |

## Key numbers

- model: C-alpha RMSD 0.43 Å to 1Q5K chain A, 98.6 % Ramachandran favoured, GA341 1.00
- pocket: P2Rank rank 1, score 9.54, probability 0.51 — the ATP cleft
- box: centre (22.328, 24.372, 9.159) Å, 26 × 26 × 26 Å, 1Q5K frame
- docking: Vina 1.2.5, exhaustiveness 32, 12 modes, seed 20241121
- **redocking control fails**: the top pose of AR-A014418 is the reversed
  orientation; the crystal pose is sampled (1.3 Å) only with the conserved
  waters present, and is ranked 12th. See §5.1 of the report.

## Provenance

Structures from RCSB PDB; sequence from UniProt; activities and structures
of all ligands from ChEMBL (accessed November 2026). Nothing in the ligand
set was drawn by hand.
