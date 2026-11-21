# GSK-3β structure-based virtual screen — CADD mini project

Group "Breaking Bad". Assignment 1 mechanism: the amyloid–tau cascade in
Alzheimer's disease. Target: **GSK-3β (UniProt P49841)**, node 4 of that
cascade — the kinase that turns an upstream amyloid signal into the tau
lesion that actually tracks with dementia.

## What is here

| path | contents |
|---|---|
| `CADD mini project.pdf` | the assignment brief |
| `Assignment1_Alzheimers_Amyloid_Tau.pptx` | the Assignment 1 deck this builds on |
| `project/07_report/` | **the report** (`CADD_miniproject_GSK3B_report.docx` / `.pdf`, 10 pages) |
| `project/scripts/` | the whole pipeline, 24 scripts numbered in run order |
| `project/01_target/` … `project/06_analysis/` | every intermediate: structures, alignments, models, ligands, docking runs, analysis |
| `project/README.md` | folder map and how to reproduce each stage |

Not in the repo: `tools/` (AutoDock Vina, P2Rank, MODELLER installer, the
PyMOL wheel — 582 MB of third-party binaries) and the submission zip, which
`project/scripts/13_package.py` rebuilds from the tracked files.

## Headline results

- **Model** — three GSK-3β templates (1Q5K, 1J1B, 4ACC), MODELLER 10.7, 5 models.
  Best model: Cα RMSD **0.429 Å** to 1Q5K, **98.6 %** of residues in favoured
  Ramachandran regions, GA341 1.00.
- **Site** — P2Rank's top-ranked pocket (score 9.54) is the ATP cleft.
  Cavity 397 Å³, 74 % of the co-crystal ligand buried, 70 % carbon lining.
  Box: centre (22.328, 24.372, 9.159) Å, 26 × 26 × 26 Å.
- **Screen** — 138 ligands (11 known actives + 80 CNS-like approved drugs +
  48 natural products) into both the model and the 1Q5K crystal, AutoDock
  Vina 1.2.5, exhaustiveness 32, seed 20241121.
- **Candidates** — tofacitinib, hesperetin, luteolin.

## Two results the report does not bury

**The redocking control fails.** Vina ranks AR-A014418 end-for-end reversed.
Locally minimising the crystallographic pose and scoring it gives
−6.58 kcal/mol against −7.67 for the reversed pose, so the error is in the
scoring function, not the search. With the ten conserved pocket waters the
crystal orientation is sampled (1.33 Å) but ranked twelfth.

**Enrichment is weak** — ROC-AUC 0.645 on the model, 0.640 on the crystal,
and the score does not track potency across three orders of magnitude
(Spearman ρ = 0.155). The screen is a triage filter, not a pose predictor or
an affinity model, and §5–7 of the report are written on that basis.

## Reproducing

Needs Python 3.13, MODELLER 10.7, AutoDock Vina 1.2.5, P2Rank 2.5 and
PyMOL (open source). `project/README.md` lists the scripts in order and the
exact settings.
