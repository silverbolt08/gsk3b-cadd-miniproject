"""Sections 5-8 of the report (screen, analysis, candidates, limitations)."""
from docx.shared import Pt
from docx.enum.text import WD_ALIGN_PARAGRAPH


def tail(R, D):
    box = D["box"]; redock = D["redock"]; enrich = D["enrich"]
    ranking = D["ranking"]; cands = D["cands"]; hinge = D["hinge"]
    rank_by_id = {r["id"]: r for r in ranking}

    # ----------------------------------------------------------- Stage 5 --
    R.h1("5  The screen")

    R.h2("5.1  Protocol validation: redocking AR-A014418")
    rc, rm = redock["receptors"]["crystal"], redock["receptors"]["model"]
    R.p("Before docking anything unknown we docked something known. AR-A014418 is resolved "
        "in 1Q5K chain A; we rebuilt it from SMILES through our own ligand pipeline and "
        "docked it back into both receptors with the production settings, then measured the "
        "symmetry-corrected heavy-atom RMSD of each pose against the crystallographic one.")
    R.table(["receptor", "top-pose score (kcal/mol)", "top-pose RMSD (Å)",
             "best RMSD over 12 modes (Å)"],
            [["1Q5K crystal", f"{rc['top_pose_score']:.2f}", f"{rc['top_pose_rmsd']:.2f}",
              f"{rc['best_rmsd_any_mode']:.2f} (mode {rc['best_rmsd_mode']})"],
             ["homology model", f"{rm['top_pose_score']:.2f}", f"{rm['top_pose_rmsd']:.2f}",
              f"{rm['best_rmsd_any_mode']:.2f} (mode {rm['best_rmsd_mode']})"]],
            caption="Redocking the co-crystal ligand. The usual pass mark is a top pose "
                    "under 2 Å. Our protocol does not reach it.",
            widths=[3.6, 4.4, 4.0, 4.0])
    R.p("**This test fails, and the failure is informative.** Three diagnostics locate it:")
    R.bullet("*Not the box or the receptor.* Every pose lands in the cleft — the top "
             "pose’s centroid is 2.9 Å from the crystal ligand’s. The "
             "molecule is in the right place, **end-for-end reversed**: the methoxyphenyl "
             "sits where the nitrothiazole should be, 14 Å from its crystallographic "
             "position, while the internal conformation is right (1.5–1.9 Å "
             "after superposition).")
    R.bullet("*Not the search.* Raising exhaustiveness from 16 to 64 with 20 modes does "
             "not find the crystal pose. Adding back the ten conserved pocket waters does: "
             "the crystal orientation then appears at "
             f"{rc['best_rmsd_any_mode']:.2f} Å on the crystal receptor and "
             f"{rm['best_rmsd_any_mode']:.2f} Å on the model. The waters were kept for "
             "this reason.")
    R.bullet("*It is the scoring function.* Locally minimising the crystallographic pose "
             "itself with Vina (it moves 0.39 Å) and scoring it gives −6.58 "
             "kcal/mol against −7.67 for the reversed pose: Vina prefers the wrong "
             "orientation by ≈1.1 kcal/mol. Its scoring is charge-blind and rewards "
             "buried hydrophobic surface, and the real pose leaves the methoxyphenyl "
             "half-exposed near Arg141. Rescoring with Vinardo is worse.")
    R.figure("fig9_redock.png",
             "The redocking test. Crystallographic AR-A014418 in yellow; Vina’s "
             "top-ranked pose in cyan, running end-for-end the wrong way; the best of the "
             "12 modes in magenta, which does reproduce the crystal pose but is ranked "
             "last. Hinge residues Asp133 and Val135 labelled.", width_cm=9.5)
    R.p("The consequence is carried through the rest of the report: **this protocol is a "
        "ranking and triage filter, not a pose predictor.** The binding modes below are "
        "hypotheses about which residues a molecule can reach. What has to be validated "
        "instead is whether the ranking puts real actives near the top — Section 5.3.")

    R.h2("5.2  Docking run")
    R.p(f"AutoDock Vina 1.2.5, rigid receptor, exhaustiveness 32, 12 modes, energy range "
        f"4 kcal/mol, fixed seed 20241121, one CPU per ligand, 14 in parallel. "
        f"{enrich['n_compounds']} ligands were docked into the model (production) and the "
        "same set into the 1Q5K crystal receptor (cross-check); one config and one log "
        "file per ligand are in `05_docking/<receptor>/`.")

    R.h2("5.3  Does the ranking recover the known actives?")
    e1, e5, e10 = (enrich["EF1pct_model"], enrich["EF5pct_model"], enrich["EF10pct_model"])
    R.p(f"Scoring the {enrich['n_actives']} positive controls against the "
        f"{enrich['n_compounds'] - enrich['n_actives']} library compounds as background "
        f"gives **ROC-AUC {enrich['roc_auc_model']}** on the model "
        f"(BEDROCα=20 {enrich['bedroc20_model']}), with enrichment factors "
        f"EF5% = {e5['EF']} ({e5['actives_found']} of {enrich['n_actives']} actives in the "
        f"top {e5['n_top']}) and EF10% = {e10['EF']} ({e10['actives_found']} in the top "
        f"{e10['n_top']}); the single top-ranked compound is not an active, so EF1% is 0. "
        "**That is weak enrichment.** An AUC of 0.5 is a coin toss and 0.8+ is what a "
        "usable retrospective screen looks like. Four controls land in the top 16 "
        "(indirubin-3′-monoxime 4th, tideglusib 5th, SB-415286 8th, CHIR-99021 16th), "
        "but the two most potent compounds in the set — LY-2090314 at 1.1 nM and "
        "alsterpaullone at 4 nM — sit at ranks 71 and 68. The ranking is a coarse "
        "triage filter and is treated as one below.")
    if enrich.get("roc_auc_crystal") is not None:
        R.p(f"Running the identical protocol against the 1Q5K crystal receptor gives "
            f"ROC-AUC {enrich['roc_auc_crystal']}, so **the homology model is neither "
            "better nor worse as a virtual-screening receptor than the crystal it was "
            "built from** — which is the result that matters for Stage 2: the model "
            "is not the weak link. Across all compounds the two score sets correlate at "
            f"Pearson r = {enrich['model_vs_crystal_pearson']}, Spearman ρ = "
            f"{enrich['model_vs_crystal_spearman']}, with the model scoring on average "
            f"{abs(enrich['mean_score_shift_model_minus_crystal']):.2f} kcal/mol "
            f"{'weaker' if enrich['mean_score_shift_model_minus_crystal'] > 0 else 'stronger'}.")
    R.p(f"What the score does **not** do is track potency. Across the "
        f"{enrich['actives_n_with_potency']} controls with a measured value, spanning "
        "1.1 nM to 1.7 µM — three orders of magnitude — the correlation "
        f"between Vina score and pChEMBL is Pearson r = "
        f"{enrich['score_vs_pchembl_pearson']}, Spearman ρ = "
        f"{enrich['score_vs_pchembl_spearman']}. The screen can tell a plausible binder "
        "from an implausible one; it cannot rank binders by how good they are, and no "
        "statement below treats a 0.3 kcal/mol gap as meaningful.")
    R.p("Two caveats keep this honest: the library compounds are a **background set, not "
        "confirmed inactives**, so even this AUC is optimistic; and the controls are "
        "chemically biased towards flat heteroaromatics.")
    R.figure("fig4_score_distribution.png",
             "Score distribution of the library, with the known actives marked as vertical "
             "lines.", width_cm=10.5)
    R.figure_pair("fig5_roc.png",
                  "ROC curves for recovery of the 11 known actives, model versus crystal. "
                  "Both hug the diagonal: the enrichment is real but small.",
                  "fig6_model_vs_crystal.png",
                  "Per-compound scores on the two receptors: the model is a faithful "
                  "surrogate for the crystal.", width_cm=7.4)

    R.h2("5.4  Top of the ranking")
    rows = []
    for r in ranking[:12]:
        kind = {"known_active": "control", "A_approved_drug": "drug",
                "B_natural_product": "NP"}[r["set"]]
        rows.append([r["rank"], r["name"][:26], kind, r["best_score"],
                     r["ligand_efficiency"], r["cns_mpo4"], r["crystal_rank"],
                     r["mw"]])
    R.table(["#", "compound", "class", "Vina (kcal/mol)", "LE", "CNS-MPO (0-4)",
             "rank on crystal", "MW"], rows,
            caption="Top 12 of the screen on the homology model. LE = ligand efficiency "
                    "(−score / heavy atom); CNS-MPO is four of the six Wager "
                    "desirability terms (MW, clogP, TPSA, HBD), so 4 is the best possible. "
                    "Full ranking: 06_analysis/ranking_model.csv.",
            widths=[0.9, 4.6, 1.5, 2.3, 1.3, 2.1, 1.5, 1.4])

    # ----------------------------------------------------------- Stage 6 --
    R.h1("6  Interaction analysis and the three candidates")

    R.h2("6.1  Do the hits use the same residues as the actives?")
    R.p("For the 11 controls and the top 20 library hits the best pose was written into a "
        "complex with the prepared receptor, inspected visually and analysed with "
        "**PLIP**, giving a residue × compound interaction fingerprint. The poses "
        "were inspected and rendered in PyMOL; the dashed contacts in every figure are "
        "drawn from the PLIP geometry rather than from PyMOL’s own distance "
        "criteria, so figure and text cannot disagree. The 31 complexes are in "
        "`06_analysis/complexes/` and the .pml scripts in "
        "`06_analysis/figures/pymol/`.")
    R.p("The answer is the single most informative result in this report, and it is not "
        "the one we expected. The reference interaction is the **hinge**: in the 1Q5K "
        "crystal, AR-A014418 donates two hydrogen bonds to the Val135 backbone carbonyl "
        "(2.49 Å and 3.11 Å), and essentially every ATP-competitive GSK-3β "
        "inhibitor with a published co-crystal does something equivalent at Asp133/Val135. "
        f"In our docked poses, **not one of the {hinge['n_actives']} known actives makes "
        f"a hinge hydrogen bond, while {hinge['hits_hinge']} of the {hinge['n_hits']} top "
        "library hits do.**")
    R.p("Read alone that would say the library hits bind better than the known inhibitors, "
        "which is absurd. Read with Section 5.1 it is coherent: the protocol "
        "systematically mis-orients the chemotype it should handle best. The controls are "
        "flat hinge binders; Vina rotates them to bury hydrophobic surface and loses the "
        "hinge, the same failure the redocking test caught. What they are left making is "
        "hydrophobic contact with Ile62, Val70 and Leu188 plus water bridges to Asp200 "
        "— right cleft, wrong orientation. The polyphenols carry "
        "hydroxyls on every edge, so whichever way they are rotated something reaches the "
        "hinge. We therefore still require a hinge bond below, but as **prior structural "
        "knowledge from the crystallography**, not as a pattern learned from our own "
        "docked controls.")
    R.figure("fig8_fingerprint.png",
             "PLIP interaction fingerprint. Rows: known actives (red, above the line) then "
             "the top library hits, each labelled with its Vina score. Columns: ATP-site "
             "residues. Dark = hydrogen bond, water bridge or salt bridge; light = "
             "hydrophobic or π contact. Hinge residues are labelled in red.",
             width_cm=11.5)


    R.h2("6.2  Selection criteria")
    R.p("Four requirements were applied mechanically to every library compound "
        "(`06_analysis/shortlist.csv` records pass/fail on each for all 127): "
        "**(i)** a Vina score at least as good as the median positive control "
        "(−7.14 kcal/mol); **(ii)** at least one hydrogen bond to the hinge "
        "(Asp133/Val135), the anchor the crystal structures use; **(iii)** ligand "
        "efficiency ≥ 0.30 kcal/mol per heavy atom, so the score is not just "
        "molecular weight; **(iv)** a plausible route to the brain (CNS-MPO ≥ 3 of 4) "
        "and, for a repurposing library, an existing human safety record. Exactly **one** "
        "compound passes all four; the other two each miss one, and we say which rather "
        "than quietly moving the threshold.")
    R.p(cands.get("rejected", ""))

    R.h2("6.3  The candidates")
    for i, c in enumerate(cands["candidates"], 1):
        r = rank_by_id[c["id"]]
        R.h2(f"Candidate {i}: {c['name']}  —  rank {r['rank']}, "
             f"{float(r['best_score']):.2f} kcal/mol, LE {r['ligand_efficiency']}")
        R.p(f"*{c['what_it_is']}*")
        R.p("**Docking result.** " + c["docking"])
        R.p("**Key interactions (PLIP).** " + c["interactions"])
        R.p("**What it would do to the Assignment 1 mechanism.** " + c["mechanism"])
        R.p("**Against it.** " + c["against"])
        R.figure(c["figure"], c["figure_caption"], width_cm=8.0)

    # ----------------------------------------------------------- Stage 7 --
    R.h1("7  Limitations")
    for b in cands["limitations"]:
        R.bullet(b)

    R.h1("8  What we would do next")
    R.p(cands["next_steps"])

    R.h1("References")
    R.columns(2)
    for i, ref in enumerate(cands["references"], 1):
        p = R.p(f"{i}.  {ref}", size=7.8, space_after=1)
        p.paragraph_format.left_indent = Pt(18)
        p.paragraph_format.first_line_indent = Pt(-18)

    R.columns(1)
    R.h1("Appendix  —  software and data")
    R.table(["step", "tool / version", "key settings"],
            [["template search", "NCBI BLASTp (db pdb) + RCSB sequence search", "E < 1e-5"],
             ["comparative model", "MODELLER 10.7",
              "3 templates, AutoModel, 5 models, md_level = refine.slow, DOPE/GA341"],
             ["model validation", "MODELLER DOPE profile; Top8000 Ramachandran density; "
                                  "Biopython superposition", "15-residue smoothing"],
             ["pocket / prep", "P2Rank 2.5; PDB2PQR 3.x + PROPKA; Meeko",
              "default model; AMBER, pH 7.4; rigid receptor, Gasteiger"],
             ["ligand prep", "Open Babel (pH 7.4) → RDKit ETKDGv3 + MMFF94s → Meeko",
              "seed 0xF00D, 500 minimisation steps"],
             ["docking / analysis", "AutoDock Vina 1.2.5; PLIP; RDKit filter catalogs",
              f"box {int(box['size'][0])}³ Å, exhaustiveness 32, 12 modes, "
              f"seed 20241121; PAINS + Brenk"],
             ["figures", "PyMOL 3.2 (open source), matplotlib", "ray-traced, 1800 px; interaction geometry from PLIP"]],
            caption="Software and settings. Every script is in the data folder under "
                    "`scripts/`, numbered in run order.",
            widths=[3.2, 5.4, 7.4])
