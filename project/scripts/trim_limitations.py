"""Shorten the limitations and next-steps text in 07_report/candidates.json.

Same content, fewer words - the report has a page budget.
"""
import json, os

P = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "07_report",
                 "candidates.json")

LIMITATIONS = [
 "**The redocking control failed, and that sets the ceiling on everything else.** Vina "
 "ranks the reversed orientation of AR-A014418 first; the crystal pose, which the search "
 "does find once the conserved waters are present, comes twelfth. Scoring the crystal pose "
 "directly (−6.58 kcal/mol) against the reversed one (−7.67) locates the error in "
 "the scoring function, not the search. Every pose here is a hypothesis about reachable "
 "residues, not a claim about how the molecule sits.",

 "**Score is not affinity, and the enrichment is weak.** Across eleven controls spanning "
 "three orders of magnitude of potency, the score barely tracks pChEMBL, and the ROC-AUC "
 "of 0.65 is only modestly better than chance. Differences of a few tenths of a kcal/mol "
 "between neighbouring compounds carry no information.",

 "**The background is not a set of inactives.** Enrichment was measured against 127 "
 "library compounds that have simply never been tested on GSK-3β, and no "
 "property-matched decoy set (DUD-E, DEKOIS) was used, so even that modest AUC is an "
 "optimistic estimate.",

 "**The receptor is rigid and inhibitor-shaped.** The glycine-rich loop and the "
 "Lys85–Glu97 pair move substantially between GSK-3β structures. A model built "
 "on three inhibitor complexes is pre-organised for inhibitor-like ligands and biased "
 "against anything needing an induced fit; no ensemble and no flexible side chains.",

 "**The model is comparative, not predictive.** Target and templates are the same protein, "
 "so MODELLER completed two disordered loops rather than crossing a fold gap. The 0.43 "
 "Å RMSD says the rebuild is clean, not that homology modelling solved anything hard.",

 "**Solvation, tautomers, metals.** One microspecies per ligand at pH 7.4, no tautomer "
 "enumeration, no Mg²⁺ (present in physiological Mg-ATP), and only ten waters, "
 "transplanted from a structure solved with a different ligand bound.",

 "**Selectivity was not assessed at all.** GSK-3α shares essentially the whole ATP "
 "site and the CMGC kinases are close behind. Nothing here separates a GSK-3β "
 "inhibitor from a promiscuous kinase binder — one of our own controls, "
 "staurosporine, is exactly that.",

 "**Two controls cannot be scored by this protocol by construction.** Tideglusib and "
 "TDZD-8 are non-ATP-competitive, so docking them into the ATP site is meaningless; their "
 "ranks (5th and 113th) should be read that way, and tideglusib's high rank is a warning "
 "about reading any single rank too hard.",

 "**Repurposing and natural-product libraries bring their own traps.** Eighteen of the top "
 "forty compounds carry PAINS or Brenk alerts — mostly catechols and Michael "
 "acceptors — and such compounds are frequent hitters in real assays "
 "(`06_analysis/structural_alerts.csv`).",

 "**Nothing here is experimental, and brain exposure is only estimated.** No MD, no "
 "free-energy calculation, no assay. The CNS-MPO-like score uses four of the six Wager "
 "terms and says nothing about P-glycoprotein efflux, which is what actually kills most "
 "CNS candidates.",
]

NEXT = ("The honest next step is not more docking. First, fix the known defect: rescore the "
        "retained poses with a function that is not charge-blind — MM-GBSA on short "
        "restrained minimisations, or a trained scoring model — and check whether "
        "AR-A014418’s crystal orientation then ranks first. A protocol that cannot "
        "pass its own positive control should not prioritise purchases. Second, re-run "
        "against an ensemble of three or four GSK-3β conformations and a "
        "property-matched decoy set, so the enrichment number means something. Third, "
        "counter-screen the three candidates against GSK-3α, CDK2, CDK5 and DYRK1A, "
        "because selectivity is the real liability of an ATP-site programme. Only then is "
        "a kinase assay on the three candidates — all commercially available — "
        "worth bench time, followed by a cellular tau-Ser396 phosphorylation readout, "
        "which is the measurement that connects back to node 4 of the cascade.")


def main():
    d = json.load(open(P, encoding="utf8"))
    d["limitations"] = LIMITATIONS
    d["next_steps"] = NEXT
    json.dump(d, open(P, "w", encoding="utf8"), indent=1, ensure_ascii=False)
    print("limitations:", sum(len(x.split()) for x in LIMITATIONS), "words;",
          "next steps:", len(NEXT.split()), "words")


if __name__ == "__main__":
    main()
