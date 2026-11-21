"""Fill in the three candidate write-ups in 07_report/candidates.json.

The numbers in these paragraphs are read off 06_analysis/ranking_model.csv,
06_analysis/plip/interactions.csv and 06_analysis/shortlist.csv; the
judgements are ours and live here so they are visible and editable rather
than buried in the report builder.
"""
import json, os

HERE = os.path.dirname(os.path.abspath(__file__))
P = os.path.join(HERE, "..", "07_report", "candidates.json")

CANDIDATES = [
 {
  "id": "LIG068",
  "name": "Tofacitinib",
  "figure": "fig10_tofacitinib.png",
  "what_it_is": "An approved JAK1/JAK3 inhibitor (tofacitinib citrate): an ATP-competitive "
      "pyrrolo[2,3-d]pyrimidine, MW 312, clogP 1.5, TPSA 89 Å². It is the only "
      "compound in the library that passes all four selection criteria.",
  "docking": "Rank 9 of 138 at −8.48 kcal/mol, ligand efficiency 0.369 and CNS-MPO 3.8 "
      "of 4 — the best combination of score, efficiency and brain-likeness in the "
      "screen. It scores better than seven of the eleven positive controls. It is also the "
      "least surprising hit in the set: it is a kinase inhibitor, so an ATP-site pose is "
      "what it was designed for, and that is a point in its favour rather than against it.",
  "interactions": "A single, almost ideal hinge hydrogen bond — the Val135 backbone NH "
      "donates to the ligand at 2.80 Å with a 169° angle, the same anchor "
      "AR-A014418 uses in the 1Q5K crystal. Two further hydrogen bonds reach Asp200 of the "
      "DFG motif (3.14 Å ligand-donated, 3.61 Å protein-donated), and the "
      "scaffold is sandwiched between Ile62 and Leu188 with additional contacts to Lys85 "
      "and Tyr134. This is the canonical kinase-inhibitor pattern, reproduced without any "
      "of the reversal seen in the redocking control.",
  "mechanism": "Two nodes of the Assignment 1 cascade at once. At **node 4**, inhibiting "
      "GSK-3β lowers tau phosphorylation at the PHF epitopes, so tau stays bound to "
      "microtubules and tangle formation and trans-synaptic spread slow. At **node 6**, "
      "which our own slide listed as a drug node, JAK/STAT blockade damps microglial "
      "priming and the chronic cytokine loop that strips synapses. A molecule that touched "
      "both the kinase driving tangle formation and the inflammatory amplifier would act "
      "exactly at the junction our mechanism says matters.",
  "against": "Brain exposure is the problem: tofacitinib is a P-glycoprotein substrate with "
      "low CNS penetration, so the CNS-MPO score flatters it. Chronic systemic JAK "
      "inhibition in elderly patients carries boxed warnings for serious infection, "
      "thrombosis and malignancy, which makes straight repurposing implausible. And its "
      "GSK-3β activity here is predicted, not measured — published kinase panels "
      "show tofacitinib is quite JAK-selective, so the prior probability of potent "
      "GSK-3β inhibition is low. The honest reading is that the pyrrolopyrimidine "
      "hinge binder is a validated starting scaffold, not a drug to repurpose as it stands.",
  "figure_caption": "Tofacitinib (cyan) in the ATP site of the model. Blue dashes are "
      "hydrogen bonds, grey dashes hydrophobic contacts; small cyan spheres are the "
      "retained crystallographic waters. The Val135 hinge hydrogen bond "
      "(2.80 Å, 169°) is the anchor.",
 },
 {
  "id": "LIG093",
  "name": "Hesperetin",
  "figure": "fig11_hesperetin.png",
  "what_it_is": "A citrus flavanone, the aglycone of hesperidin, with a sizeable literature "
      "in Alzheimer models. MW 302, clogP 2.5, TPSA 96 Å². Alone among the nine "
      "flavonoids that reached the top 20 it carries no PAINS or Brenk structural alert "
      "— its B ring is 3′-hydroxy-4′-methoxy, not a catechol.",
  "docking": "Rank 11 at −8.38 kcal/mol with ligand efficiency 0.381, better than every "
      "positive control except indirubin-3′-monoxime. Its CNS-MPO of 2.96 misses our "
      "3.0 cut-off by a rounding error; we kept it because every flavonoid that clears the "
      "score cut-off instead fails the structural-alert check.",
  "interactions": "The second most complete hinge engagement in the screen: a bidentate "
      "hinge contact with Asp133 (2.85 Å) and Val135 (2.88 Å), both "
      "ligand-donated, plus hydrogen bonds to Asn186 (2.76 Å) and Asp200 (2.95 Å) "
      "and a water-mediated bridge to Asp200 through one of the retained crystallographic "
      "waters. Hydrophobic contacts to Ile62, Val70, Tyr134 and Leu188 complete the "
      "sandwich. This chromone-on-the-hinge geometry is the mode flavonoids adopt in "
      "published kinase co-crystals, which is weak independent support for the pose.",
  "mechanism": "**Node 4**, partially. A modest, chronic reduction in GSK-3β output is "
      "exactly the partial inhibition the target argues for, and a dietary flavanone is the "
      "one chemical class that could plausibly be dosed for years. Its reported antioxidant "
      "activity would also touch node 3, where metal-catalysed oligomer formation sits. We "
      "deliberately claim no more: the realistic potency of a flavanone on GSK-3β is "
      "micromolar, which would blunt the node rather than block it.",
  "against": "Potency is the weak point — flavanones generally inhibit kinases in the "
      "micromolar range, and this screen cannot tell micromolar from nanomolar. Hesperetin "
      "is heavily glucuronidated on first pass, so free brain concentrations after oral "
      "dosing are very low, and flavonoids as a class bind many kinases, so selectivity is "
      "unlikely. Its TPSA and three donors put it at the edge of CNS chemical space.",
  "figure_caption": "Hesperetin (cyan) in the ATP site. The chromone makes a bidentate hinge "
      "contact with Asp133 and Val135, and a retained water bridges the ligand to Asp200.",
 },
 {
  "id": "LIG084",
  "name": "Luteolin",
  "figure": "fig12_luteolin.png",
  "what_it_is": "A dietary flavone, and the compound the Alzheimer natural-product "
      "literature claims most often. MW 286, clogP 2.3. It has the highest ligand "
      "efficiency of any hinge-binding compound in the screen — and it carries a "
      "catechol, which both the PAINS and the Brenk filters flag.",
  "docking": "Rank 15 at −8.30 kcal/mol, ligand efficiency 0.395: the best of any "
      "compound that also makes a hinge hydrogen bond, and better than ten of the eleven "
      "controls. Twelve PLIP interactions, seven of them hydrogen bonds — the densest "
      "network in the top 20.",
  "interactions": "Bidentate hinge again — Asp133 at 2.81 Å and Val135 at "
      "2.83 Å (ligand-donated) and 3.10 Å (protein-donated) — plus Asn186 "
      "(2.72 Å, near-linear at 170°), Asp200 (2.96 Å) and Asn64 on the "
      "glycine-rich loop (2.81 Å). Hydrophobic contacts run to Ile62, Val70, Ala83 and "
      "Leu188. Every one of those residues is in the P2Rank pocket and in the 1Q5K contact "
      "list, so whatever else is true, the molecule fills the cleft the way ATP does.",
  "mechanism": "**Node 4** by the same logic as hesperetin, with the addition that "
      "luteolin’s claimed anti-inflammatory activity would reach node 6 (microglial "
      "priming) and its metal-chelating catechol would touch node 3 (metal-catalysed "
      "aggregation). It is the candidate whose claimed polypharmacology best matches the "
      "shape of our cascade — and, for exactly that reason, the one whose claims most "
      "need a selective assay rather than another disease model.",
  "against": "The catechol is the problem. Catechols are redox-active, oxidise to quinones, "
      "chelate metals and turn up as hits in a large fraction of screening campaigns "
      "regardless of target; a docking score on a catechol is weak evidence on its own. "
      "Luteolin also has poor oral bioavailability and is rapidly conjugated. We include it "
      "because it is the literature’s favourite and because the screen turns a vague "
      "claim into a concrete, falsifiable pose — but of the three it is the one most "
      "likely to be an artefact, and it needs a redox-interference counter-screen before "
      "anything else.",
  "figure_caption": "Luteolin (cyan) in the ATP site, with the densest hydrogen-bond network "
      "of the top 20: bidentate hinge (Asp133, Val135) plus Asn186, Asp200 and Asn64 on the "
      "glycine-rich loop.",
 },
]

REJECTED = (
    "Three compounds that scored better were rejected, and the reasons matter as much as "
    "the selections. **Withaferin A** (rank 1, −9.30) has a ligand efficiency of only "
    "0.273 — it buys its score with 34 heavy atoms — and carries a reactive "
    "epoxide and an enone; a withanolide that alkylates cysteines is not a GSK-3β "
    "inhibitor hypothesis, it is a reactivity hypothesis. **Niraparib** (rank 2, "
    "−9.22) and **capmatinib** (rank 3, −9.20) are both excellent-scoring "
    "ATP-site-shaped oncology drugs, but neither makes a hinge hydrogen bond in its pose, "
    "which is the one structural feature we are prepared to trust; both are also poor CNS "
    "agents with toxicity profiles no chronic neurodegeneration indication could absorb."
)


def main():
    d = json.load(open(P, encoding="utf8"))
    d["candidates"] = CANDIDATES
    d["rejected"] = REJECTED
    json.dump(d, open(P, "w", encoding="utf8"), indent=1, ensure_ascii=False)
    print("wrote", len(CANDIDATES), "candidates ->", os.path.abspath(P))


if __name__ == "__main__":
    main()
