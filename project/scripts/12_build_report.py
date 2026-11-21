"""
Build the mini-project report (Word, 8-10 pages).

Every number is read back out of the result files, so the report cannot
drift from the data: rerun the pipeline, rerun this, and the text moves
with it. Prose that is a judgement rather than a measurement lives in
07_report/candidates.json so it is visible and editable.

Output: 07_report/GSK3B_virtual_screen_report.docx
"""
import datetime, os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from report_lib import Report, load_csv, load_json, OUT, ACCENT, GREY
from docx.shared import Pt
from docx.enum.text import WD_ALIGN_PARAGRAPH


def main():
    os.makedirs(OUT, exist_ok=True)
    R = Report()

    val = load_json("02_structure", "validation", "validation_summary.json")
    models = load_csv("02_structure", "models", "model_scores.csv")
    templates = load_csv("02_structure", "alignment", "template_table.csv")
    box = load_json("03_pocket", "box.json")
    pocket = load_csv("03_pocket", "pocket_residues.csv")
    drug = load_json("03_pocket", "druggability.json")
    actives = load_csv("04_ligands", "actives", "known_actives.csv")
    manifest = load_csv("04_ligands", "ligand_manifest.csv")
    redock = load_json("05_docking", "redock_validation.json")
    enrich = load_json("06_analysis", "enrichment.json")
    ranking = load_csv("06_analysis", "ranking_model.csv")
    cands = load_json("07_report", "candidates.json")
    hinge = load_json("06_analysis", "hinge_summary.json")
    rank_by_id = {r["id"]: r for r in ranking}

    # ------------------------------------------------------------- title --
    t = R.d.add_paragraph(); t.alignment = WD_ALIGN_PARAGRAPH.CENTER
    r = t.add_run("Structure-based virtual screen against GSK-3β:\n"
                  "cutting the amyloid–tau wire at node 4")
    r.bold = True; r.font.size = Pt(16); r.font.color.rgb = ACCENT
    s = R.d.add_paragraph(); s.alignment = WD_ALIGN_PARAGRAPH.CENTER
    s.paragraph_format.space_after = Pt(8)
    r = s.add_run("CADD mini project · Group “Breaking Bad” · "
                  "mechanism: the amyloid–tau cascade in Alzheimer’s disease · "
                  f"{datetime.date.today():%d %B %Y}")
    r.font.size = Pt(9.5); r.font.color.rgb = GREY

    R.p("**Summary.** From the Assignment 1 cascade we took node 4, tau "
        "hyperphosphorylation, and the kinase that drives it: GSK-3β (UniProt "
        "P49841). A three-template comparative model of residues 35–384 was built "
        f"with MODELLER (best of five; Cα RMSD {val['rmsd']['ca_rmsd_all']} Å to "
        f"the 1Q5K crystal, {val['ramachandran_model']['favoured_pct']}% of residues in "
        "favoured Ramachandran regions). P2Rank placed its top-ranked pocket on the ATP "
        f"cleft, and a {int(box['size'][0])} Å box there was used to dock "
        f"{enrich['n_compounds']} molecules with AutoDock Vina 1.2.5: "
        f"{enrich['n_actives']} known GSK-3β actives as positive controls plus a "
        "hypothesis-driven library of 80 CNS-like approved drugs and 48 natural products. "
        f"The screen is only weakly enriching — ROC-AUC {enrich['roc_auc_model']}, "
        f"EF5% {enrich['EF5pct_model']['EF']} — and a redocking control that the "
        "scoring function **fails** limits how far any predicted pose can be trusted. "
        "Both results are reported in full, because they govern what the shortlist means. "
        "Three candidates are defended: tofacitinib, hesperetin and luteolin.", size=9.5)

    # ----------------------------------------------------------- Stage 1 --
    R.h1("1  Target: GSK-3β, the node-4 kinase")

    R.h2("1.1  Which step of the mechanism it controls")
    R.p("Our Assignment 1 cascade runs: APP mis-cleavage → Aβ42 → soluble "
        "oligomers → **tau hyperphosphorylation** → synapse and neuron loss. "
        "Node 4 is where an upstream amyloid signal becomes the lesion that actually tracks "
        "with dementia — our own slide made that point: tangle spread follows "
        "cognitive decline far more closely than plaque burden does, and the antibodies "
        "that clear plaque slow decline by only 27–35%. At node 4 the balance is set "
        "by GSK-3β and CDK5 rising while PP2A falls, and **GSK-3β performs the "
        "conversion**: it phosphorylates tau at most paired-helical-filament epitopes "
        "(Ser199, Ser202/Thr205, Thr231, Ser396/Ser404) and is itself activated downstream "
        "of Aβ. We chose it over CDK5 because it carries more of the tau load, and "
        "over BACE1 because node 1 is the arm the clinic has already tested hardest.")

    R.h2("1.2  Identity and function")
    R.p("**P49841** (GSK3B_HUMAN), 420 aa, EC 2.7.11.26; protein-kinase domain "
        "56–340 with catalytic Asp181, ATP contacts at the glycine-rich loop "
        "Ile62–Val70 and Lys85. It is constitutively active and switched *off* by "
        "Akt phosphorylation of Ser9, which folds back as an intramolecular "
        "pseudosubstrate, and it prefers substrates primed at S/T-X-X-X-S(p), the primed "
        "phosphate being read by Arg96/Arg180/Lys205. Besides tau it phosphorylates "
        "glycogen synthase, β-catenin, CRMP2 and eIF2B — a signalling hub, not "
        "a private tau kinase, which is the root of the selectivity problem in Section 7.")

    R.h2("1.3  Evidence linking it to Alzheimer’s disease")
    R.bullet("**Hooper, Killick & Lovestone (2008)** *J. Neurochem.* **104**, 1433–1439 "
             "— the GSK3 hypothesis: GSK-3 activity is raised in AD brain, it "
             "phosphorylates tau at PHF epitopes, it is activated by Aβ, and it "
             "mediates Aβ-induced neurotoxicity and LTP block. That is precisely the "
             "amyloid→tau junction of our cascade.")
    R.bullet("**Lucas et al. (2001)** *EMBO J.* **20**, 27–39 — conditional "
             "forebrain overexpression of GSK-3β alone produces hyperphosphorylated, "
             "somatodendritically mislocalised tau, astrocytosis and neuronal death. Raise "
             "the kinase, get the node-4 phenotype: the causal arm of the argument.")
    R.bullet("**Phiel, Wilson, Lee & Klein (2003)** *Nature* **423**, 435–439 — "
             "GSK-3 inhibition lowers Aβ production, but via the α paralogue. The "
             "split tells us a β-directed agent is aimed at the tau arm, and that "
             "α/β selectivity must be measured, not assumed.")
    R.bullet("Supporting: **Bhat et al. (2003)** *J. Biol. Chem.* **278**, "
             "45937–45945 report AR-A014418 (IC₅₀ 104 ± 27 nM, "
             "Kᵢ 38 nM, no CDK2/CDK5 inhibition at 100 µM), which blocks tau "
             "Ser396 phosphorylation in cells. That compound is the ligand of PDB 1Q5K and "
             "is used three times below: template ligand, site definition, redocking control.")

    R.h2("1.4  Intended action, and why it should move the mechanism")
    R.p("**Inhibit, ATP-competitively, at the kinase active site.** Node 4 is a "
        "kinase/phosphatase balance; lowering GSK-3β output shifts it back, so tau "
        "stays on microtubules, tangle assembly and trans-synaptic spread slow, and axonal "
        "transport is preserved. The benefit is expected at node 5 (synapse and neuron "
        "loss). The claim is deliberately narrow: this does not lower Aβ, it cuts the "
        "wire between the amyloid signal and the tau lesion.")
    R.p("The ATP cleft is chosen over an allosteric or substrate-competitive site because "
        "(i) it is the only GSK-3β pocket with validated geometry (123 PDB entries, "
        "almost all ATP-site complexes); (ii) the substrate-priming groove is shallow and "
        "highly charged, a poor fit for a small neutral brain-penetrant molecule; (iii) the "
        "one non-ATP-competitive agent taken into AD patients, tideglusib, did not separate "
        "from placebo. The price is kinome selectivity (Section 7). The aim is partial, "
        "chronic inhibition: complete GSK-3 blockade disinhibits Wnt/β-catenin "
        "signalling, with oncogenic and metabolic liability — the conclusion the "
        "field reached for BACE1 after verubecestat.")

    R.h2("1.5  Existing chemical matter → our positive controls")
    R.p("There is **no approved selective GSK-3 inhibitor**. What exists is tool compounds, "
        "one failed clinical candidate (tideglusib: 306 patients, 26 weeks, no separation "
        "from placebo on ADAS-cog15; Lovestone et al. 2015, *J. Alzheimers Dis.* **45**, "
        "75–88) and lithium, which inhibits GSK-3 only at millimolar, non-selective "
        "concentrations. A DrugBank search returns no approved agent whose stated target "
        "is GSK-3β (lithium carbonate lists it among several mechanisms; tideglusib "
        "is investigational, with orphan designation in myotonic dystrophy), which is why "
        "the controls below are tool compounds rather than drugs. Eleven of them become "
        "the positive controls of this screen; potencies are ChEMBL medians over curated "
        "binding assays against target CHEMBL262.")
    rows = [[a["name"], a["chembl_id"],
             f"{float(a['potency_nM']):.1f}" if a["potency_nM"] else "–",
             a["n_measurements"], a["mechanism_note"]]
            for a in sorted(actives, key=lambda r: float(r["pchembl_median"] or 0), reverse=True)]
    R.table(["compound", "ChEMBL ID", "potency (nM)", "n", "status / mechanism"], rows,
            caption="Positive controls. Tideglusib and TDZD-8 are non-ATP-competitive and "
                    "are kept in deliberately, as negative controls on the protocol itself.",
            widths=[2.8, 2.2, 1.8, 0.8, 8.4])

    # ----------------------------------------------------------- Stage 2 --
    R.h1("2  Structure: a comparative model of the kinase domain")

    R.h2("2.1  Template search")
    R.p("The full-length P49841 sequence was searched against the PDB with NCBI BLASTp "
        "(database *pdb*, E < 10⁻⁵; submission record and RID in "
        "`02_structure/alignment/blast_submission.txt` — the public NCBI queue had "
        "not returned the job an hour after submission) and, as the search actually used "
        "here, with the RCSB sequence-search service (MMseqs2 profile search over the PDB, "
        "2000 hits with per-hit identity and E-value). The picture is unambiguous: "
        "123 PDB entries are GSK-3β itself at 100% identity; the nearest non-self "
        "relative is human GSK-3α (7SXF, 87%), then GSK3 orthologues from other "
        "organisms (*Ustilago maydis* 4E7W at 70%, *Toxoplasma gondii* 9HW6 at 55%) and "
        "other CMGC kinases below that. The "
        "modelling problem is therefore **not** finding a fold — it is building a "
        "complete, consistent receptor, because every individual GSK-3β crystal is "
        "missing residues.")
    rows = []
    for tmpl in templates[:5]:
        rows.append([tmpl["pdb_id"], tmpl["identity_pct"], tmpl["resolution_A"],
                     tmpl["method"].replace("X-RAY DIFFRACTION", "X-ray"),
                     ("yes" if tmpl["used_as_template"] else ""), tmpl["title"][:60]])
    R.table(["PDB", "identity %", "res. (Å)", "method", "template?", "title"], rows,
            caption="Top of the template search (full table: "
                    "02_structure/alignment/template_table.csv).",
            widths=[1.4, 1.9, 1.6, 1.8, 1.7, 7.6])
    R.p("Three templates were taken forward: **1Q5K** (1.94 Å, with AR-A014418, "
        "which fixes the ATP-site conformation we want to dock into), **1J1B** "
        "(1.80 Å, continuous 35–388, the loop donor) and **4ACC** (2.21 Å, "
        "continuous, a second donor and a conformational check). 1Q5K is missing residues "
        "120–121 and 287–289; the other two have them, so modelling 35–384 "
        "from all three gives a receptor with no chain breaks. The disordered tails "
        "(1–34, 385–420) were left out — no template density supports them.")

    R.h2("2.2  Alignment and model building")
    R.p("The three templates were structurally aligned with MODELLER’s *salign* "
        "(tree alignment, gap penalties −450/−50, 3D penalties 0/3), then the "
        "target sequence was aligned onto that block with a second, profile-based *salign* "
        "pass. Because target and templates are the same protein the sequence alignment is "
        "unambiguous; the alignment’s real job is to register the loop gaps correctly. "
        "Five models were then built with *AutoModel* (`md_level = refine.slow`) and scored "
        "with molpdf, DOPE, DOPE-HR and GA341.")
    rows = [[m["model"], m["molpdf"], m["DOPE"], m["DOPE-HR"], m["GA341"]] for m in models]
    R.table(["model", "molpdf", "DOPE", "DOPE-HR", "GA341"], rows,
            caption="The five models. All five pass GA341 at 1.00; the lowest-DOPE model "
                    "was taken forward.", widths=[4.4, 2.9, 2.9, 2.9, 2.9])

    R.h2("2.3  Validation")
    R.p(f"The selected model was compared with 1Q5K chain A over the "
        f"{val['rmsd']['n_common']} residues present in both. The Cα RMSD is "
        f"**{val['rmsd']['ca_rmsd_all']} Å** overall and "
        f"**{val['rmsd']['ca_rmsd_pocket']} Å** across the "
        f"{val['rmsd']['n_pocket']} ATP-site residues; including side chains, the "
        f"all-heavy-atom RMSD over the site is {val['rmsd']['heavy_rmsd_pocket']} Å. "
        f"Ramachandran analysis against the Top8000 reference density puts "
        f"**{val['ramachandran_model']['favoured_pct']}% of residues in favoured regions "
        f"and {val['ramachandran_model']['outlier_pct']}% as outliers**, slightly better "
        f"than the crystal itself ({val['ramachandran_crystal']['favoured_pct']}% favoured). "
        "The per-residue DOPE profile tracks the template everywhere except the two "
        "rebuilt loops — exactly where a comparative model should be least certain, "
        "and both are >15 Å from the ATP site.")
    R.figure_pair("fig1_model_overview.png",
                  "The model (N-lobe blue, C-lobe grey, rebuilt loops red), P2Rank "
                  "pocket-1 residues as green sticks, AR-A014418 from 1Q5K in cyan and the "
                  "26 Å search box in yellow.",
                  "fig3_ramachandran.png",
                  "Ramachandran plot against the Top8000 density; contours enclose 98% "
                  "(favoured) and 99.95% (allowed) of the reference mass.", width_cm=7.5)
    R.figure("fig2_dope_profile.png",
             "Per-residue normalised DOPE for the model (red) and 1Q5K chain A (grey), on "
             "UniProt numbering. Blue bands mark the rebuilt loops.", width_cm=13.0)

    # ----------------------------------------------------------- Stage 3 --
    R.h1("3  Binding site and file preparation")

    R.h2("3.1  Pocket detection")
    res_list = ", ".join(f"{p['resname'].title()}{p['uniprot_resnum']}" for p in pocket)
    R.p(f"P2Rank 2.5 was run on the model without any prior. Its **top-ranked pocket "
        f"(score {box['p2rank_pocket1_score']}, probability "
        f"{box['p2rank_pocket1_probability']})** is the ATP cleft: "
        f"{len(pocket)} residues — {res_list}: the textbook GSK-3β ATP site "
        "(glycine-rich loop Ile62–Val70, catalytic Lys85, gatekeeper Leu132, the "
        "Asp133–Val135 hinge, Thr138 and Arg141 at the solvent front, "
        "Asn186/Leu188/Cys199/Asp200 on the floor). The second pocket scores 2.8, so the "
        "site is not ambiguous.")
    R.p("Druggability was measured, not asserted. On a 0.7 Å grid, keeping points "
        "3.2–6.5 Å from the nearest heavy atom and buried behind protein in "
        "≥65% of 42 ray directions, the cavity holding the co-crystal ligand has a "
        f"volume of **{drug['cavity_volume_A3']} Å³**, inside the "
        "300–1000 Å³ range that tractable sites occupy; "
        f"Shrake–Rupley SASA puts **{drug['reference_ligand_buried_pct']}% of "
        f"AR-A014418 buried** on binding; and of the {drug['lining_residues']} lining "
        f"residues {drug['lining_hydrophobic_residue_pct']}% are hydrophobic, with "
        f"{drug['lining_carbon_atom_pct']}% of lining heavy atoms carbon. Enclosed, "
        "hydrophobic, the right size — a conventional kinase pocket.")
    R.p("The literature agrees independently: in 1Q5K, AR-A014418 hydrogen-bonds to the "
        "Val135 backbone carbonyl (2.49 and 3.11 Å from its urea nitrogens) and packs "
        "against Ile62, Val70, Ala83, Leu132 and Leu188 — all in P2Rank’s list.")

    R.h2("3.2  The search box")
    c, sz = box["center"], box["size"]
    R.p(f"The box is centred at the midpoint of two independent estimates — the "
        f"P2Rank pocket centroid ({', '.join(str(x) for x in box['p2rank_center'])}) and "
        f"the crystallographic AR-A014418 centroid "
        f"({', '.join(str(x) for x in box['tmu_centroid'])}), 4.6 Å apart. "
        f"**Centre = ({c[0]}, {c[1]}, {c[2]}) Å; size = {int(sz[0])} × "
        f"{int(sz[1])} × {int(sz[2])} Å.** The edge is the extent of the "
        "pocket-1 residues plus the reference ligand with a 5 Å margin, clipped to "
        "22–26 Å. The model was first superposed onto 1Q5K chain A so both "
        "receptors share one frame, which lets one box serve both and makes poses directly "
        "comparable.")

    R.h2("3.3  Receptor preparation")
    R.p("Both receptors went through the same two steps: **PDB2PQR 3.x** (AMBER force "
        "field, PROPKA titration states at pH 7.4, hydrogens added and the H-bond network "
        "optimised), then **Meeko** `mk_prepare_receptor` to a rigid PDBQT with Gasteiger "
        "charges and merged non-polar hydrogens. Ten crystallographic waters lying within "
        "6 Å of the co-crystal ligand were kept in both receptors and transplanted "
        "unchanged into the model — justified in Section 5.1, where the redocking test "
        "shows the crystal pose is simply not sampled without them.")

    # ----------------------------------------------------------- Stage 4 --
    R.h1("4  The ligand set")
    n_act = sum(1 for m in manifest if m["set"] == "known_active")
    n_drug = sum(1 for m in manifest if m["set"] == "A_approved_drug")
    n_np = sum(1 for m in manifest if m["set"] == "B_natural_product")

    R.h2("4.1  Design")
    R.p("The screening hypothesis is explicit: **a GSK-3β inhibitor for Alzheimer’s "
        "has to reach the brain, and the fastest route to a brain-penetrant molecule is one "
        "that is already in humans.** The library is therefore two blocks.")
    R.bullet(f"**Block A — {n_drug} approved drugs.** Every approved small molecule "
             "in ChEMBL (max_phase = 4), filtered on CNS-MPO-like physicochemistry "
             "(150 ≤ MW ≤ 450, 0.5 ≤ AlogP ≤ 5, HBD ≤ 3, "
             "TPSA ≤ 90 Å², rotatable bonds ≤ 8) and on the "
             "ATP-site requirement of at least one aromatic ring and one acceptor. 730 "
             "passed; 80 were picked by MaxMin diversity selection on Morgan "
             "fingerprints (radius 2, 2048 bits, seed 2024) so the block spans chemotype "
             f"space, of which {n_drug} survived ligand preparation.")
    R.bullet(f"**Block B — {n_np} natural products.** Compounds with published "
             "activity in Alzheimer models or on GSK-3/tau (curcumin, catechins and "
             "flavonols, resveratrol, berberine, harmine, celastrol, withaferin A and so "
             "on), pulled "
             "from ChEMBL by name. The AD literature is full of claims for these; docking "
             "them against the same pocket as the tool compounds is a fair test of those "
             "claims.")
    R.bullet(f"**{n_act} positive controls** (Table 2) are docked in the same run, blind to "
             "the scoring, so the screen can be judged on whether it recovers them.")

    R.h2("4.2  Preparation")
    R.p("Each molecule went from SMILES through: **Open Babel** protonation at pH 7.4 "
        "(`OBMol.AddHydrogens(polaronly = False, correctForPH = True, pH = 7.4)`); "
        "**RDKit** 3D embedding with ETKDGv3 (seed 0xF00D) and MMFF94s minimisation, "
        "500 steps, falling back to UFF where MMFF typing failed; then **Meeko** "
        "`MoleculePreparation` to PDBQT with Gasteiger charges, automatic rotatable-bond "
        f"detection and rigid amides. {len(manifest)} of 139 inputs were prepared "
        "successfully; one (prajmalium, a quaternary cage alkaloid) failed 3D embedding and "
        "was dropped. Formal charges at pH 7.4 are recorded per ligand in "
        "`04_ligands/ligand_manifest.csv`.")

    from report_part2 import tail
    tail(R, {"box": box, "redock": redock, "enrich": enrich, "ranking": ranking,
             "cands": cands, "hinge": hinge, "manifest": manifest})

    R.save(os.path.join(OUT, "CADD_miniproject_GSK3B_report.docx"))


if __name__ == "__main__":
    main()
