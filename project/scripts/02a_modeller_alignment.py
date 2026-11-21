"""
Stage 2a - multiple-template alignment for the GSK-3beta comparative model.

Templates (chosen from the BLAST/PDB search, see 02_structure/alignment/):
    1q5k chain A   1.94 A   GSK-3beta + AR-A014418 (defines the ATP-site
                            conformation we want to dock into)
    1j1b chain A   1.80 A   GSK-3beta + AMP-PNP, continuous 35-388, used to
                            rebuild the loops that are disordered in 1q5k
                            (120-121 and 287-289)
    4acc chain A   2.21 A   GSK-3beta + inhibitor, continuous 35-385, second
                            loop donor / conformational check

Target: UniProt P49841 residues 35-384 (kinase domain 56-340 plus the ordered
flanks). The disordered N-terminal 1-34 and C-terminal 385-420 tails are not
modelled - there is no template density for them and modelling them ab initio
would only add noise around the binding site.

Run:  python 02a_modeller_alignment.py
Produces: 02_structure/alignment/gsk3b-mult.ali , .pap
"""
import os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import modeller_env  # noqa: F401  (sets up MODELLER paths)
import os, sys

HERE = os.path.dirname(os.path.abspath(__file__))
ALN = os.path.join(HERE, "..", "02_structure", "alignment")
TPL = os.path.join(HERE, "..", "02_structure", "templates")
SEQFILE = os.path.join(HERE, "..", "01_target", "P49841.seq")

TEMPLATES = [("1q5k", "A"), ("1j1b", "A"), ("4acc", "A")]
START, END = 35, 384


def write_target_pir():
    seq = open(SEQFILE).read().strip()
    frag = seq[START - 1:END]              # 1-based UniProt numbering
    lines = [">P1;GSK3B",
             f"sequence:GSK3B:{START}: :{END}: ::::"]
    for i in range(0, len(frag), 60):
        lines.append(frag[i:i + 60])
    lines[-1] += "*"
    path = os.path.join(ALN, "GSK3B.ali")
    open(path, "w").write("\n".join(lines) + "\n")
    print(f"target sequence {START}-{END} ({len(frag)} aa) -> {path}")
    return path


def main():
    from modeller import Environ, log, alignment, model
    log.verbose()
    env = Environ()
    env.io.atom_files_directory = [os.path.abspath(TPL), "."]
    os.chdir(ALN)

    # ---- 1. structural (template-template) alignment with SALIGN -----------
    aln = alignment(env)
    for code, chain in TEMPLATES:
        mdl = model(env, file=os.path.join(TPL, code.upper() + ".pdb"),
                    model_segment=("FIRST:" + chain, "LAST:" + chain))
        aln.append_model(mdl, atom_files=code.upper() + ".pdb",
                         align_codes=code + chain)
    aln.salign(rms_cutoff=3.0, normalize_pp_scores=False,
               rr_file="$(LIB)/as1.sim.mat", overhang=30,
               gap_penalties_1d=(-450, -50), gap_penalties_3d=(0, 3),
               gap_gap_score=0, gap_residue_score=0,
               alignment_type="tree", output="ALIGNMENT QUALITY")
    aln.write(file=os.path.join(ALN, "templates.ali"), alignment_format="PIR")
    aln.write(file=os.path.join(ALN, "templates.pap"), alignment_format="PAP",
              alignment_features="INDICES HELIX BETA STRAIGHTNESS ACCESSIBILITY CONSERVATION")
    aln.id_table(matrix_file=os.path.join(ALN, "template_identities.mat"))
    env.dendrogram(matrix_file=os.path.join(ALN, "template_identities.mat"),
                   cluster_cut=-1.0)

    # ---- 2. align the target onto the template block -----------------------
    write_target_pir()
    aln2 = alignment(env)
    aln2.append(file=os.path.join(ALN, "templates.ali"), align_codes="all")
    aln2.append(file=os.path.join(ALN, "GSK3B.ali"), align_codes="GSK3B")
    aln2.salign(output="", max_gap_length=20, gap_function=True,
                alignment_type="PAIRWISE", align_block=len(TEMPLATES),
                feature_weights=(1., 0., 0., 0., 0., 0.), overhang=0,
                gap_penalties_1d=(-450, 0), gap_penalties_2d=(0.35, 1.2, 0.9, 1.2, 0.6, 8.6, 1.2, 0., 0.),
                similarity_flag=True)
    aln2.write(file=os.path.join(ALN, "gsk3b-mult.ali"), alignment_format="PIR")
    aln2.write(file=os.path.join(ALN, "gsk3b-mult.pap"), alignment_format="PAP")
    print("\nalignment written to", os.path.join(ALN, "gsk3b-mult.ali"))


if __name__ == "__main__":
    main()
