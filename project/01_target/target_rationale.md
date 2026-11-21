# Stage 1 - Target selection and justification

## 1.1 The protein and the step it controls

**Target: glycogen synthase kinase-3 beta (GSK-3β)**
UniProt **P49841** (GSK3B_HUMAN), EC 2.7.11.26, 420 aa, protein-kinase domain
residues 56-340, catalytic aspartate (proton acceptor) Asp181, ATP-binding
residues Ile62-Val70 and Lys85.

In the Assignment 1 cascade GSK-3β sits at **node 4, "tau
hyperphosphorylation"**. The cascade runs: APP mis-cleavage → Aβ42 →
soluble oligomers → **tau hyperphosphorylation** → synapse and neuron loss.
Node 4 is the step at which an upstream amyloid signal is converted into the
cytoskeletal lesion that actually tracks with dementia. GSK-3β is the kinase
that performs that conversion: it phosphorylates tau at the majority of the
paired-helical-filament epitopes (Ser199, Ser202/Thr205 (AT8), Thr231,
Ser396/Ser404 (PHF-1)), and the balance at node 4 is set by GSK-3β and CDK5
rising while PP2A falls.

Choosing node 4 rather than node 1 or 2 is deliberate. Our own Assignment 1
slide made the point that **tangle spread tracks cognitive decline far more
closely than plaque burden does** (Braak staging), and that the anti-amyloid
antibodies that clear plaque only slow decline by 27-35%. If the readout that
matters is the tau arm, the kinase that drives the tau arm is the rational
place to intervene.

## 1.2 Function

GSK-3β is a constitutively active serine/threonine kinase, unusual in that it
is switched *off* rather than on by signalling: insulin/PI3K-Akt
phosphorylates Ser9, which folds the N-terminal tail back into the
substrate groove as a pseudosubstrate. It prefers substrates that have
already been "primed" by another kinase at the S/T-X-X-X-S(p) motif, and the
primed phosphate is read by a positively charged pocket (Arg96, Arg180,
Lys205) next to the catalytic site. Its substrates include tau, glycogen
synthase, β-catenin, CRMP2 and the translation factor eIF2B, so it is a hub
rather than a private tau kinase - which is the core of the selectivity
problem discussed in the limitations.

## 1.3 Evidence linking GSK-3β to Alzheimer's disease

1. **Hooper C., Killick R., Lovestone S. (2008)** *The GSK3 hypothesis of
   Alzheimer's disease.* J. Neurochem. **104**, 1433-1439.
   Collects the evidence that GSK-3 activity is increased in AD brain, that
   it phosphorylates tau at PHF epitopes, that it is activated by Aβ, and
   that it mediates Aβ-induced neurotoxicity and LTP impairment - i.e. it
   places GSK-3 exactly at the amyloid→tau junction of our cascade.

2. **Lucas J.J. et al. (2001)** *Decreased nuclear β-catenin, tau
   hyperphosphorylation and neurodegeneration in GSK-3β conditional
   transgenic mice.* EMBO J. **20**, 27-39.
   Conditional forebrain overexpression of GSK-3β alone is sufficient to
   produce hyperphosphorylated, somatodendritically mislocalised tau,
   reactive astrocytosis and neuronal death. This is the causal arm of the
   argument: raise the kinase, get the node-4 phenotype.

3. **Phiel C.J., Wilson C.A., Lee V.M.-Y., Klein P.S. (2003)** *GSK-3α
   regulates production of Alzheimer's disease amyloid-β peptides.* Nature
   **423**, 435-439.
   Lithium and GSK-3 knockdown lower Aβ production via γ-secretase-mediated
   APP cleavage. The paralogue specificity (α over β for Aβ production) is
   itself informative: it means a GSK-3β-directed agent is aimed at the tau
   arm, and tells us that paralogue selectivity has to be measured, not
   assumed.

Supporting: **Bhat R. et al. (2003)** J. Biol. Chem. **278**, 45937-45945,
which reports the ATP-competitive inhibitor AR-A014418 (IC50 104 ± 27 nM,
Ki 38 nM, no inhibition of CDK2/CDK5 at 100 µM) and shows it blocks tau
Ser396 phosphorylation in cells and protects neurons from Aβ toxicity. The
same compound is the ligand in PDB 1Q5K, which we use below as the primary
template, as the source of the binding-site definition, and as the redocking
control.

## 1.4 Intended action

**Inhibit — ATP-competitive, at the kinase active site.**

Rationale within the cascade: node 4 is a kinase/phosphatase balance
(GSK-3β and CDK5 up, PP2A down). Lowering GSK-3β catalytic output shifts the
balance back, so tau stays bound to microtubules, tangle formation and
trans-synaptic spread are slowed, and axonal transport is preserved. The
expected downstream effect is on node 5 (synapse and neuron loss), not on
nodes 1-3: this intervention does not claim to lower Aβ. It claims to cut
the wire between the amyloid signal and the tau lesion.

Why the ATP site rather than an allosteric or substrate-competitive site:
(i) the ATP cleft is the only pocket on GSK-3β with a validated, druggable,
well-characterised geometry (123 PDB entries, almost all ATP-site complexes);
(ii) the substrate-priming pocket is a shallow, highly charged groove that
is poorly suited to a small neutral molecule that must also cross the
blood-brain barrier; (iii) the clinically tested non-ATP-competitive agent
tideglusib failed to separate from placebo. The cost of this choice is
selectivity against the ~500-member kinome, which is treated explicitly in
the limitations.

A partial, not complete, inhibition is the aim. Complete GSK-3 blockade
disinhibits Wnt/β-catenin signalling and carries an oncogenic and metabolic
liability; the therapeutic hypothesis is a modest, chronic reduction in
kinase output, exactly as the field concluded for BACE1 after verubecestat.

## 1.5 Existing drugs and tool compounds (positive controls)

There is **no approved selective GSK-3 inhibitor**. The chemical matter that
exists is tool compounds and one failed clinical candidate, and these become
our positive controls (all potencies are ChEMBL medians over curated
binding assays against CHEMBL262 = human GSK-3β; see
`04_ligands/actives/known_actives.csv`):

| compound | ChEMBL | median potency | status / mechanism |
|---|---|---|---|
| LY-2090314 | CHEMBL362558 | 1.1 nM | clinical (oncology), ATP-competitive maleimide |
| Alsterpaullone | CHEMBL50894 | 4.0 nM | tool, ATP-competitive paullone |
| CHIR-99021 (laduviglusib) | CHEMBL412142 | 17.4 nM | tool, highly selective aminopyrimidine |
| Kenpaullone | CHEMBL296586 | 22.9 nM | tool, paullone |
| Staurosporine | CHEMBL388978 | 24.0 nM | pan-kinase reference |
| SB-216763 | CHEMBL102714 | 33.9 nM | tool, maleimide |
| Indirubin-3'-monoxime | CHEMBL1496098 | 55.0 nM | natural-product-derived bisindole |
| SB-415286 | CHEMBL322970 | 77.6 nM | tool, maleimide |
| Tideglusib | CHEMBL3545157 | 102.3 nM | phase II in AD (ARGO), **non**-ATP-competitive |
| AR-A014418 | CHEMBL259850 | 104.7 nM | tool, ATP-competitive, co-crystal 1Q5K |
| TDZD-8 | CHEMBL284861 | 1.74 µM | tool, non-ATP-competitive thiadiazolidinedione |

Clinical precedent and its lesson: **Lovestone S. et al. (2015)** *A phase II
trial of tideglusib in Alzheimer's disease.* J. Alzheimers Dis. **45**,
75-88. 306 patients, 26 weeks, no significant difference from placebo on
ADAS-cog15; dose-dependent reversible transaminase rises. Lithium, which is
approved but inhibits GSK-3 only at millimolar, non-selective concentrations,
is the other precedent. Two of our controls (tideglusib, TDZD-8) are
non-ATP-competitive and are therefore expected to behave *badly* in an
ATP-site docking run - they are kept in deliberately, as a negative control
on the protocol itself.
