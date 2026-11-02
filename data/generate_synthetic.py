"""
Synthetic RAW molecular dataset generator for OncoKG.

Simulates 600 patients profiled at three laboratories (SITE_A/B/C) that each export
data in their own conventions -- different marker names, units, scoring systems,
HGVS styles, VAF scales, copy-number metrics and RNA gene identifiers -- plus
realistic noise: missing tests, low tumour purity, CHIP-like low-VAF variants,
VUS, low-read fusion artefacts and IHC/RNA discordance.

Outputs (data/raw/):
  clinical.csv, ihc.csv, ish.csv, ngs_variants.csv, ngs_fusions.csv, ngs_cnv.csv,
  ngs_summary.csv, serum.csv, rnaseq_counts_SITE_A.csv, _SITE_B.csv, _SITE_C.csv
and data/ground_truth.csv (hidden labels, used only for evaluation).
"""
import re
import json
from pathlib import Path
import numpy as np
import pandas as pd
from rdflib import Graph, Namespace, RDF
from rdflib.namespace import SKOS

ROOT = Path(__file__).resolve().parents[1]
RAW = ROOT / "data" / "raw"
RAW.mkdir(parents=True, exist_ok=True)
rng = np.random.default_rng(20261102)
SITES = ["SITE_A", "SITE_B", "SITE_C"]

# ---------------------------------------------------------------- RNA panel from ontology
OKG = Namespace("http://example.org/oncokg#")
G = Graph().parse(ROOT / "ontology" / "oncokg.ttl")
GENES = {}
for gnode in G.subjects(RDF.type, OKG.Gene):
    sym = str(G.value(gnode, OKG.hgncSymbol))
    ln = G.value(gnode, OKG.effectiveLengthBp)
    if ln is None:
        continue
    aliases = [str(a) for a in G.objects(gnode, SKOS.altLabel)]
    hk = (gnode, RDF.type, OKG.HousekeepingGene) in G
    GENES[sym] = dict(length=int(ln), alias=sorted(aliases)[0] if aliases else sym, hk=hk)
HK_LEVEL = {"ACTB": 13, "GAPDH": 13, "RPLP0": 12, "PGK1": 9, "HPRT1": 7, "TBP": 6}  # mean = 10
SITE_GENE_BIAS = {g: rng.normal(0, 0.25) for g in GENES}  # SITE_C probe bias

AA3 = dict(A="Ala", R="Arg", N="Asn", D="Asp", C="Cys", E="Glu", Q="Gln", G="Gly", H="His", I="Ile",
           L="Leu", K="Lys", M="Met", F="Phe", P="Pro", S="Ser", T="Thr", W="Trp", Y="Tyr", V="Val")
AA3["*"] = "Ter"


def U(a, b):
    return float(rng.uniform(a, b))


def p(x):
    return rng.random() < x


def pick(options, weights=None):
    i = rng.choice(len(options), p=None if weights is None else np.array(weights) / sum(weights))
    return options[i]


# ---------------------------------------------------------------- variant catalog
TP53_MIS = ["p.R175H", "p.R248Q", "p.R273H", "p.R282W", "p.Y220C", "p.G245S"]
TP53_TRUNC = ["p.R213*", "p.R342*", "p.E298*"]
TRUNC = {
    "BRCA1": ["p.E23Vfs*17", "p.Q1756Pfs*74", "p.R1443*"],
    "BRCA2": ["p.S1982Rfs*22", "p.N1784Kfs*3", "p.R3128*"],
    "ATM": ["p.E1978*", "p.R2443*"], "CDK12": ["p.Q1032*", "p.K765Rfs*12"], "PALB2": ["p.Q775*"],
    "RB1": ["p.R320*", "p.R579*", "p.Q702*"], "PTEN": ["p.R130*", "p.R233*"],
    "STK11": ["p.Q37*", "p.E57Gfs*7"], "KEAP1": ["p.G333C", "p.R470C"],
}
VUS = [("BRCA2", "p.T1915M"), ("ATM", "p.S707P"), ("BRCA1", "p.R1347G"), ("PALB2", "p.L939W"), ("CDK12", "p.A1015T")]
CHIP = [("DNMT3A", "p.R882H"), ("TET2", "p.Q1034*"), ("TP53", "p.R248W")]
EX19 = ["p.E746_A750del", "p.L747_P753delinsS", "p.E746_S752delinsV"]


def var(gene, prot, cls="Pathogenic", vtype=None, cdna="", vaf=None):
    if vtype is None:
        vtype = "frameshift" if "fs" in prot else "nonsense" if prot.endswith("*") else \
            "inframe_deletion" if "del" in prot else "missense"
    return dict(gene=gene, prot=prot, cls=cls, vtype=vtype, cdna=cdna, vaf=vaf)


def tp53():
    return var("TP53", pick(TP53_MIS + TP53_TRUNC))


# ---------------------------------------------------------------- patient profile
def base(organ, sex=None):
    P = dict(organ=organ, sex=sex or ("F" if organ == "breast" else "M" if organ == "prostate" else pick(["M", "F"], [55, 45])),
             age=int(np.clip(rng.normal(64, 10), 30, 90)), site=pick(SITES), purity=U(0.3, 0.85),
             specimen=None, history="", pack_years=None, ihc={}, her2=None, her2_ish=None, alk_fish=None,
             variants=[], fusions=[], cn={}, fga=U(0.12, 0.55), tmb=U(1, 4), rna={}, serum={},
             truth=dict(lineage=organ, histology=None, molecular="", modifiers=[], edge="", note=""))
    if p(0.05):
        P["purity"] = U(0.1, 0.19)
        P["truth"]["note"] += "low purity; "
    return P


def status_markers(P, spec):
    """spec: {marker: (p_pos, p_weak)}; others negative."""
    for m, (pp, pw) in spec.items():
        r = rng.random()
        P["ihc"][m] = "pos" if r < pp else "weak" if r < pp + pw else "neg"


ALL_STATUS = ["TTF1", "NAPSA", "P40", "CK56", "GATA3", "TRPS1", "MGB", "PSA", "NKX31", "PSMA", "AMACR", "HOXB13",
              "SYP", "CHGA", "INSM1"]


def default_neg(P):
    for m in ALL_STATUS:
        P["ihc"].setdefault(m, "neg")


def specimen(P, primary_names, met_names, p_met=0.15):
    if p(p_met):
        P["specimen"] = pick(met_names)
        P["truth"]["note"] += "metastatic specimen; "
    else:
        P["specimen"] = pick(primary_names)


LUNG_SITES = ["Lung, RUL core", "Lung LLL biopsy", "Bronchial biopsy", "lung wedge resection"]
LUNG_METS = ["Pleural fluid cell block", "LN, mediastinal (EBUS)", "Liver", "Brain"]
BREAST_SITES = ["Breast, left, core", "R breast core bx", "Breast lumpectomy"]
BREAST_METS = ["LN axillary", "Liver", "Bone (iliac)", "Brain"]
PROSTATE_SITES = ["Prostate core bx", "prostate needle biopsy", "Prostatectomy"]
PROSTATE_METS = ["Bone (iliac)", "LN, pelvic", "Liver"]


# ---- LUNG -------------------------------------------------------------------
def lung(P, kind):
    T = P["truth"]
    smoker = kind in ("LUSC", "SCLC", "KRASG12C", "MET", "NOS") or p(0.35)
    if kind in ("EGFR", "ALK", "ROS1", "RET", "NTRK"):
        smoker = p(0.25)
        P["age"] = int(np.clip(rng.normal(58, 11), 28, 88))
    P["pack_years"] = int(U(20, 70)) if smoker else (0 if p(0.8) else int(U(1, 10)))
    if smoker:
        P["tmb"] = U(6, 22)
    specimen(P, LUNG_SITES, LUNG_METS)
    P["ihc"]["Ki67"] = U(10, 60)
    r = rng.random()
    P["ihc"]["PDL1"] = 0 if r < 0.3 else int(U(1, 49)) if r < 0.7 else int(U(50, 95))
    T["modifiers"].append("MOD_PDL1_NEG" if P["ihc"]["PDL1"] < 1 else "MOD_PDL1_LOW" if P["ihc"]["PDL1"] < 50 else "MOD_PDL1_HIGH")
    P["serum"]["CEA"] = U(6, 60) if p(0.5) else U(1, 4.5)
    P["serum"]["CYFRA21_1"] = U(1, 3)
    P["serum"]["ProGRP"] = U(20, 60)
    P["serum"]["NSE"] = U(8, 15)

    if kind == "SCLC":
        T["lineage"], T["histology"] = "lung", "SmallCellLungCancer"
        status_markers(P, {"TTF1": (0.85, 0.05), "SYP": (0.95, 0.03), "CHGA": (0.8, 0.1), "INSM1": (0.9, 0.05)})
        P["rna"]["ASCL1"] = "on" if p(0.8) else "off"
        P["ihc"]["Ki67"] = U(60, 95)
        P["ihc"].pop("PDL1")
        T["modifiers"] = [m for m in T["modifiers"] if not m.startswith("MOD_PDL1")]
        if p(0.9):
            T["modifiers"].append("MOD_RB1")
            if p(0.4):
                P["cn"]["RB1"] = 0
            else:
                P["variants"].append(var("RB1", pick(TRUNC["RB1"])))
        if p(0.9):
            P["variants"].append(tp53())
        P["serum"]["ProGRP"] = U(100, 3000) if p(0.7) else U(30, 60)
        P["serum"]["NSE"] = U(20, 150) if p(0.6) else U(9, 15)
        P["pack_years"] = int(U(30, 80))
    elif kind == "LUSC":
        T["histology"] = "LungSquamousCellCarcinoma"
        status_markers(P, {"P40": (0.98, 0.0), "CK56": (0.9, 0.05), "TTF1": (0.0, 0.03)})
        P["rna"]["SOX2"] = "on"
        P["serum"]["CYFRA21_1"] = U(4, 30) if p(0.6) else U(1, 3)
        if p(0.8): P["variants"].append(tp53())
        if p(0.15): P["variants"].append(var("PIK3CA", "p.E545K"))
        if p(0.15): P["variants"].append(var("KEAP1", pick(TRUNC["KEAP1"]))); T["modifiers"].append("MOD_KEAP1")
    else:  # adenocarcinoma family
        T["histology"] = "LungAdenocarcinoma"
        if p(0.04):
            status_markers(P, {"TTF1": (0, 0), "NAPSA": (0, 0)})
            T["edge"] = "TTF-1/Napsin A double-negative adenocarcinoma (expect NSCLC-NOS)"
        else:
            status_markers(P, {"TTF1": (0.85, 0.1), "NAPSA": (0.75, 0.1), "P40": (0.0, 0.05)})
        if P["ihc"].get("TTF1") == "neg" and P["ihc"].get("NAPSA") == "neg" and not T["edge"]:
            P["ihc"]["TTF1"] = "pos"
        if p(0.05): P["ihc"]["SYP"] = "weak"
        if kind != "ALK" and p(0.5): P["ihc"]["ALK"] = "neg"
        P["rna"]["SFTPC"] = "on" if p(0.6) else "off"
        mol = {"EGFR": "EGFRMutantNSCLC", "ALK": "ALKRearrangedNSCLC", "ROS1": "ROS1RearrangedNSCLC",
               "RET": "RETRearrangedNSCLC", "NTRK": "NTRKFusionNSCLC", "MET": "METEx14SkippingNSCLC",
               "BRAF": "BRAFV600EMutantNSCLC", "KRASG12C": "KRASG12CMutantNSCLC"}.get(kind, "")
        T["molecular"] = mol
        if kind == "EGFR":
            prot = pick(["p.L858R", pick(EX19), pick(["p.G719S", "p.L861Q"])], [45, 45, 10])
            P["variants"].append(var("EGFR", prot))
            if p(0.08):
                P["variants"].append(var("EGFR", "p.T790M")); T["modifiers"].append("MOD_EGFR_T790M")
            if p(0.5): P["variants"].append(tp53())
        elif kind == "ALK":
            P["ihc"]["ALK"] = "pos" if p(0.9) else None
            P["fusions"].append(("EML4", "ALK", int(U(15, 300))))
            if p(0.6): P["alk_fish"] = U(18, 70)
            P["rna"]["ALK"] = "on"
            P["ihc"]["PDL1"] = int(U(0, 30)); T["modifiers"] = [m for m in T["modifiers"] if not m.startswith("MOD_PDL1")]
            T["modifiers"].append("MOD_PDL1_NEG" if P["ihc"]["PDL1"] < 1 else "MOD_PDL1_LOW")
        elif kind == "ROS1":
            P["fusions"].append((pick(["CD74", "SLC34A2", "EZR", "SDC4"]), "ROS1", int(U(10, 200))))
        elif kind == "RET":
            P["fusions"].append((pick(["KIF5B", "CCDC6"]), "RET", int(U(10, 150))))
        elif kind == "NTRK":
            P["fusions"].append(pick([("ETV6", "NTRK3"), ("TPM3", "NTRK1"), ("SQSTM1", "NTRK1")]) + (int(U(10, 120)),))
        elif kind == "MET":
            P["age"] = int(np.clip(rng.normal(74, 6), 55, 92))
            P["variants"].append(var("MET", "exon14skip", vtype="splice", cdna=pick(["c.3082+1G>T", "c.3028+2T>C", "c.3082+2T>C"])))
        elif kind == "BRAF":
            P["variants"].append(var("BRAF", "p.V600E"))
        elif kind == "KRASG12C":
            P["variants"].append(var("KRAS", "p.G12C"))
            if p(0.3): P["variants"].append(var("STK11", pick(TRUNC["STK11"]))); T["modifiers"].append("MOD_STK11")
            if p(0.2): P["variants"].append(var("KEAP1", pick(TRUNC["KEAP1"]))); T["modifiers"].append("MOD_KEAP1")
            if p(0.4): P["variants"].append(tp53())
        elif kind == "NOS":
            if p(0.4):
                P["variants"].append(var("KRAS", pick(["p.G12D", "p.G12V", "p.G13C", "p.Q61H"]))); T["modifiers"].append("MOD_KRAS_OTHER")
            if p(0.5): P["variants"].append(tp53())
        if p(0.06) and kind != "ALK":  # sub-threshold artefact read
            P["fusions"].append(("EML4", "ALK", int(U(1, 4))))
    default_neg(P)


# ---- BREAST -----------------------------------------------------------------
def her2_negative_profile(P, weights=(20, 20, 40, 20)):
    k = pick(["Null", "Ultralow", "Low1", "Low2"], list(weights))
    if k == "Null":
        P["her2"] = dict(ihc=0, faint=None)
        P["truth"]["modifiers"].append("MOD_HER2_NULL")
    elif k == "Ultralow":
        P["her2"] = dict(ihc=0, faint=int(U(2, 10)))
        P["truth"]["modifiers"].append("MOD_HER2_ULTRALOW")
    elif k == "Low1":
        P["her2"] = dict(ihc=1, faint=None)
        P["truth"]["modifiers"].append("MOD_HER2_LOW")
    else:
        P["her2"] = dict(ihc=2, faint=None)
        grp = pick(["G5", "G2", "G4"], [80, 10, 10])
        P["her2_ish"] = {"G5": (U(1.0, 1.7), U(1.8, 3.6)), "G2": (U(2.05, 2.6), U(2.6, 3.8)),
                         "G4": (U(1.2, 1.8), U(4.2, 5.8))}[grp]
        P["truth"]["modifiers"].append("MOD_HER2_LOW")
        if grp != "G5":
            P["truth"]["note"] += f"HER2 ISH group {grp[1]}; "
    P["cn"]["ERBB2"] = 2 if P["her2"]["ihc"] < 2 else int(pick([2, 3]))


def her2_positive_profile(P):
    if p(0.8):
        P["her2"] = dict(ihc=3, faint=None)
        if p(0.2): P["her2_ish"] = (U(3, 8), U(8, 20))
    else:
        P["her2"] = dict(ihc=2, faint=None)
        if p(0.85):
            P["her2_ish"] = (U(2.2, 4.5), U(4.5, 10))
        else:
            P["her2_ish"] = (U(1.6, 1.95), U(6.2, 8.0))
            P["truth"]["note"] += "HER2 ISH group 3; "
    P["cn"]["ERBB2"] = int(U(8, 25))
    P["rna"]["ERBB2"] = "amp"


def breast(P, kind):
    T = P["truth"]
    T["lineage"] = "breast"
    specimen(P, BREAST_SITES, BREAST_METS)
    if P["specimen"] in BREAST_METS and p(0.6):
        P["history"] = pick(["Hx breast ca 2019, s/p lumpectomy", "prior breast carcinoma (2018)", "history of breast cancer"])
    P["serum"]["CEA"] = U(1, 4.5) if p(0.7) else U(6, 30)
    P["serum"]["CA15_3"] = U(40, 200) if P["specimen"] in BREAST_METS else U(10, 28)
    lum = {"GATA3": (0.95, 0.03), "TRPS1": (0.95, 0.03), "MGB": (0.7, 0.1)}
    tn = {"GATA3": (0.55, 0.15), "TRPS1": (0.85, 0.05), "MGB": (0.25, 0.1)}
    if p(0.12): P["variants"].append(var("GATA3", "p.D336Gfs*17", cls="Likely pathogenic"))

    def pik(prob):
        if p(prob):
            P["variants"].append(var("PIK3CA", pick(["p.H1047R", "p.E545K", "p.E542K", "p.H1047L"]))); T["modifiers"].append("MOD_PIK3CA")

    def tp(prob):
        if p(prob): P["variants"].append(tp53()); T["modifiers"].append("MOD_TP53")

    def brca(prob, gene="BRCA1"):
        if p(prob):
            P["variants"].append(var(gene, pick(TRUNC[gene]))); T["modifiers"].append(f"MOD_{gene}")

    if kind == "LUMA":
        T["histology"] = "LuminalALike"
        status_markers(P, lum)
        P["ihc"].update(ER=U(50, 100), PR=U(25, 100), Ki67=U(3, 18), AR=U(30, 90))
        her2_negative_profile(P); pik(0.45); tp(0.1)
        if p(0.12):
            P["fga"] = U(0.02, 0.08); T["note"] += "near-diploid; "
    elif kind == "LUMB":
        T["histology"] = "LuminalBLikeHER2Negative"
        status_markers(P, lum)
        if p(0.6):
            P["ihc"].update(ER=U(30, 100), PR=U(0, 100), Ki67=U(22, 50))
        else:
            P["ihc"].update(ER=U(30, 100), PR=U(1, 15), Ki67=U(8, 18))
        P["ihc"]["AR"] = U(20, 80)
        her2_negative_profile(P); pik(0.35); tp(0.3)
        if p(0.06): P["variants"].append(var("ESR1", pick(["p.Y537S", "p.D538G"]))); T["modifiers"].append("MOD_ESR1")
    elif kind == "LUMB_HER2P":
        T["histology"] = "LuminalBLikeHER2Positive"
        status_markers(P, lum)
        P["ihc"].update(ER=U(20, 100), PR=U(0, 60), Ki67=U(20, 60), AR=U(10, 60))
        her2_positive_profile(P); pik(0.25); tp(0.5)
    elif kind == "HER2E":
        T["histology"] = "HER2Enriched"
        status_markers(P, {"GATA3": (0.85, 0.1), "TRPS1": (0.9, 0.05), "MGB": (0.4, 0.1)})
        P["ihc"].update(ER=0.0, PR=0.0, Ki67=U(30, 70), AR=U(10, 70))
        her2_positive_profile(P); tp(0.7)
    elif kind in ("TNBC", "LAR"):
        T["histology"] = "TripleNegativeBreastCancer"
        status_markers(P, tn)
        P["ihc"].update(ER=0.0 if p(0.8) else 0.5, PR=0.0, Ki67=U(30, 90))
        her2_negative_profile(P, (40, 25, 25, 10))
        P["fga"] = U(0.25, 0.65)
        P["tmb"] = U(2, 6)
        if kind == "LAR":
            T["molecular"] = "TNBCLuminalAndrogenReceptor"
            P["ihc"]["AR"] = U(30, 90); P["ihc"]["Ki67"] = U(15, 40)
            status_markers(P, {"GATA3": (0.9, 0.05), "TRPS1": (0.9, 0.05), "MGB": (0.5, 0.1)})
            pik(0.5); tp(0.4)
        else:
            P["ihc"]["AR"] = U(0, 5)
            tp(0.85)
            if p(0.25):
                brca(1.0, "BRCA1"); T["note"] += "BRCA1-associated TNBC; "
            elif p(0.06):
                brca(1.0, "BRCA2")
            if p(0.1):
                P["cn"]["RB1"] = 0; T["modifiers"].append("MOD_RB1")
    elif kind == "ERLOW":
        T["histology"] = "ERLowHER2NegativeBreastCancer"
        status_markers(P, tn)
        P["ihc"].update(ER=U(2, 9), PR=U(0, 5), Ki67=U(25, 70), AR=U(0, 30))
        T["modifiers"].append("MOD_ER_LOW")
        her2_negative_profile(P, (35, 25, 30, 10)); tp(0.6); brca(0.1)
    elif kind == "ERNEG_PRPOS":
        T["histology"] = "BreastCancer"
        T["edge"] = "ER-/PR+ discordant phenotype (expect flag + retest)"
        status_markers(P, lum)
        P["ihc"].update(ER=0.0, PR=U(15, 60), Ki67=U(15, 40), AR=U(5, 40))
        her2_negative_profile(P); tp(0.3)
    elif kind == "HER2_EQUIVOCAL":
        T["histology"] = "BreastCancer"
        T["edge"] = "HER2 IHC 2+ with ISH pending (expect HER2-equivocal flag)"
        status_markers(P, lum)
        P["ihc"].update(ER=U(40, 100), PR=U(10, 80), Ki67=U(15, 40), AR=U(20, 70))
        P["her2"] = dict(ihc=2, faint=None)
        P["cn"]["ERBB2"] = int(pick([2, 3, 7]))
    if T["histology"] in ("LuminalALike", "LuminalBLikeHER2Negative", "LuminalBLikeHER2Positive"):
        pass
    default_neg(P)


# ---- PROSTATE ---------------------------------------------------------------
def prostate(P, kind):
    T = P["truth"]
    T["lineage"] = "prostate"
    P["age"] = int(np.clip(rng.normal(68, 8), 45, 92))
    P["tmb"] = U(0.5, 3)
    P["fga"] = U(0.1, 0.4)
    specimen(P, PROSTATE_SITES, PROSTATE_METS, p_met=0.25 if kind != "NEPC" else 0.0)
    P["serum"]["PSA"] = float(np.exp(rng.normal(np.log(14), 1.0))) + 4.2
    if kind in ("ADENO", "HRR"):
        T["histology"] = "ProstateAdenocarcinoma"
        status_markers(P, {"PSA": (0.95, 0.03), "NKX31": (0.95, 0.03), "PSMA": (0.85, 0.05), "AMACR": (0.9, 0.05),
                           "HOXB13": (0.9, 0.05)})
        P["ihc"].update(AR=U(60, 100), Ki67=U(3, 25))
        if p(0.45):
            P["fusions"].append(("TMPRSS2", "ERG", int(U(10, 300)))); P["rna"]["ERG"] = "on"; T["modifiers"].append("MOD_TMPRSS2_ERG")
        elif p(0.18):
            P["variants"].append(var("SPOP", pick(["p.F133V", "p.W131G", "p.Y87C"]))); T["modifiers"].append("MOD_SPOP")
        if p(0.2):
            P["cn"]["PTEN"] = 0; T["modifiers"].append("MOD_PTEN")
        if p(0.15):
            P["variants"].append(tp53()); T["modifiers"].append("MOD_TP53")
        if kind == "HRR":
            T["molecular"] = "HRRDeficientProstateCancer"
            gene = pick(["BRCA2", "ATM", "CDK12", "PALB2", "BRCA1"], [55, 20, 15, 5, 5])
            P["variants"].append(var(gene, pick(TRUNC[gene])))
            T["modifiers"].append(f"MOD_{gene}" if gene in ("BRCA1", "BRCA2", "ATM", "CDK12") else "")
            T["modifiers"] = [m for m in T["modifiers"] if m]
    elif kind == "NEPC":
        T["histology"] = "NeuroendocrineProstateCancer"
        P["history"] = pick(["prior prostate adenocarcinoma (2021), ADT + enzalutamide",
                             "Hx prostate ca 2019 Gleason 4+4, on ADT", "mCRPC on abiraterone, prostate adenocarcinoma 2020"])
        P["specimen"] = pick(["Liver", "Prostate core bx", "LN, pelvic", "Bone (iliac)"], [40, 30, 20, 10])
        status_markers(P, {"SYP": (0.95, 0.03), "CHGA": (0.8, 0.1), "INSM1": (0.85, 0.05), "PSA": (0.0, 0.2),
                           "NKX31": (0.0, 0.3), "TTF1": (0.3, 0.1), "HOXB13": (0.5, 0.15), "PSMA": (0.1, 0.1)})
        P["rna"]["ASCL1"] = "on" if p(0.8) else "off"
        P["ihc"].update(AR=U(0, 6), Ki67=U(50, 90))
        P["serum"]["PSA"] = U(0.05, 3.5)
        P["serum"]["NSE"] = U(20, 120) if p(0.6) else U(8, 15)
        P["tmb"] = U(2, 5)
        if p(0.8):
            T["modifiers"].append("MOD_RB1")
            if p(0.5): P["cn"]["RB1"] = 0
            else: P["variants"].append(var("RB1", pick(TRUNC["RB1"])))
        if p(0.7): P["variants"].append(tp53()); T["modifiers"].append("MOD_TP53")
        if p(0.45):
            P["fusions"].append(("TMPRSS2", "ERG", int(U(10, 200)))); P["rna"]["ERG"] = "on"; T["modifiers"].append("MOD_TMPRSS2_ERG")
        if p(0.4): P["cn"]["MYCN"] = int(U(8, 30)); T["modifiers"].append("MOD_MYCN")
        if p(0.4): P["cn"]["AURKA"] = int(U(7, 20)); T["modifiers"].append("MOD_AURKA")
    default_neg(P)


# ---- BENIGN CONTROLS --------------------------------------------------------
def benign(P, organ):
    T = P["truth"]
    T["lineage"], T["histology"] = organ, "NoMalignancyDetected"
    P["purity"] = None
    P["fga"], P["tmb"] = U(0.0, 0.03), U(0.0, 0.8)
    P["ihc"]["Ki67"] = U(1, 5)
    if organ == "lung":
        P["specimen"] = pick(LUNG_SITES)
        status_markers(P, {"TTF1": (0.95, 0.05), "NAPSA": (0.9, 0.05)})
        P["rna"]["SFTPC"] = "on"
        P["pack_years"] = int(U(0, 40))
        P["serum"].update(CEA=U(1, 4.5), CYFRA21_1=U(1, 3), ProGRP=U(20, 60), NSE=U(8, 15))
        T["note"] = "benign lung nodule (granuloma/hamartoma)"
    elif organ == "breast":
        P["specimen"] = pick(BREAST_SITES)
        status_markers(P, {"GATA3": (0.95, 0.05), "TRPS1": (0.95, 0.05), "MGB": (0.4, 0.2), "CK56": (0.6, 0.3)})
        P["ihc"].update(ER=U(12, 40), PR=U(10, 40), AR=U(5, 30))
        P["her2"] = dict(ihc=0, faint=None)
        P["serum"].update(CEA=U(1, 4.5), CA15_3=U(10, 28))
        T["note"] = "benign breast (fibroadenoma / fibrocystic)"
    else:
        P["specimen"] = pick(PROSTATE_SITES[:2])
        status_markers(P, {"PSA": (0.98, 0.02), "NKX31": (0.98, 0.02), "PSMA": (0.3, 0.3), "HOXB13": (0.9, 0.05)})
        P["ihc"]["AR"] = U(50, 95)
        P["serum"]["PSA"] = U(4.2, 12)
        T["note"] = "benign prostatic hyperplasia with elevated serum PSA"
    default_neg(P)


COHORT = [
    ("lung", "EGFR", 40), ("lung", "KRASG12C", 30), ("lung", "ALK", 12), ("lung", "ROS1", 6), ("lung", "BRAF", 6),
    ("lung", "MET", 8), ("lung", "RET", 5), ("lung", "NTRK", 3), ("lung", "NOS", 30), ("lung", "LUSC", 40), ("lung", "SCLC", 30),
    ("breast", "LUMA", 75), ("breast", "LUMB", 45), ("breast", "LUMB_HER2P", 25), ("breast", "HER2E", 20),
    ("breast", "TNBC", 35), ("breast", "LAR", 10), ("breast", "ERLOW", 10), ("breast", "ERNEG_PRPOS", 3),
    ("breast", "HER2_EQUIVOCAL", 5),
    ("prostate", "ADENO", 80), ("prostate", "HRR", 25), ("prostate", "NEPC", 15),
    ("benign", "lung", 12), ("benign", "breast", 14), ("benign", "prostate", 14),
]

patients = []
for organ, kind, n in COHORT:
    for _ in range(n):
        if organ == "benign":
            P = base(kind); benign(P, kind)
        else:
            P = base(organ)
            if organ == "breast" and kind == "LUMA" and sum(1 for q in patients if q["sex"] == "M" and q["organ"] == "breast") < 2 and p(0.05):
                P["sex"] = "M"; P["truth"]["note"] += "male breast cancer; "
            {"lung": lung, "breast": breast, "prostate": prostate}[organ](P, kind)
        P["scenario"] = f"{organ}:{kind}"
        if any(v["gene"] == "TP53" and v["cls"] != "VUS" for v in P["variants"]):
            P["truth"]["modifiers"].append("MOD_TP53")
        patients.append(P)

order = rng.permutation(len(patients))
patients = [patients[i] for i in order]
for i, P in enumerate(patients):
    P["pid"] = f"PT{i + 1:04d}"

# ensure some lung cases are labelled with their histology for "NOS" kind
for P in patients:
    T = P["truth"]
    if T["histology"] is None:
        T["histology"] = "LungAdenocarcinoma"


# ---------------------------------------------------------------- renderers (site-specific raw formats)
STAIN_NAMES = {
    "TTF1": ["TTF-1", "TTF1", "NKX2-1"], "NAPSA": ["Napsin A", "NAPSA", "napsin-a"], "P40": ["p40", "TP63", "deltaNp63"],
    "CK56": ["CK5/6", "CK5/6", "KRT5"], "GATA3": ["GATA3", "GATA-3", "GATA3"], "TRPS1": ["TRPS1"] * 3,
    "MGB": ["Mammaglobin", "SCGB2A2", "MGB1"], "ER": ["ER", "Estrogen receptor", "ESR1"], "PR": ["PR", "PgR", "Progesterone receptor"],
    "HER2": ["HER2", "HER-2/neu", "c-erbB-2"], "Ki67": ["Ki-67", "MIB-1", "Ki67"], "AR": ["AR", "Androgen receptor", "AR"],
    "PSA": ["PSA", "Prostate specific antigen", "KLK3"], "NKX31": ["NKX3.1", "NKX3-1", "NKX3.1"], "PSMA": ["PSMA", "FOLH1", "PSMA"],
    "AMACR": ["AMACR", "P504S", "AMACR"], "HOXB13": ["HOXB13"] * 3, "SYP": ["Synaptophysin", "SYP", "Synaptophysin"],
    "CHGA": ["Chromogranin A", "CgA", "CHGA"], "INSM1": ["INSM1"] * 3, "PDL1": ["PD-L1 (22C3)", "PD-L1", "CD274"],
    "ALK": ["ALK (D5F3)", "ALK D5F3", "ALK"],
}
STATUS_TXT = {"pos": ["Positive", "POS (diffuse)", "positive, strong diffuse"],
              "weak": ["Weak positive", "focal +", "focally positive, weak"],
              "neg": ["Negative", "NEG", "no staining"]}


def si(site):
    return SITES.index(site)


def allred(pct):
    ps = 0 if pct == 0 else 1 if pct < 1 else 2 if pct <= 10 else 3 if pct <= 33 else 4 if pct <= 66 else 5
    inten = 0 if ps == 0 else int(pick([1, 2, 3], [20, 40, 40]))
    return ps, inten


def fmt_receptor(pct, site, name):
    s = si(site)
    if s == 0:
        return f"{pct:.0f}%" if pct >= 1 or pct == 0 else "<1%"
    if s == 1:
        ps, it = allred(pct)
        return f"Allred PS{ps} IS{it} (TS {ps + it})"
    if pct < 1:
        return "Negative (<1%)" if pct > 0 else "Negative, 0%"
    if pct <= 10:
        return f"Low positive, {pct:.0f}%"
    return f"Positive, {pct:.0f}%, {pick(['strong', 'moderate'])}"


def fmt_ihc_value(m, v, site, P):
    s = si(site)
    if m in ("ER", "PR"):
        return fmt_receptor(v, site, m)
    if m == "Ki67":
        return [f"{v:.0f}%", f"MIB-1 index {v / 100:.2f}", f"Ki-67: {v:.0f}% (hot spot)"][s]
    if m == "AR":
        if s == 2:
            return f"AR negative ({v:.0f}%)" if v < 1 else f"AR positive, {v:.0f}%"
        return [f"{v:.0f}%", f"{v:.0f}% nuclear"][s]
    if m == "PDL1":
        if v < 1:
            return ["TPS <1%", "TPS <1%", "22C3 TPS: <1%"][s]
        return [f"TPS {v:.0f}%", f"{v:.0f}", f"22C3 TPS: {v:.0f}%"][s]
    return STATUS_TXT[v][s]


def fmt_her2(h, site):
    s = si(site)
    sc, faint = h["ihc"], h["faint"]
    if sc == 0:
        if faint:
            return [f"0 (faint incomplete, {faint}% cells)", "Score 0, ultralow (<=10% faint)",
                    f"HER2 IHC 0 - faint/barely perceptible incomplete staining in {faint}%"][s]
        return ["0", "Score 0 (null)", "HER2 IHC 0 - no staining"][s]
    tag = {1: "", 2: " (equivocal)", 3: " (positive)"}[sc]
    return [f"{sc}+", f"Score {sc}+{tag}", f"HER2 IHC {sc}+"][s]


ihc_rows, ish_rows, var_rows, fus_rows, cnv_rows, sum_rows, ser_rows, clin_rows, truth_rows = ([] for _ in range(9))
rna_by_site = {s: {} for s in SITES}
CN_PANEL = ["ERBB2", "RB1", "PTEN", "MYCN", "AURKA", "CDKN2A", "MYC"]
STATUS_TO_GENE = {"TTF1": "NKX2-1", "NAPSA": "NAPSA", "P40": "TP63", "CK56": "KRT5", "GATA3": "GATA3", "TRPS1": "TRPS1",
                  "MGB": "SCGB2A2", "PSA": "KLK3", "NKX31": "NKX3-1", "PSMA": "FOLH1", "AMACR": "AMACR", "HOXB13": "HOXB13",
                  "SYP": "SYP", "CHGA": "CHGA", "INSM1": "INSM1", "ALK": "ALK"}


def fmt_vaf(v, site):
    return [f"{v:.3f}", f"{v * 100:.1f}%", f"{v * 100:.1f}"][si(site)]


def fmt_prot(prot, site):
    if prot == "exon14skip":
        return ["p.?", "p.?", "METex14 skipping"][si(site)]
    if si(site) == 1:
        return "p." + re.sub(r"[A-Z*]", lambda m: AA3[m.group(0)], prot[2:])
    if si(site) == 2:
        return prot[2:]
    return prot


def fmt_cls(c, site):
    return {"Pathogenic": ["Pathogenic", "P", "pathogenic"], "Likely pathogenic": ["Likely pathogenic", "LP", "likely_pathogenic"],
            "VUS": ["VUS", "VUS", "uncertain significance"]}[c][si(site)]


for P in patients:
    pid, site, T = P["pid"], P["site"], P["truth"]
    s = si(site)
    # ---------- clinical
    smoke = ""
    if P["pack_years"] is not None:
        py = P["pack_years"]
        smoke = (pick(["never", "Never smoker", "non-smoker"]) if py == 0 else
                 pick([f"{py} pack-years", f"former smoker, {py} py", f"current {py} PY"]))
    pur = P["purity"]
    pur_txt = "" if pur is None else [f"{pur * 100:.0f}%", f"{pur:.2f}", f"{pur * 100:.0f} %"][s]
    sex_txt = {"M": ["M", "Male", "male"], "F": ["F", "Female", "female"]}[P["sex"]][s]
    clin_rows.append(dict(patient_id=pid, site=site, sex=sex_txt, age=P["age"], specimen_site=P["specimen"],
                          smoking_history=smoke, clinical_history=P["history"], tumor_cellularity=pur_txt))

    # ---------- IHC (8% per-stain missingness)
    for m, v in P["ihc"].items():
        if v is None or p(0.08):
            continue
        ihc_rows.append(dict(patient_id=pid, site=site, stain=STAIN_NAMES[m][s], result=fmt_ihc_value(m, v, site, P)))
    if P["her2"] is not None:
        ihc_rows.append(dict(patient_id=pid, site=site, stain=STAIN_NAMES["HER2"][s], result=fmt_her2(P["her2"], site)))
    # ---------- ISH / FISH
    if P["her2_ish"] is not None:
        r, c = P["her2_ish"]
        txt = [f"HER2/CEP17 ratio={r:.2f}; HER2 copies/cell={c:.1f}", f"Ratio {r:.2f}, mean HER2 signals {c:.1f}",
               f"ERBB2:CEP17 {r:.2f} | avg ERBB2 {c:.1f}"][s]
        ish_rows.append(dict(patient_id=pid, site=site, test=["HER2 dual-probe ISH", "HER2 SISH", "ERBB2 FISH"][s], result_text=txt))
    if P["alk_fish"] is not None:
        a = P["alk_fish"]
        txt = [f"ALK break-apart: {a:.0f}% cells rearranged", f"ALK FISH (Vysis) positive, {a:.0f}% split signals",
               f"ALK BA FISH positive ({a:.0f}%)"][s]
        ish_rows.append(dict(patient_id=pid, site=site, test="ALK FISH", result_text=txt))
    elif P["organ"] == "lung" and p(0.3):
        a = U(0, 8)
        ish_rows.append(dict(patient_id=pid, site=site, test="ALK FISH",
                             result_text=[f"ALK break-apart: {a:.0f}% cells rearranged", f"ALK FISH (Vysis) negative, {a:.0f}% split signals",
                                          f"ALK BA FISH negative ({a:.0f}%)"][s]))

    # ---------- NGS DNA (95%)
    has_ngs = p(0.95)
    fusion_done = has_ngs and p(1.0 if P["organ"] == "lung" else 0.85)
    purity = pur if pur is not None else 0.0
    if has_ngs:
        variants = list(P["variants"])
        if p(0.15):
            g_, pr_ = CHIP[int(rng.integers(len(CHIP)))]
            variants.append(dict(var(g_, pr_), vaf=U(0.01, 0.035), chip=True))
        if p(0.2):
            g_, pr_ = VUS[int(rng.integers(len(VUS)))]
            variants.append(dict(var(g_, pr_, cls="VUS"), vaf=U(0.4, 0.55)))
        for v in variants:
            vaf = v["vaf"]
            if vaf is None:
                loh = v["gene"] in ("TP53", "RB1", "BRCA1", "BRCA2", "PTEN")
                vaf = purity * (U(0.75, 1.0) if loh else U(0.35, 0.6))
                vaf = max(vaf, 0.0)
            depth = int(U(150, 1200)) if not p(0.03) else int(U(40, 95))
            gene = v["gene"]
            var_rows.append(dict(patient_id=pid, site=site, gene=gene, hgvs_c=v["cdna"], hgvs_p=fmt_prot(v["prot"], site),
                                 variant_type=[v["vtype"], v["vtype"].upper(), v["vtype"]][s], vaf=fmt_vaf(min(vaf, 0.99), site),
                                 depth=depth, classification=fmt_cls(v["cls"], site)))
        for gene in CN_PANEL:
            cn = P["cn"].get(gene, 2)
            ratio = (purity * cn + 2 * (1 - purity)) / 2 * float(np.exp(rng.normal(0, 0.04)))
            metric, val = [("log2R", np.log2(ratio)), ("CN", 2 * ratio), ("FC", ratio)][s]
            gname = "HER2" if (gene == "ERBB2" and s == 2) else gene
            cnv_rows.append(dict(patient_id=pid, site=site, gene=gname, metric=metric, value=round(float(val), 3)))
        tmb, fga = P["tmb"], P["fga"]
        sum_rows.append(dict(patient_id=pid, site=site, tmb=[f"{tmb:.1f} mut/Mb", f"{tmb:.1f}", f"TMB: {tmb:.1f}"][s],
                             fga=[f"{fga * 100:.0f}%", f"{fga:.2f}", f"{fga * 100:.1f}"][s], msi=["MSS", "stable", "MSS"][s],
                             rna_fusion_panel=["performed", "yes", "Y"][s] if fusion_done else ["not performed", "no", "N"][s]))
    # ---------- RNA fusion panel (lung 100%, others 85%)
    if fusion_done:
        for five, three, reads in P["fusions"]:
            fus_rows.append(dict(patient_id=pid, site=site, fusion=[f"{five}::{three}", f"{five}-{three}", f"{five}/{three} (e{int(U(2, 20))}:e{int(U(2, 20))})"][s],
                                 supporting_reads=reads, assay="RNA fusion panel"))
    # ---------- serum (60%)
    if P["serum"] and p(0.6):
        for k, v in P["serum"].items():
            name = {"PSA": ["PSA", "Total PSA", "tPSA"], "CEA": ["CEA", "CEA", "Carcinoembryonic Ag"], "CA15_3": ["CA 15-3", "CA15-3", "CA15.3"],
                    "CYFRA21_1": ["CYFRA 21-1", "CYFRA21-1", "Cyfra"], "ProGRP": ["ProGRP", "Pro-GRP", "ProGRP"], "NSE": ["NSE", "NSE", "Neuron specific enolase"]}[k][s]
            unit_canon = {"PSA": "ng/mL", "CEA": "ng/mL", "CA15_3": "U/mL", "CYFRA21_1": "ng/mL", "ProGRP": "pg/mL", "NSE": "ng/mL"}[k]
            unit, factor = unit_canon, 1.0
            if s == 1:
                unit = {"ng/mL": "µg/L", "U/mL": "kU/L", "pg/mL": "ng/L"}[unit_canon]
            if s == 2 and k == "ProGRP":
                unit, factor = "ng/mL", 0.001
            val = v * factor
            txt = "<0.1" if (k == "PSA" and val < 0.1) else (f"{val:.3f}" if factor < 1 else f"{val:.2f}")
            ser_rows.append(dict(patient_id=pid, site=site, analyte=name, value=txt, unit=unit))
    # ---------- RNA expression (70%)
    if p(0.7):
        state = {}
        for m, v in P["ihc"].items():
            if m in STATUS_TO_GENE and isinstance(v, str):
                state[STATUS_TO_GENE[m]] = {"pos": "on", "weak": "mid", "neg": "off"}[v]
        for g_, st in P["rna"].items():
            state[g_] = st
        ihc = P["ihc"]
        if "ER" in ihc: state["ESR1"] = "on" if ihc["ER"] > 10 else "mid" if ihc["ER"] >= 1 else "off"
        if "PR" in ihc: state["PGR"] = "on" if ihc["PR"] > 10 else "mid" if ihc["PR"] >= 1 else "off"
        if "AR" in ihc: state["AR"] = "on" if ihc["AR"] >= 10 else "mid" if ihc["AR"] >= 1 else "off"
        if "PDL1" in ihc: state["CD274"] = "on" if ihc["PDL1"] >= 50 else "mid" if ihc["PDL1"] >= 1 else "off"
        ki = ihc.get("Ki67", 10)
        prolif = 5.5 + ki / 100 * 5
        lib = float(np.exp(rng.normal(0, 0.4))) * {"SITE_A": 1.0, "SITE_B": 2.2, "SITE_C": 0.5}[site]
        col = {}
        for g_, info in GENES.items():
            if info["hk"]:
                L = HK_LEVEL[g_] + rng.normal(0, 0.15)
            elif g_ in ("MKI67", "TOP2A", "BIRC5", "CCNB1"):
                L = prolif + rng.normal(0, 0.35)
            elif g_ == "ERBB2":
                h = P["her2"]["ihc"] if P["her2"] else 0
                L = (12.5 if state.get("ERBB2") == "amp" else {0: 6, 1: 7.2, 2: 8.2, 3: 12.5}[h]) + rng.normal(0, 0.4)
            else:
                st = state.get(g_, "off")
                if p(0.02):  # random discordance
                    st = {"on": "off", "off": "on", "mid": "off"}[st]
                L = {"on": rng.normal(10, 0.6), "mid": rng.normal(7.6, 0.4), "off": rng.normal(3, 0.8)}[st]
            if site == "SITE_C":
                L += SITE_GENE_BIAS[g_]
            lam = (2 ** L) * lib * info["length"] / 1000
            gid = info["alias"] if site == "SITE_C" else g_
            col[gid] = int(rng.poisson(lam))
        rna_by_site[site][pid] = col

    truth_rows.append(dict(patient_id=pid, scenario=P["scenario"], true_lineage=T["lineage"], true_histology=T["histology"],
                           true_molecular=T["molecular"], true_modifiers=";".join(sorted(set(T["modifiers"]))),
                           edge_case=T["edge"], note=T["note"].strip("; ")))

pd.DataFrame(clin_rows).to_csv(RAW / "clinical.csv", index=False)
pd.DataFrame(ihc_rows).to_csv(RAW / "ihc.csv", index=False)
pd.DataFrame(ish_rows).to_csv(RAW / "ish.csv", index=False)
pd.DataFrame(var_rows).to_csv(RAW / "ngs_variants.csv", index=False)
pd.DataFrame(fus_rows).to_csv(RAW / "ngs_fusions.csv", index=False)
pd.DataFrame(cnv_rows).to_csv(RAW / "ngs_cnv.csv", index=False)
pd.DataFrame(sum_rows).to_csv(RAW / "ngs_summary.csv", index=False)
pd.DataFrame(ser_rows).to_csv(RAW / "serum.csv", index=False)
for s_, cols in rna_by_site.items():
    pd.DataFrame(cols).rename_axis("gene_id").to_csv(RAW / f"rnaseq_counts_{s_}.csv")
pd.DataFrame(truth_rows).to_csv(ROOT / "data" / "ground_truth.csv", index=False)

print(f"patients={len(patients)}")
for f in sorted(RAW.glob("*.csv")):
    print(f"  {f.name:30s} {sum(1 for _ in open(f)) - 1:6d} rows")
