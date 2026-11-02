"""
Normalize raw multi-site molecular data into canonical, ontology-aligned findings.

Every vocabulary mapping, unit conversion, reference limit and threshold is read
from ontology/oncokg.ttl. Steps:

  clinical  : sex, age, pack-years, specimen site, history flags, tumour purity -> fraction
  IHC       : stain synonyms -> HGNC biomarker; %, Allred PS/IS, fraction, free text -> percent
              + guideline status (ER neg/low/pos, PR, AR, PD-L1 TPS category, HER2 IHC score
              incl. ultralow/null); qualitative stains -> Positive / Weak / Negative
  ISH/FISH  : HER2/CEP17 ratio + copies/cell -> ASCO/CAP ISH group; ALK % split cells
  HER2      : integrated ASCO/CAP 2018 (affirmed 2023) + NCCN 2026 null/ultralow/low status
  variants  : gene aliases, 3-letter -> 1-letter HGVS, textual exon-14 skipping, VAF scale
              inference per site, ACMG/AMP class harmonisation, depth/VAF QC, alteration matching
  fusions   : separator harmonisation, read-support QC, fusion alteration matching
  CNV       : log2R / CN / fold-change -> ratio -> purity-adjusted copy number -> amp / homdel
  summary   : TMB (mut/Mb), fraction genome altered (scale inference)
  serum     : analyte synonyms, unit conversion to canonical unit, censored values, x ULN
  RNA       : gene-ID harmonisation, length normalisation, housekeeping-centred log2 (nlog2),
              per-gene high/low calls, proliferation score
  consolidation : per-marker status from IHC (preferred) or RNA, with discordance flags

Outputs (output/): normalized_findings.csv (long, with provenance),
                   normalized_features.csv + features.json (wide, one row per patient),
                   normalization_report.json
"""
import re
import json
import math
from collections import defaultdict
import numpy as np
import pandas as pd
from common import ROOT, OKG, q, params, local

RAW = ROOT / "data" / "raw"
OUT = ROOT / "output"
OUT.mkdir(exist_ok=True)
PRM = params()
report = defaultdict(int)
findings = []          # long-format provenance rows
F = defaultdict(dict)  # wide features per patient
FLAGS = defaultdict(list)


def add(pid, feature, value, unit="", source="", raw="", method=""):
    F[pid][feature] = value
    findings.append(dict(patient_id=pid, feature=feature, value=value, unit=unit, source=source,
                         raw_value=raw, method=method))


# ------------------------------------------------------------------ ontology lookups
IHC_NAME = {}
for name, sym in q("SELECT ?n ?s WHERE { ?b a okg:ProteinExpressionBiomarker ; okg:ihcName ?n ; okg:hgncSymbol ?s }"):
    IHC_NAME[str(name).strip().lower()] = str(sym)

GENE_ALIAS = {}
GENE_INFO = {}
for gnode, sym in q("SELECT ?g ?s WHERE { ?g a okg:Gene ; okg:hgncSymbol ?s }"):
    GENE_ALIAS[str(sym).upper()] = str(sym)
for a, sym in q("SELECT ?a ?s WHERE { ?g a okg:Gene ; okg:hgncSymbol ?s ; skos:altLabel ?a }"):
    GENE_ALIAS[str(a).upper()] = str(sym)
for sym, ln, thr in q("""SELECT ?s ?l ?t WHERE { ?g a okg:Gene ; okg:hgncSymbol ?s ; okg:effectiveLengthBp ?l
                         OPTIONAL { ?g okg:rnaHighThreshold ?t } }"""):
    GENE_INFO[str(sym)] = dict(length=int(ln), thr=float(thr) if thr is not None else None, hk=False)
for (sym,) in q("SELECT ?s WHERE { ?g a okg:HousekeepingGene ; okg:hgncSymbol ?s }"):
    GENE_INFO[str(sym)]["hk"] = True

SERUM = {}
for b, key, uln, cu in q("""SELECT ?b ?k ?u ?cu WHERE { ?b a okg:SerumBiomarker ; okg:featureKey ?k ;
                           okg:upperLimitNormal ?u ; okg:canonicalUnit ?c . ?c okg:ucumCode ?cu }"""):
    SERUM[str(key).split(".")[1]] = dict(uln=float(uln), unit=str(cu))
SERUM_SYN = {str(a).strip().lower(): str(k).split(".")[1] for a, k in
             q("SELECT ?a ?k WHERE { ?b a okg:SerumBiomarker ; skos:altLabel ?a ; okg:featureKey ?k }")}
UNIT_SYN = {str(a): local(u) for u, a in q("SELECT ?u ?a WHERE { ?u a okg:Unit ; skos:altLabel ?a }")}
CONV = {(str(k), local(fu)): float(f) for k, fu, f in
        q("SELECT ?k ?fu ?f WHERE { ?c a okg:UnitConversion ; okg:analyteKey ?k ; okg:fromUnit ?fu ; okg:conversionFactor ?f }")}

SEQ_ALTS, FUS_ALTS, CNA_ALTS = [], [], []
for a, key, gene, kind, pathreq in q("""SELECT ?a ?k ?g ?t ?pr WHERE { ?a okg:featureKey ?k ; okg:matchGene ?g ; a ?t .
        FILTER(?t IN (okg:SequenceVariantAlteration, okg:GeneFusionAlteration, okg:CopyNumberAlteration))
        OPTIONAL { ?a okg:requiresPathogenic ?pr } }"""):
    G = None
    from common import ontology
    G = ontology()
    rec = dict(key=str(key), gene=str(gene), req_path=bool(pathreq.toPython()) if pathreq is not None else False,
               prot=[re.compile(str(x)) for x in G.objects(a, OKG.matchProteinRegex)],
               cdna=[re.compile(str(x)) for x in G.objects(a, OKG.matchCdnaRegex)])
    t = local(kind)
    if t == "SequenceVariantAlteration":
        SEQ_ALTS.append(rec)
    elif t == "GeneFusionAlteration":
        rec["three"] = str(G.value(a, OKG.fusion3Prime))
        five = G.value(a, OKG.fusion5Prime)
        rec["five"] = str(five) if five is not None else None
        FUS_ALTS.append(rec)
    else:
        rec["cna"] = str(G.value(a, OKG.cnaType))
        CNA_ALTS.append(rec)
ALL_ALT_KEYS = sorted({r["key"] for r in SEQ_ALTS + FUS_ALTS + CNA_ALTS})

AA1 = {"Ala": "A", "Arg": "R", "Asn": "N", "Asp": "D", "Cys": "C", "Glu": "E", "Gln": "Q", "Gly": "G", "His": "H",
       "Ile": "I", "Leu": "L", "Lys": "K", "Met": "M", "Phe": "F", "Pro": "P", "Ser": "S", "Thr": "T", "Trp": "W",
       "Tyr": "Y", "Val": "V", "Ter": "*"}


def num(s):
    m = re.search(r"-?\d+(?:\.\d+)?", str(s))
    return float(m.group(0)) if m else None


# ------------------------------------------------------------------ clinical
clin = pd.read_csv(RAW / "clinical.csv", dtype=str).fillna("")
PIDS = list(clin.patient_id)
for r in clin.itertuples():
    pid = r.patient_id
    add(pid, "site", r.site, source="clinical")
    sex = "M" if r.sex.strip().lower() in ("m", "male") else "F" if r.sex.strip().lower() in ("f", "female") else None
    add(pid, "clin.sex", sex, source="clinical", raw=r.sex, method="sex vocabulary")
    add(pid, "clin.age", int(r.age), "years", "clinical", r.age)
    sm = r.smoking_history.lower()
    py = 0.0 if re.search(r"never|non-smoker", sm) else num(sm) if sm else None
    if py is not None:
        add(pid, "clin.pack_years", py, "pack-years", "clinical", r.smoking_history, "free text -> pack-years")
    add(pid, "specimen_site", r.specimen_site.lower(), source="clinical", raw=r.specimen_site)
    hist = r.clinical_history.lower()
    add(pid, "history.prostate_cancer", bool(re.search(r"prostat", hist)), source="clinical", raw=r.clinical_history, method="keyword")
    add(pid, "history.breast_cancer", bool(re.search(r"breast", hist)), source="clinical", raw=r.clinical_history, method="keyword")
    if r.tumor_cellularity:
        v = num(r.tumor_cellularity)
        frac = v / 100 if ("%" in r.tumor_cellularity or v > 1) else v
        add(pid, "purity", round(frac, 3), "fraction", "clinical", r.tumor_cellularity, "% or fraction -> fraction")
        report["purity_converted_from_percent"] += int("%" in r.tumor_cellularity or v > 1)
        if frac < PRM["LOW_PURITY_FLAG"]:
            FLAGS[pid].append(f"Low tumour cellularity ({frac:.0%}): sensitivity for variants/CNAs reduced")

# ------------------------------------------------------------------ IHC
ALLRED_RANGE = {0: (0, 0), 1: (0, 0.99), 2: (1, 10), 3: (11, 33), 4: (34, 66), 5: (67, 100)}


def parse_percent_text(t):
    """returns (pct, lo, hi, method)"""
    tl = t.lower()
    m = re.search(r"allred\s*ps\s*(\d)", tl)
    if m:
        lo, hi = ALLRED_RANGE[int(m.group(1))]
        return (lo + hi) / 2, lo, hi, f"Allred proportion score {m.group(1)} -> {lo}-{hi}%"
    m = re.search(r"(<)\s*(\d+(?:\.\d+)?)\s*%", tl)
    if m:
        return float(m.group(2)) / 2, 0.0, float(m.group(2)) - 0.01, "'<x%' -> below threshold"
    m = re.search(r"(\d+(?:\.\d+)?)\s*%", tl)
    if m:
        v = float(m.group(1)); return v, v, v, "percent"
    m = re.search(r"index\s*(\d*\.\d+)", tl)
    if m:
        v = float(m.group(1)) * 100; return v, v, v, "fraction index -> percent"
    v = num(tl)
    if v is not None:
        if v <= 1 and "." in tl:
            return v * 100, v * 100, v * 100, "fraction -> percent"
        return v, v, v, "bare number -> percent"
    if re.search(r"neg|no staining", tl):
        return 0.0, 0.0, 0.0, "qualitative negative -> 0%"
    return None, None, None, "unparsed"


def parse_status(t):
    tl = t.lower().strip()
    if re.search(r"\bneg|no staining|^-$|not detected", tl):
        return "Negative"
    if re.search(r"weak|focal|equivocal", tl):
        return "Weak"
    if re.search(r"pos|strong|diffuse|\+", tl):
        return "Positive"
    return None


def parse_her2_ihc(t):
    tl = t.lower()
    m = re.search(r"([0-3])\s*\+", tl)
    score = int(m.group(1)) if m else (0 if re.search(r"(score|ihc)?\s*\b0\b", tl) else None)
    sub = None
    if score == 0:
        sub = "Ultralow" if re.search(r"ultralow|faint|incomplete|barely", tl) else "Null"
    return score, sub


ihc = pd.read_csv(RAW / "ihc.csv", dtype=str).fillna("")
IHC_STATUS = defaultdict(dict)
unmapped = set()
for r in ihc.itertuples():
    pid, txt = r.patient_id, r.result
    sym = IHC_NAME.get(r.stain.strip().lower())
    if sym is None:
        unmapped.add(r.stain); continue
    report["ihc_stain_synonyms_mapped"] += int(r.stain != sym)
    src = f"IHC:{r.stain}"
    if sym in ("ESR1", "PGR", "AR"):
        pct, lo, hi, meth = parse_percent_text(txt)
        if pct is None:
            continue
        key = {"ESR1": "ER", "PGR": "PR", "AR": "AR"}[sym]
        add(pid, f"{key}_pct", round(pct, 1), "%", src, txt, meth)
        F[pid][f"{key}_pct__lo"], F[pid][f"{key}_pct__hi"] = lo, hi
        report["allred_converted" if "Allred" in meth else "receptor_percent_parsed"] += 1
        if key == "ER":
            st = "Negative" if hi < PRM["ER_LOW_MIN_PCT"] else "Low" if hi <= PRM["ER_LOW_MAX_PCT"] else "Positive"
            add(pid, "ER_status", st, "", src, txt, "ASCO/CAP 2020: <1% neg, 1-10% low, >10% pos")
            IHC_STATUS[pid]["ESR1"] = {"Positive": "Positive", "Low": "Weak", "Negative": "Negative"}[st]
        elif key == "PR":
            st = "Negative" if hi < PRM["PR_POS_MIN_PCT"] else "Positive"
            add(pid, "PR_status", st, "", src, txt, "ASCO/CAP 2020: >=1% positive")
            IHC_STATUS[pid]["PGR"] = st
        else:
            IHC_STATUS[pid]["AR"] = "Positive" if pct >= PRM["AR_POS_MIN_PCT"] else "Weak" if pct >= 1 else "Negative"
    elif sym == "MKI67":
        pct, _, _, meth = parse_percent_text(txt)
        add(pid, "Ki67_pct", round(pct, 1), "%", src, txt, meth)
        report["ki67_fraction_converted"] += int("fraction" in meth)
    elif sym == "CD274":
        pct, lo, hi, meth = parse_percent_text(txt)
        add(pid, "PDL1_TPS", round(pct, 1), "%", src, txt, meth)
        st = "High" if lo >= PRM["PDL1_HIGH_TPS"] else "Low" if lo >= PRM["PDL1_POS_TPS"] else "Negative"
        add(pid, "PDL1_status", st, "", src, txt, "TPS <1 / 1-49 / >=50")
    elif sym == "ERBB2":
        score, sub = parse_her2_ihc(txt)
        add(pid, "HER2_ihc", score, "score", src, txt, "IHC score parse")
        if sub:
            add(pid, "HER2_ihc0_subtype", sub, "", src, txt, "NCCN 2026 null vs ultralow")
    else:
        st = parse_status(txt)
        if st:
            IHC_STATUS[pid][sym] = st
            add(pid, f"ihc.{sym}", st, "", src, txt, "qualitative stain vocabulary")
report["ihc_unmapped_stains"] = len(unmapped)

# ------------------------------------------------------------------ ISH / FISH
ish = pd.read_csv(RAW / "ish.csv", dtype=str).fillna("")
for r in ish.itertuples():
    pid, t = r.patient_id, r.result_text
    if "alk" in r.test.lower():
        pct = num(re.search(r"(\d+(?:\.\d+)?)\s*%", t).group(1))
        add(pid, "ALK_FISH_pct", pct, "% cells", f"ISH:{r.test}", t, "% split signals")
        continue
    ratio = float(re.search(r"(?:ratio|CEP17)\s*[=:]?\s*(\d+(?:\.\d+)?)", t, re.I).group(1))
    copies = float(re.search(r"(?:copies/cell|signals|avg ERBB2)\s*[=:]?\s*(\d+(?:\.\d+)?)", t, re.I).group(1))
    rc, cl, ch = PRM["HER2_ISH_RATIO_CUT"], PRM["HER2_ISH_COPY_LOW"], PRM["HER2_ISH_COPY_HIGH"]
    grp = (1 if ratio >= rc and copies >= cl else 2 if ratio >= rc else 3 if copies >= ch else 4 if copies >= cl else 5)
    add(pid, "HER2_ish_ratio", ratio, "ratio", f"ISH:{r.test}", t, "regex")
    add(pid, "HER2_ish_copies", copies, "copies/cell", f"ISH:{r.test}", t, "regex")
    add(pid, "HER2_ish_group", grp, "", f"ISH:{r.test}", t, "ASCO/CAP 2018 ISH groups 1-5")

# HER2 integrated status
for pid in PIDS:
    f = F[pid]
    s, grp = f.get("HER2_ihc"), f.get("HER2_ish_group")
    st, why = None, ""
    if s == 3:
        st, why = "Positive", "IHC 3+"
    elif s == 2:
        if grp is None:
            st, why = "Equivocal", "IHC 2+ without ISH"
            FLAGS[pid].append("HER2 IHC 2+ (equivocal) - reflex ISH required before subtype assignment")
        elif grp in (1, 3):
            st, why = "Positive", f"IHC 2+ / ISH group {grp}"
        else:
            st, why = "Low", f"IHC 2+ / ISH group {grp} (not amplified)"
    elif s == 1:
        st, why = "Low", "IHC 1+"
    elif s == 0:
        st, why = f.get("HER2_ihc0_subtype", "Null"), "IHC 0 (" + f.get("HER2_ihc0_subtype", "Null").lower() + ")"
    elif grp is not None:
        st, why = ("Positive" if grp == 1 else "Negative"), f"ISH only, group {grp}"
    if st:
        add(pid, "HER2_status", st, "", "IHC+ISH", why, "ASCO/CAP 2018 (affirmed 2023) + NCCN 2026")

# ------------------------------------------------------------------ variants
var = pd.read_csv(RAW / "ngs_variants.csv", dtype=str).fillna("")
SITE_VAF_PERCENT = {}
for site, d in var.groupby("site"):
    vals = [num(v) for v in d.vaf]
    SITE_VAF_PERCENT[site] = any("%" in v for v in d.vaf) or max(vals) > 1.0
report["vaf_scale_by_site"] = {k: ("percent" if v else "fraction") for k, v in SITE_VAF_PERCENT.items()}


def norm_protein(pr):
    s = pr.strip()
    if re.search(r"ex(on)?\s*14", s, re.I):
        return "exon14skip"
    if s in ("", "p.?"):
        return s
    if not s.startswith("p."):
        s = "p." + s
    s2 = re.sub("|".join(AA1), lambda m: AA1[m.group(0)], s)
    return s2.replace("X", "*")


def norm_class(c):
    c = c.strip().lower()
    if c in ("p", "pathogenic", "path"):
        return "Pathogenic"
    if c in ("lp", "likely pathogenic", "likely_pathogenic"):
        return "Likely pathogenic"
    if "uncertain" in c or c == "vus":
        return "VUS"
    if "benign" in c or c in ("b", "lb"):
        return "Benign"
    return "Unknown"


NGS_DONE = set(pd.read_csv(RAW / "ngs_summary.csv", dtype=str).patient_id)
ALT_HITS = defaultdict(set)
DRIVERS = defaultdict(int)
for r in var.itertuples():
    pid = r.patient_id
    gene = GENE_ALIAS.get(r.gene.upper(), r.gene.upper())
    prot = norm_protein(r.hgvs_p)
    report["hgvs_3letter_converted"] += int(bool(re.search("|".join(AA1), r.hgvs_p)))
    report["hgvs_prefix_added"] += int(bool(r.hgvs_p) and not r.hgvs_p.startswith("p."))
    v = num(r.vaf)
    vaf = v / 100 if SITE_VAF_PERCENT[r.site] else v
    depth = int(r.depth)
    cls_ = norm_class(r.classification)
    vt = r.variant_type.lower()
    qc = "PASS"
    if depth < PRM["MIN_VARIANT_DEPTH"]:
        qc = "LOW_DEPTH"
    elif vaf < PRM["MIN_VAF"]:
        qc = "LOW_VAF"
    report[f"variant_qc_{qc}"] += 1
    hits = []
    if qc == "PASS":
        for a in SEQ_ALTS:
            if a["gene"] != gene:
                continue
            if a["prot"] or a["cdna"]:
                ok = any(rx.search(prot) for rx in a["prot"]) or any(rx.search(r.hgvs_c) for rx in a["cdna"])
            else:
                ok = cls_ in ("Pathogenic", "Likely pathogenic") if a["req_path"] else True
            if ok:
                hits.append(a["key"]); ALT_HITS[pid].add(a["key"])
        if hits or cls_ in ("Pathogenic", "Likely pathogenic"):
            DRIVERS[pid] += 1
    elif cls_ in ("Pathogenic", "Likely pathogenic") and qc == "LOW_DEPTH":
        FLAGS[pid].append(f"{gene} {prot} reported P/LP but failed depth QC ({depth}x) - confirm")
    findings.append(dict(patient_id=pid, feature=f"variant.{gene}", value=f"{prot} VAF={vaf:.3f} {cls_} [{qc}]",
                         unit="", source=f"NGS:{r.site}", raw_value=f"{r.gene} {r.hgvs_p} {r.hgvs_c} vaf={r.vaf} {r.classification}",
                         method="alias->HGNC; 3->1 letter HGVS; VAF scale; ACMG class; QC; matched=" + ",".join(hits)))

# ------------------------------------------------------------------ fusions
fus = pd.read_csv(RAW / "ngs_fusions.csv", dtype=str).fillna("")
FUSION_DONE = set()
clin_organ_guess = None
FUSIONS = defaultdict(int)
for r in fus.itertuples():
    pid = r.patient_id
    core = re.sub(r"\(.*?\)", "", r.fusion).strip()
    parts = [x for x in re.split(r"::|--|/|-", core) if x]
    five, three = (GENE_ALIAS.get(parts[0].upper(), parts[0].upper()), GENE_ALIAS.get(parts[1].upper(), parts[1].upper()))
    reads = int(r.supporting_reads)
    qc = "PASS" if reads >= PRM["MIN_FUSION_READS"] else "LOW_SUPPORT"
    report[f"fusion_qc_{qc}"] += 1
    hits = []
    if qc == "PASS":
        for a in FUS_ALTS:
            if a["three"] == three and (a["five"] is None or a["five"] == five):
                hits.append(a["key"]); ALT_HITS[pid].add(a["key"])
        if hits:
            FUSIONS[pid] += 1
    findings.append(dict(patient_id=pid, feature="fusion", value=f"{five}::{three} reads={reads} [{qc}]", unit="",
                         source="RNA fusion panel", raw_value=r.fusion, method="separator harmonisation; matched=" + ",".join(hits)))
# the fusion export lists detected fusions only; whether the panel was run comes from the NGS summary manifest
_s = pd.read_csv(RAW / "ngs_summary.csv", dtype=str).fillna("")
FUSION_DONE = set(_s.patient_id[_s.rna_fusion_panel.str.lower().isin(["performed", "yes", "y"])])
FUS_KEYS = {a["key"] for a in FUS_ALTS}

# ------------------------------------------------------------------ CNV
cnv = pd.read_csv(RAW / "ngs_cnv.csv", dtype=str)
CNAS = defaultdict(int)
for r in cnv.itertuples():
    pid = r.patient_id
    gene = GENE_ALIAS.get(r.gene.upper(), r.gene.upper())
    val = float(r.value)
    ratio = {"log2R": 2 ** val, "CN": val / 2, "FC": val}[r.metric]
    pur = F[pid].get("purity", 1.0)
    if pur < PRM["MIN_PURITY_FOR_CN"]:
        report["cn_not_called_low_purity"] += 1
        continue
    cn = (2 * ratio - 2 * (1 - pur)) / pur
    cn = max(cn, 0.0)
    call = "amp" if cn >= PRM["CN_AMP_MIN"] else "homdel" if cn <= PRM["CN_HOMDEL_MAX"] else "neutral"
    report[f"cnv_metric_{r.metric}"] += 1
    if call != "neutral":
        for a in CNA_ALTS:
            if a["gene"] == gene and a["cna"] == call:
                ALT_HITS[pid].add(a["key"])
        CNAS[pid] += 1
    if call != "neutral" or gene == "ERBB2":
        add(pid, f"cn.{gene}", round(cn, 2), "copies (purity-adjusted)", f"CNV:{r.site}", f"{r.metric}={r.value}",
            f"{r.metric} -> ratio -> (2r - 2(1-p))/p, p={pur:.2f}; call={call}")

# ------------------------------------------------------------------ NGS summary
summ = pd.read_csv(RAW / "ngs_summary.csv", dtype=str)
SITE_FGA_PERCENT = {s: any("%" in v for v in d.fga) or max(num(v) for v in d.fga) > 1 for s, d in summ.groupby("site")}
for r in summ.itertuples():
    add(r.patient_id, "tmb", num(r.tmb), "mut/Mb", "NGS summary", r.tmb, "number extraction")
    fv = num(r.fga)
    add(r.patient_id, "fga", round(fv / 100 if SITE_FGA_PERCENT[r.site] else fv, 3), "fraction", "NGS summary", r.fga,
        "site scale inference (% vs fraction)")

for pid in PIDS:
    if pid in NGS_DONE:
        for k in ALL_ALT_KEYS:
            if k in FUS_KEYS and pid not in FUSION_DONE:
                continue  # fusion panel not run -> unknown, not negative
            F[pid][k] = k in ALT_HITS[pid]
        h = ALT_HITS[pid]
        F[pid]["alt.RB1_ANY"] = bool({"alt.RB1_TRUNC", "alt.RB1_HOMDEL"} & h)
        F[pid]["alt.PTEN_ANY"] = bool({"alt.PTEN_TRUNC", "alt.PTEN_HOMDEL"} & h)
        if pid in FUSION_DONE:
            F[pid]["alt.NTRK_ANY"] = bool({"alt.NTRK1_FUSION", "alt.NTRK2_FUSION", "alt.NTRK3_FUSION"} & h)
        F[pid]["alt.HRR_ANY"] = bool({f"alt.{g}_PATHOGENIC" for g in ("BRCA1", "BRCA2", "ATM", "CDK12", "PALB2")} & h)
        F[pid]["n_somatic_drivers"] = DRIVERS[pid]
        if pid in FUSION_DONE:
            F[pid]["n_driver_fusions"] = FUSIONS[pid]
        F[pid]["n_focal_cna"] = CNAS[pid]
        for k in sorted(h):
            findings.append(dict(patient_id=pid, feature=k, value=True, unit="", source="NGS", raw_value="", method="ontology alteration match"))
    # ALK integrated
    f = F[pid]
    ev = []
    if f.get("alt.ALK_FUSION"): ev.append("RNA fusion")
    if f.get("ALK_FISH_pct", 0) >= PRM["ALK_FISH_POS_PCT"]: ev.append("FISH >=15%")
    if IHC_STATUS[pid].get("ALK") == "Positive": ev.append("D5F3 IHC")
    tested = pid in FUSION_DONE or "ALK_FISH_pct" in f or "ALK" in IHC_STATUS[pid]
    if tested:
        add(pid, "ALK_rearranged", bool(ev), "", "NGS/FISH/IHC", "; ".join(ev), "any of NGS fusion, FISH>=15%, IHC D5F3")
        if IHC_STATUS[pid].get("ALK") == "Positive" and "ALK_FISH_pct" in f and f["ALK_FISH_pct"] < PRM["ALK_FISH_POS_PCT"]:
            FLAGS[pid].append("ALK IHC positive but FISH negative - discordant, orthogonal confirmation advised")

# ------------------------------------------------------------------ serum
ser = pd.read_csv(RAW / "serum.csv", dtype=str)
for r in ser.itertuples():
    key = SERUM_SYN.get(r.analyte.strip().lower())
    unit = UNIT_SYN.get(r.unit.strip())
    if key is None or unit is None:
        report["serum_unmapped"] += 1; continue
    factor = CONV.get((key, unit))
    censored = r.value.strip().startswith("<")
    v = num(r.value)
    if censored:
        v = v / 2
    val = v * factor
    report["serum_unit_converted"] += int(factor != 1.0 or r.unit.strip() != SERUM[key]["unit"])
    add(r.patient_id, f"serum.{key}", round(val, 3), SERUM[key]["unit"], f"serum:{r.site}", f"{r.value} {r.unit}",
        f"x{factor} -> {SERUM[key]['unit']}" + ("; left-censored -> LOD/2" if censored else ""))
    add(r.patient_id, f"serum.{key}_xULN", round(val / SERUM[key]["uln"], 2), "x ULN", f"serum:{r.site}", "",
        f"ULN {SERUM[key]['uln']} {SERUM[key]['unit']}")

# ------------------------------------------------------------------ RNA
RNA_CALL = defaultdict(dict)
HK = [g for g, i in GENE_INFO.items() if i["hk"]]
PROLIF = ["MKI67", "TOP2A", "BIRC5", "CCNB1"]
rna_qc = []
for site_file in sorted(RAW.glob("rnaseq_counts_*.csv")):
    m = pd.read_csv(site_file, index_col=0)
    mapped = [GENE_ALIAS.get(str(i).upper(), str(i)) for i in m.index]
    report["rna_gene_ids_remapped"] += sum(a != b for a, b in zip(mapped, m.index))
    m.index = mapped
    lengths = np.array([GENE_INFO[g]["length"] for g in m.index]) / 1000.0
    rpk = m.div(lengths, axis=0)
    l2 = np.log2(rpk + 1)
    hk_mean = l2.loc[HK].mean(axis=0)
    nl2 = l2 - hk_mean + PRM["RNA_NORM_OFFSET"]
    for pid in m.columns:
        lib = int(m[pid].sum())
        rna_qc.append(dict(patient_id=pid, site=site_file.stem.split("_", 2)[-1], library_counts=lib,
                           hk_log2_mean=round(float(hk_mean[pid]), 2)))
        for gsym in m.index:
            if GENE_INFO[gsym]["hk"]:
                continue
            val = round(float(nl2.at[gsym, pid]), 2)
            thr = GENE_INFO[gsym]["thr"]
            F[pid][f"rna.{gsym}"] = val
            if thr is not None and gsym not in PROLIF:
                RNA_CALL[pid][gsym] = "High" if val >= thr else "Low"
        ps = round(float(nl2.loc[PROLIF, pid].mean()), 2)
        add(pid, "prolif_score", ps, "nlog2", "RNA", "", "mean nlog2 of MKI67/TOP2A/BIRC5/CCNB1")
        add(pid, "prolif_high", ps >= PRM["PROLIF_HIGH_NLOG2"], "", "RNA", "", f">= {PRM['PROLIF_HIGH_NLOG2']}")
        findings.append(dict(patient_id=pid, feature="rna.panel", value=f"{len(m.index) - len(HK)} genes", unit="nlog2",
                             source=f"RNA:{site_file.stem}", raw_value=f"library={lib}",
                             method="counts / (length kb) -> log2(+1) - mean(HK) + 10"))
pd.DataFrame(rna_qc).to_csv(OUT / "rna_qc.csv", index=False)

# ------------------------------------------------------------------ consolidation
MARKERS = [str(s) for (s,) in q("SELECT ?s WHERE { ?b a okg:ProteinExpressionBiomarker ; okg:hgncSymbol ?s }")]
NE = ["SYP", "CHGA", "INSM1", "ASCL1"]
for pid in PIDS:
    ih, rc = IHC_STATUS[pid], RNA_CALL[pid]
    for sym in MARKERS:
        if sym in ("MKI67", "CD274", "ERBB2", "ALK") and sym not in ih:
            continue
        i_st, r_st = ih.get(sym), rc.get(sym)
        if i_st:
            st, src = i_st, "IHC"
            if r_st and ((i_st == "Positive" and r_st == "Low") or (i_st == "Negative" and r_st == "High")):
                FLAGS[pid].append(f"{sym}: IHC {i_st} vs RNA {r_st} - discordant (IHC used)")
                report["ihc_rna_discordance"] += 1
        elif r_st:
            thr = GENE_INFO[sym]["thr"] if sym in GENE_INFO else None
            v = F[pid].get(f"rna.{sym}")
            st = "Negative" if r_st == "Low" else ("Positive" if v >= thr + PRM["RNA_STRONG_DELTA"] else "Weak")
            src = "RNA"
            report["status_from_rna_only"] += 1
        else:
            continue
        add(pid, f"status.{sym}", st, "", src, f"IHC={i_st} RNA={r_st}", "IHC preferred; RNA fallback")
    # ER / PR fallback to RNA when IHC is missing (flagged; ER-low cannot be inferred from RNA)
    for key, sym in (("ER", "ESR1"), ("PR", "PGR")):
        if f"{key}_status" not in F[pid] and f"rna.{sym}" in F[pid] and F[pid].get("HER2_status"):
            st = F[pid].get(f"status.{sym}")
            if st in ("Positive", "Negative"):
                add(pid, f"{key}_status", st, "", "RNA", f"{sym} nlog2={F[pid][f'rna.{sym}']}", "IHC missing -> RNA fallback")
                FLAGS[pid].append(f"{key} status inferred from {sym} RNA (IHC missing) - confirm by IHC")
                report["receptor_status_from_rna"] += 1
            else:
                FLAGS[pid].append(f"{key} IHC missing and {sym} RNA intermediate - cannot classify (possible {key}-low)")
    known = [F[pid].get(f"status.{s}") for s in NE]
    if any(k is not None for k in known):
        add(pid, "NE_count", sum(k in ("Positive", "Weak") for k in known), "markers", "IHC/RNA", "",
            "count of SYP/CHGA/INSM1/ASCL1 positive")
    if "ER_status" in F[pid] and "PR_status" in F[pid] and F[pid]["ER_status"] == "Negative" and F[pid]["PR_status"] == "Positive":
        FLAGS[pid].append("ER-negative / PR-positive: rare discordant phenotype - repeat ER/PR testing advised")
    F[pid]["flags"] = FLAGS[pid]

# ------------------------------------------------------------------ write
fd = pd.DataFrame(findings)
fd.to_csv(OUT / "normalized_findings.csv", index=False)
wide = pd.DataFrame.from_dict({p: {k: v for k, v in F[p].items() if k != "flags" and not k.startswith("rna.")} for p in PIDS},
                              orient="index").rename_axis("patient_id")
wide.to_csv(OUT / "normalized_features.csv")
json.dump({p: F[p] for p in PIDS}, open(OUT / "features.json", "w"), default=str)
report = {k: (dict(v) if isinstance(v, dict) else v) for k, v in report.items()}
json.dump(report, open(OUT / "normalization_report.json", "w"), indent=1)
print(f"patients={len(PIDS)} findings={len(fd)} wide_cols={wide.shape[1]} unmapped_stains={sorted(unmapped)}")
print(json.dumps(report, indent=1))
