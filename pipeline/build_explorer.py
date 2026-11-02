"""Bundle ontology, normalization provenance, reasoning traces and evaluation into one self-contained HTML explorer."""
import json
from collections import defaultdict
import pandas as pd
from common import ROOT, OKG, q, local, label, ontology
from reason import CLASS_INFO, PARENT, GROUPS, MODS, LINEAGES

OUT = ROOT / "output"
RAW = ROOT / "data" / "raw"
O = ontology()
ev = json.load(open(OUT / "evaluation.json"))
rep = json.load(open(OUT / "normalization_report.json"))
traces = json.load(open(OUT / "evidence_trace.json"))
dx = pd.read_csv(OUT / "patient_diagnoses.csv", dtype=str).fillna("")
det = pd.read_csv(OUT / "evaluation_detail.csv", dtype=str, keep_default_na=False).set_index("patient_id")
fnd = pd.read_csv(OUT / "normalized_findings.csv", dtype=str).fillna("")
clin = pd.read_csv(RAW / "clinical.csv", dtype=str).fillna("").set_index("patient_id")

# ---------------- ontology views
def anc(c):
    out, st = {c}, [c]
    while st:
        for p in PARENT[st.pop()]:
            if p not in out: out.add(p); st.append(p)
    return out
spec_counts_true = det.true_molecular.where(det.true_molecular != "", det.true_histology).value_counts().to_dict()
spec_counts_pred = det.molecular.where(det.molecular != "", det.histology).value_counts().to_dict()
classes = []
for c, info in CLASS_INFO.items():
    if "Disease" not in anc(c) or c == "Disease":
        continue
    par = [p for p in PARENT[c]]
    classes.append(dict(id=c, label=info["label"], parent=par[0] if par else None, rare=info["rare"], tier=info["tier"],
                        prev=info["prevalence"], therapy=info["therapy"], nt=int(spec_counts_true.get(c, 0)),
                        np=int(spec_counts_pred.get(c, 0))))

lin_rows = defaultdict(lambda: {"LungLineage": 0, "BreastLineage": 0, "ProstateLineage": 0})
for e, fk, lin, w, lab in q("""SELECT ?e ?fk ?l ?w ?lab WHERE { ?e a okg:LineageEvidence ; okg:featureKey ?fk ; okg:forLineage ?l ;
                               okg:weight ?w ; rdfs:label ?lab }"""):
    key = str(fk).replace("status.", "")
    lin_rows[key][local(lin)] = float(w)
lineage_matrix = [dict(marker=k, **v) for k, v in sorted(lin_rows.items(), key=lambda kv: -max(kv[1].values()))]

rules = [dict(id=g["id"], label=g["label"], slot=g["slot"], applies=CLASS_INFO.get(g["applies"], {}).get("label", g["applies"]),
              rules=[dict(id=r["id"], prio=r["priority"], assigns=CLASS_INFO[r["assigns"]]["label"], rare=CLASS_INFO[r["assigns"]]["rare"],
                          gl=r["guideline"], crit=[c["label"] for c in r["criteria"]] or ["(default when nothing more specific applies)"])
                     for r in g["rules"]]) for g in GROUPS]
mods = [dict(id=m["id"], label=m["label"], applies=CLASS_INFO.get(m["applies"], {}).get("label", m["applies"]), tx=m["therapy"]) for m in MODS]
params = [dict(name=str(n), value=float(v), gl=label(g), note=str(nt)) for n, v, g, nt in
          q("SELECT ?n ?v ?g ?nt WHERE { ?p a okg:Parameter ; rdfs:label ?n ; okg:paramValue ?v ; okg:hasGuideline ?g ; okg:clinicalNote ?nt } ORDER BY ?n")]
alts = [dict(label=str(l), gene=str(g)) for l, g in q("SELECT ?l ?g WHERE { ?a okg:matchGene ?g ; rdfs:label ?l } ORDER BY ?g")]
n_triples = len(O)

# ---------------- normalization examples
def ex(mask, domain, n=1):
    rows = fnd[mask].drop_duplicates("raw_value").head(n)
    return [dict(domain=domain, raw=r.raw_value, norm=f"{r.feature} = {r.value} {r.unit}".strip(), method=r.method) for r in rows.itertuples()]
E = []
E += ex(fnd.method.str.contains("Allred") & (fnd.feature == "ER_pct"), "IHC · ER (Allred score)")
E += ex(fnd.method.str.contains("ASCO/CAP 2020: <1%") & fnd.raw_value.str.contains("Allred PS2"), "IHC · ER-low from Allred")
E += ex(fnd.method.str.contains("fraction index"), "IHC · Ki-67 as fraction")
E += ex(fnd.raw_value.str.contains("faint", case=False) & (fnd.feature == "HER2_ihc0_subtype"), "IHC · HER2 ultralow")
E += ex(fnd.raw_value.str.contains("ISH group") & (fnd.feature == "HER2_status"), "IHC + ISH · HER2 integrated", 2)
E += ex(fnd.feature.eq("HER2_ish_group") & fnd.raw_value.str.contains("ERBB2:CEP17"), "ISH · ratio + copies")
E += ex(fnd.raw_value.str.contains("22C3 TPS: <1%"), "IHC · PD-L1 TPS")
E += ex(fnd.feature.str.startswith("variant.") & fnd.raw_value.str.contains("Leu858Arg|Glu746|Gly12Cys"), "NGS · 3-letter HGVS, % VAF", 2)
E += ex(fnd.feature.str.startswith("variant.") & fnd.raw_value.str.contains("METex14"), "NGS · MET exon 14 text")
E += ex(fnd.feature.str.startswith("variant.") & fnd.value.str.contains("LOW_VAF") & fnd.raw_value.str.contains("DNMT3A|TET2"), "NGS · CHIP-like variant filtered")
E += ex(fnd.feature.str.startswith("variant.") & fnd.value.str.contains("VUS"), "NGS · VUS not counted")
E += ex(fnd.feature.eq("fusion") & fnd.raw_value.str.contains(r"/"), "RNA fusion · separator")
E += ex(fnd.feature.eq("fusion") & fnd.value.str.contains("LOW_SUPPORT"), "RNA fusion · low read support")
E += ex(fnd.feature.str.startswith("cn.") & fnd.raw_value.str.startswith("FC") & ~fnd.method.str.contains("neutral"), "CNV · fold change → purity-adjusted CN")
E += ex(fnd.feature.str.startswith("cn.") & fnd.raw_value.str.startswith("log2R") & fnd.method.str.contains("homdel"), "CNV · log2R → homozygous deletion")
E += ex(fnd.feature.eq("serum.ProGRP") & fnd.raw_value.str.contains("ng/mL"), "Serum · ProGRP ng/mL → pg/mL")
E += ex(fnd.feature.eq("serum.PSA") & fnd.raw_value.str.contains("µg/L"), "Serum · µg/L → ng/mL")
E += ex(fnd.feature.eq("serum.PSA") & fnd.raw_value.str.startswith("<"), "Serum · below detection limit")
E += ex(fnd.feature.eq("clin.pack_years") & fnd.raw_value.str.contains("former"), "Clinical · smoking text")
E += ex(fnd.feature.eq("purity") & fnd.raw_value.str.contains(" %"), "Clinical · tumour cellularity")
E += ex(fnd.feature.eq("ER_status") & fnd.method.str.contains("RNA fallback"), "RNA · ER fallback when IHC missing")
E += ex(fnd.feature.eq("status.NKX2-1") & fnd.raw_value.str.contains("IHC=None"), "RNA · TTF-1 (NKX2-1) status from RNA")

# ---------------- patients
raw_tables = {}
for name, cols in [("ihc", ["stain", "result"]), ("ngs_variants", ["gene", "hgvs_p", "hgvs_c", "vaf", "depth", "classification"]),
                   ("ngs_fusions", ["fusion", "supporting_reads"]), ("serum", ["analyte", "value", "unit"]),
                   ("ish", ["test", "result_text"]), ("ngs_cnv", ["gene", "metric", "value"])]:
    t = pd.read_csv(RAW / f"{name}.csv", dtype=str).fillna("")
    raw_tables[name] = {pid: d[cols].values.tolist() for pid, d in t.groupby("patient_id")}

methods = []
midx = {}
def mi(m):
    if m not in midx:
        midx[m] = len(methods); methods.append(m)
    return midx[m]
SKIP = ("rna.panel", "site", "specimen_site", "history.prostate_cancer", "history.breast_cancer")
prov = defaultdict(list)
for r in fnd.itertuples():
    if r.feature in SKIP or r.feature.startswith("alt.") or r.feature.startswith("serum.") and r.feature.endswith("_xULN"):
        continue
    prov[r.patient_id].append([r.feature, r.value, r.unit, r.raw_value, mi(r.method)])

pts = []
for r in dx.itertuples():
    pid = r.patient_id
    d, t, c = det.loc[pid], traces[pid], clin.loc[pid]
    grp = [dict(g=x["group"], fired=x["fired"], rules=[[y["rule"], CLASS_INFO[y["assigns"]]["label"], y["status"],
                                                          [[z["label"], z["observed"], z["result"]] for z in y["criteria"]]]
                                                         for y in x["rules"]]) for x in t.get("groups", [])]
    pts.append(dict(
        id=pid, site=c.site, sex=c.sex, age=c.age, spec=c.specimen_site, smoke=c.smoking_history, hist=c.clinical_history,
        th=d.true_histology, tm=d.true_molecular, tmods=d.true_modifiers, edge=d.edge_case, note=d.note,
        h=r.histology, m=r.molecular, lab=r.diagnosis_label, conf=r.confidence, rare=r.rare_subtype == "True",
        mods=r.modifiers, flags=[x for x in r.flags.split(" | ") if x], tx=[x for x in r.therapy_relevance.split(" | ") if x],
        lin=r.lineage, ho=d.hist_outcome, mo=d.mol_outcome, lo=d.lineage_outcome,
        ls=t["lineage"]["scores"], lc=t["lineage"]["contributions"], lm=t["lineage"]["margin"],
        ms=t["malignancy"]["score"], mc=t["malignancy"]["contributions"], mv=t["malignancy"]["verdict"],
        grp=grp, raw={k: raw_tables[k].get(pid, []) for k in raw_tables}, prov=prov[pid]))

gtc = pd.read_csv(ROOT / "data" / "ground_truth.csv", dtype=str)
cohort = gtc.scenario.str.split(":").str[0].value_counts().to_dict()
n_raw = sum(sum(1 for _ in open(f)) - 1 for f in RAW.glob("*.csv") if "rnaseq" not in f.name)
n_rna = sum(pd.read_csv(f, index_col=0).shape[1] for f in RAW.glob("rnaseq_counts_*.csv"))
bundle = dict(
    meta=dict(cohort=cohort, n_raw=n_raw, n_rna=n_rna, n=len(pts), triples=n_triples, n_classes=len(classes), n_rules=sum(len(g["rules"]) for g in rules),
              n_mods=len(mods), n_alts=len(alts), n_findings=len(fnd)),
    ev={k: v for k, v in ev.items() if k not in ("errors", "edge_cases")}, rep=rep,
    classes=classes, lm=lineage_matrix, rules=rules, mods=mods, params=params, alts=alts, ex=E,
    labels={k: v["label"] for k, v in CLASS_INFO.items()}, methods=methods, pts=pts,
    modlabels={m["id"]: m["label"] for m in MODS},
)
data = json.dumps(bundle, separators=(",", ":"), default=str)
html = (ROOT / "pipeline" / "explorer_template.html").read_text().replace("__DATA__", data.replace("</", "<\\/"))
(OUT / "oncokg_explorer.html").write_text(html)
print(f"explorer: {len(html) / 1e6:.2f} MB, patients={len(pts)}, examples={len(E)}")
