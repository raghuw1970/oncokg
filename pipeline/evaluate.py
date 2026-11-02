"""Evaluate reasoner output against hidden ground truth (data/ground_truth.csv)."""
import json
from collections import Counter, defaultdict
import pandas as pd
from common import ROOT
from reason import ancestors, CLASS_INFO  # reuse ontology hierarchy

OUT = ROOT / "output"
gt = pd.read_csv(ROOT / "data" / "ground_truth.csv", dtype=str).fillna("")
pr = pd.read_csv(OUT / "patient_diagnoses.csv", dtype=str, keep_default_na=False)
df = gt.merge(pr, on="patient_id")
LIN = {"lung": "LungLineage", "breast": "BreastLineage", "prostate": "ProstateLineage"}


def hist_outcome(t, p):
    if p == t:
        return "correct"
    if p in ("IndeterminateDiagnosis", ""):
        return "abstained"
    if t != "NoMalignancyDetected" and p in ancestors(t) - {"Cancer", "Disease"}:
        return "correct_less_specific"
    return "wrong"


df["lineage_ok"] = df.apply(lambda r: r.lineage == LIN[r.true_lineage], axis=1)
df["lineage_outcome"] = df.apply(lambda r: "correct" if r.lineage_ok else ("abstained" if r.lineage == "" else "wrong"), axis=1)
df["hist_outcome"] = df.apply(lambda r: hist_outcome(r.true_histology, r.histology), axis=1)
df["mol_outcome"] = df.apply(lambda r: ("correct" if r.molecular == r.true_molecular else
                                        "missed" if r.true_molecular and not r.molecular else
                                        "false_positive" if r.molecular and not r.true_molecular else
                                        "wrong") if (r.true_molecular or r.molecular) else "n/a", axis=1)

ev = {}
ev["n_patients"] = len(df)
ev["lineage"] = dict(Counter(df.lineage_outcome))
ev["histology_subtype"] = dict(Counter(df.hist_outcome))
tumors = df[df.true_histology != "NoMalignancyDetected"]
benign = df[df.true_histology == "NoMalignancyDetected"]
ev["malignancy"] = dict(
    cancers_called_malignant=int((tumors.malignancy == "malignant").sum()), cancers=len(tumors),
    cancers_called_benign=int((tumors.histology == "NoMalignancyDetected").sum()),
    benign_called_benign=int((benign.histology == "NoMalignancyDetected").sum()), benign=len(benign),
    benign_called_cancer=int((benign.malignancy == "malignant").sum()))
ev["molecular"] = dict(Counter(df.mol_outcome[df.mol_outcome != "n/a"]))

# per-class sensitivity / PPV (most-specific truth class = molecular if any else histology)
df["true_specific"] = df.true_molecular.where(df.true_molecular != "", df.true_histology)
df["pred_specific"] = df.molecular.where(df.molecular != "", df.histology)
per = []
for c in sorted(set(df.true_specific)):
    tp = int(((df.true_specific == c) & (df.pred_specific == c)).sum())
    n = int((df.true_specific == c).sum())
    npred = int((df.pred_specific == c).sum())
    info = CLASS_INFO.get(c, {})
    per.append(dict(cls=c, label=info.get("label", c), rare=info.get("rare", False), tier=info.get("tier", ""), n_true=n,
                    n_pred=npred, tp=tp, sensitivity=round(tp / n, 3) if n else None,
                    ppv=round(tp / npred, 3) if npred else None))
ev["per_class"] = per
rare_rows = [x for x in per if x["rare"]]
ev["rare_subtypes"] = dict(n_true=sum(x["n_true"] for x in rare_rows), tp=sum(x["tp"] for x in rare_rows),
                           n_pred=sum(x["n_pred"] for x in rare_rows))

# TNBC focus: anything biologically TNBC (incl. LAR) vs predicted TNBC family
tn_true = df.true_histology == "TripleNegativeBreastCancer"
tn_pred = df.histology == "TripleNegativeBreastCancer"
ev["tnbc"] = dict(true=int(tn_true.sum()), predicted=int(tn_pred.sum()), tp=int((tn_true & tn_pred).sum()),
                  sensitivity=round((tn_true & tn_pred).sum() / tn_true.sum(), 3),
                  ppv=round((tn_true & tn_pred).sum() / tn_pred.sum(), 3) if tn_pred.sum() else None,
                  missed_as=dict(Counter(df.histology[tn_true & ~tn_pred])))

# modifiers
modstats = []
allmods = sorted({m for s in df.true_modifiers for m in s.split(";") if m} | {m for s in df.modifiers for m in s.split(";") if m})
for m in allmods:
    t = df.true_modifiers.str.split(";").apply(lambda x: m in x)
    p_ = df.modifiers.str.split(";").apply(lambda x: m in x)
    tp = int((t & p_).sum())
    modstats.append(dict(modifier=m, n_true=int(t.sum()), n_pred=int(p_.sum()), tp=tp,
                         recall=round(tp / t.sum(), 3) if t.sum() else None, precision=round(tp / p_.sum(), 3) if p_.sum() else None))
ev["modifiers"] = modstats

# error review
errs = df[(df.hist_outcome == "wrong") | (df.mol_outcome.isin(["wrong", "false_positive"])) | (df.lineage_outcome == "wrong")]
ev["errors"] = errs[["patient_id", "scenario", "true_histology", "true_molecular", "histology", "molecular", "note", "flags"]].to_dict("records")
abst = df[df.hist_outcome == "abstained"]
ev["abstentions_by_reason"] = dict(Counter(
    ("lineage unresolved" if r.lineage == "" else "malignancy equivocal" if r.malignancy == "equivocal" else "other")
    for r in abst.itertuples()))
ev["edge_cases"] = df[df.edge_case != ""][["patient_id", "edge_case", "histology", "flags"]].to_dict("records")
df.to_csv(OUT / "evaluation_detail.csv", index=False)
json.dump(ev, open(OUT / "evaluation.json", "w"), indent=1, default=str)

print(json.dumps({k: ev[k] for k in ("lineage", "histology_subtype", "malignancy", "molecular", "rare_subtypes", "tnbc", "abstentions_by_reason")}, indent=1))
print("\nper-class:")
print(pd.DataFrame(per)[["cls", "rare", "n_true", "n_pred", "tp", "sensitivity", "ppv"]].to_string(index=False))
print(f"\nerrors: {len(errs)}")
print(errs[["patient_id", "scenario", "true_histology", "histology", "molecular", "note"]].to_string(index=False))
