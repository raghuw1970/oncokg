"""
Ontology-driven patient identification.

Nothing diagnostic is hard-coded here: lineage evidence, specimen-site priors, malignancy
evidence, rule groups, rules, criteria, modifiers and thresholds are all read from
ontology/oncokg.ttl with SPARQL, and evaluated against output/features.json.

Algorithm per patient
  1. Tissue of origin  : sum LineageEvidence weights whose criteria hold + specimen-site prior;
                         lineages with a requiredSex the patient lacks are excluded.
                         Assign if top >= LINEAGE_MIN_SCORE and margin >= LINEAGE_MIN_MARGIN.
  2. Malignancy        : sum MalignancyEvidence weights -> malignant / indeterminate / no malignancy.
  3. Subtype           : walk RuleGroups in groupOrder; a group applies when its target class is
                         (an ancestor of) a class already assigned; first rule (by priority)
                         whose criteria all hold assigns its class. Unknown (missing) data never
                         satisfies a criterion unless the criterion says missingIs "pass".
                         Numeric criteria on range-valued features (e.g. Allred-derived %) that
                         straddle a cut-off evaluate to 'uncertain'.
  4. Modifiers         : every MolecularModifier whose target class is assigned and criteria hold.
  5. Output            : diagnosis, rarity, therapy relevance, flags, full evidence trace, and an
                         RDF knowledge graph of patients linked to ontology classes.
"""
import json
from collections import defaultdict
import pandas as pd
from rdflib import Graph, Namespace, Literal, URIRef, RDF, RDFS, OWL, XSD
from common import ROOT, OKG, q, params, local, ontology, label

OUT = ROOT / "output"
PRM = params()
O = ontology()
feats = json.load(open(OUT / "features.json"))

# ------------------------------------------------------------------ class hierarchy
PARENT = defaultdict(set)
for c, p in q("SELECT ?c ?p WHERE { ?c rdfs:subClassOf ?p }"):
    PARENT[local(c)].add(local(p))


def ancestors(c):
    out, stack = {c}, [c]
    while stack:
        for p in PARENT[stack.pop()]:
            if p not in out:
                out.add(p); stack.append(p)
    return out


CLASS_INFO = {}
for c in set(PARENT) | {p for ps in PARENT.values() for p in ps}:
    u = OKG[c]
    rare = O.value(u, OKG.isRareSubtype)
    CLASS_INFO[c] = dict(label=label(u), rare=bool(rare.toPython()) if rare is not None else False,
                         tier=str(O.value(u, OKG.rarityTier) or ""), prevalence=str(O.value(u, OKG.approxPrevalence) or ""),
                         therapy=str(O.value(u, OKG.therapyRelevance) or ""))


# ------------------------------------------------------------------ criteria
def load_criterion(c):
    op = str(O.value(c, OKG.operator))
    d = dict(id=local(c), op=op, label=str(O.value(c, RDFS.label)))
    if op == "anyOf":
        d["any"] = [load_criterion(x) for x in O.objects(c, OKG.anyOf)]
    else:
        d["feature"] = str(O.value(c, OKG.featureKey))
        v = O.value(c, OKG.criterionValue)
        d["value"] = str(v) if v is not None else None
        d["missing"] = str(O.value(c, OKG.missingIs) or "fail")
    return d


def evaluate(c, f):
    """returns (True|False|None, observed)"""
    if c["op"] == "anyOf":
        res = [evaluate(x, f) for x in c["any"]]
        vals = [r[0] for r in res]
        out = True if True in vals else (None if None in vals else False)
        return out, " | ".join(f"{x['label']}={r[1]}" for x, r in zip(c["any"], res))
    key, op, val = c["feature"], c["op"], c["value"]
    if key not in f or f[key] is None:
        return (True if c["missing"] == "pass" else None), "missing"
    x = f[key]
    if op in ("ge", "gt", "lt", "le"):
        t = float(val)
        cmp = {"ge": lambda a: a >= t, "gt": lambda a: a > t, "lt": lambda a: a < t, "le": lambda a: a <= t}[op]
        lo, hi = f.get(key + "__lo"), f.get(key + "__hi")
        if lo is not None and hi is not None and lo != hi:
            a, b = cmp(float(lo)), cmp(float(hi))
            return (a if a == b else None), f"{x} (range {lo}-{hi})"
        return cmp(float(x)), x
    sx = str(x).lower()
    if op == "eq":
        return sx == val.lower(), x
    opts = [o.lower() for o in val.split("|")]
    return ((sx in opts) if op == "in" else (sx not in opts)), x


def criteria_of(node):
    return [load_criterion(c) for c in O.objects(node, OKG.hasCriterion)]


# ------------------------------------------------------------------ knowledge from the ontology
LINEAGES = {}
for l, root, lab in q("SELECT ?l ?r ?lab WHERE { ?l a okg:Lineage ; okg:lineageRootDisease ?r ; rdfs:label ?lab }"):
    sex = O.value(l, OKG.requiredSex)
    LINEAGES[local(l)] = dict(root=local(root), label=str(lab), sex=str(sex) if sex else None)

LIN_EVID = []
for e, fk, op, v, lin, w, lab in q("""SELECT ?e ?fk ?op ?v ?l ?w ?lab WHERE { ?e a okg:LineageEvidence ; okg:featureKey ?fk ;
        okg:operator ?op ; okg:criterionValue ?v ; okg:forLineage ?l ; okg:weight ?w ; rdfs:label ?lab }"""):
    LIN_EVID.append(dict(c=dict(id=local(e), op=str(op), feature=str(fk), value=str(v), missing="fail", label=str(lab)),
                         lineage=local(lin), w=float(w)))
SITE_PRIOR = [(str(k), local(l), float(w)) for k, l, w in
              q("SELECT ?k ?l ?w WHERE { ?s a okg:SpecimenSitePrior ; okg:siteKeyword ?k ; okg:forLineage ?l ; okg:weight ?w }")]
MAL_EVID = []
for e, fk, op, v, w, lab in q("""SELECT ?e ?fk ?op ?v ?w ?lab WHERE { ?e a okg:MalignancyEvidence ; okg:featureKey ?fk ;
        okg:operator ?op ; okg:criterionValue ?v ; okg:weight ?w ; rdfs:label ?lab }"""):
    MAL_EVID.append(dict(c=dict(id=local(e), op=str(op), feature=str(fk), value=str(v), missing="fail", label=str(lab)), w=float(w)))

RULE_BY_ID = {}
BENIGN_REQUIRES = [str(x) for x in O.objects(OKG.NoMalignancyDetected, OKG.requiresFeature)]
GROUPS = []
for grp, applies, order, slot, lab in q("""SELECT ?g ?a ?o ?s ?lab WHERE { ?g a okg:RuleGroup ; okg:groupAppliesTo ?a ;
        okg:groupOrder ?o ; okg:outputSlot ?s ; rdfs:label ?lab } ORDER BY ?o"""):
    rules = []
    for r, prio, assigns in q(f"""SELECT ?r ?p ?c WHERE {{ ?r okg:inRuleGroup <{grp}> ; okg:priority ?p ; okg:ruleAssigns ?c }} ORDER BY ?p"""):
        gl = O.value(r, OKG.hasGuideline)
        rules.append(dict(id=local(r), priority=int(prio), assigns=local(assigns), criteria=criteria_of(r),
                          guideline=label(gl) if gl is not None else ""))
    RULE_BY_ID.update({r["id"]: r for r in rules})
    GROUPS.append(dict(id=local(grp), applies=local(applies), order=int(order), slot=str(slot), label=str(lab), rules=rules))

MODS = []
for m, applies, lab in q("SELECT ?m ?a ?lab WHERE { ?m a okg:MolecularModifier ; okg:modifierAppliesTo ?a ; rdfs:label ?lab }"):
    tx = O.value(m, OKG.therapyRelevance)
    MODS.append(dict(id=local(m), applies=local(applies), label=str(lab), criteria=criteria_of(m), therapy=str(tx) if tx else ""))
MODS.sort(key=lambda d: d["id"])


# ------------------------------------------------------------------ reasoning
def reason(pid, f):
    flags = list(f.get("flags", []))
    trace = dict(patient_id=pid)
    # 1. lineage
    scores, contrib = {}, {}
    for lin, info in LINEAGES.items():
        if info["sex"] and f.get("clin.sex") != info["sex"]:
            scores[lin], contrib[lin] = None, [("excluded: requires sex " + info["sex"], 0)]
            continue
        s, cl = 0.0, []
        for e in LIN_EVID:
            if e["lineage"] == lin and evaluate(e["c"], f)[0] is True:
                s += e["w"]; cl.append((e["c"]["label"], e["w"]))
        site = f.get("specimen_site", "")
        for kw, l2, w in SITE_PRIOR:
            if l2 == lin and kw in site:
                s += w; cl.append((f"specimen site '{site}'", w)); break
        scores[lin], contrib[lin] = round(s, 2), cl
    ranked = sorted([(v, k) for k, v in scores.items() if v is not None], reverse=True)
    top_s, top_l = ranked[0]
    second = ranked[1][0] if len(ranked) > 1 else 0.0
    margin = round(top_s - second, 2)
    lineage = top_l if (top_s >= PRM["LINEAGE_MIN_SCORE"] and margin >= PRM["LINEAGE_MIN_MARGIN"]) else None
    trace["lineage"] = dict(scores=scores, contributions={k: v for k, v in contrib.items()}, assigned=lineage, margin=margin)
    if lineage is None:
        flags.append(f"Tissue of origin not resolved (top {LINEAGES[top_l]['label']} {top_s}, margin {margin}) - "
                     "consider broader lineage IHC panel / clinical correlation")

    # 2. malignancy
    ms, mc = 0.0, []
    for e in MAL_EVID:
        if evaluate(e["c"], f)[0] is True:
            ms += e["w"]; mc.append((e["c"]["label"], e["w"]))
    verdict = "malignant" if ms >= PRM["MALIGNANCY_MIN_SCORE"] else "no_malignancy" if ms <= PRM["BENIGN_MAX_SCORE"] else "equivocal"
    missing_req = [k for k in BENIGN_REQUIRES if k not in f]
    if verdict == "no_malignancy" and missing_req:
        verdict = "equivocal"
        flags.append("No malignancy evidence, but required data absent (" + ", ".join(missing_req) +
                     ") - cannot exclude malignancy")
    trace["malignancy"] = dict(score=ms, contributions=mc, verdict=verdict)

    assigned, fired, groups_trace = [], [], []
    slots = {"histology": None, "molecular": None}
    if verdict == "no_malignancy" and lineage:
        assigned = ["NoMalignancyDetected"]; slots["histology"] = "NoMalignancyDetected"
    elif verdict == "equivocal" or (verdict == "no_malignancy" and not lineage):
        assigned = ["IndeterminateDiagnosis"]; slots["histology"] = "IndeterminateDiagnosis"
        flags.append(f"Malignancy evidence equivocal (score {ms}) - histopathologic correlation required")
    elif lineage is None:
        assigned = ["Cancer"]; slots["histology"] = "IndeterminateDiagnosis"
        if f.get("NE_count", 0) >= 2:
            flags.append("High-grade neuroendocrine program present; primary site undetermined")
    else:
        assigned = [LINEAGES[lineage]["root"]]
        for grp in GROUPS:
            anc = set().union(*(ancestors(c) for c in assigned))
            if grp["applies"] not in anc:
                continue
            gt = dict(group=grp["label"], slot=grp["slot"], rules=[], fired=None)
            winners = []
            for r in grp["rules"]:
                res = [(c, *evaluate(c, f)) for c in r["criteria"]]
                ok = all(x[1] is True for x in res)
                status = "fired" if ok else ("blocked: missing/uncertain data" if any(x[1] is None for x in res) and
                                              not any(x[1] is False for x in res) else "not met")
                gt["rules"].append(dict(rule=r["id"], assigns=r["assigns"], status=status, guideline=r["guideline"],
                                        criteria=[dict(label=c["label"], observed=str(obs), result=res_) for c, res_, obs in res]))
                if ok:
                    winners.append(r)
            if winners:
                w = winners[0]
                assigned.append(w["assigns"]); fired.append(w["id"]); gt["fired"] = w["id"]
                slots[grp["slot"]] = w["assigns"]
                if grp["id"] == "RG_nsclc_driver" and len(winners) > 1:
                    flags.append("Multiple NSCLC drivers detected: " + ", ".join(x["assigns"] for x in winners) +
                                 " - verify clonality / resistance context")
            else:
                blocked = [r for r in gt["rules"] if r["status"].startswith("blocked")]
                if grp["slot"] == "histology" and blocked:
                    flags.append(f"{grp['label']}: no rule satisfied; blocked by missing/uncertain data for "
                                 + ", ".join(sorted({c['label'] for r in blocked for c in r['criteria'] if c['result'] is None})))
            groups_trace.append(gt)
        if slots["histology"] is None:
            slots["histology"] = LINEAGES[lineage]["root"]
    trace["groups"] = groups_trace

    # 4. modifiers
    anc = set().union(*(ancestors(c) for c in assigned))
    mods = []
    for m in MODS:
        if m["applies"] in anc and m["criteria"] and all(evaluate(c, f)[0] is True for c in m["criteria"]):
            mods.append(m["id"])
    # most specific classes
    specific = [c for c in assigned if not any(c != o and c in ancestors(o) for o in assigned)]
    rare = [c for c in specific if CLASS_INFO.get(c, {}).get("rare")]
    if any(c in ("TripleNegativeBreastCancer", "TNBCLuminalAndrogenReceptor", "ERLowHER2NegativeBreastCancer") for c in assigned) \
            and "MOD_BRCA1" not in mods and "MOD_BRCA2" not in mods:
        flags.append("TNBC/ER-low: recommend germline BRCA1/2 testing (PARP inhibitor eligibility)")
    # confidence
    hist = slots["histology"]
    if hist in ("IndeterminateDiagnosis",) or lineage is None:
        conf = "Low"
    elif hist in ("NonSmallCellLungCancer", "BreastCancer", "HRPositiveHER2NegativeBreastCancer") or len(flags) >= 2:
        conf = "Moderate"
    elif margin >= 4 and (ms >= 4 or hist == "NoMalignancyDetected"):
        conf = "High"
    else:
        conf = "Moderate"
    lab = CLASS_INFO.get(hist, {}).get("label", hist)
    if slots["molecular"]:
        lab += " - " + CLASS_INFO[slots["molecular"]]["label"]
    therapy = [CLASS_INFO[c]["therapy"] for c in specific if CLASS_INFO.get(c, {}).get("therapy")] + \
              [f"{m['label']}: {m['therapy']}" for m in MODS if m["id"] in mods and m["therapy"]]
    return dict(patient_id=pid, lineage=lineage or "", lineage_margin=margin, malignancy_score=ms, malignancy=verdict,
                histology=hist, molecular=slots["molecular"] or "", diagnosis_label=lab, classes=";".join(specific),
                rare_subtype=bool(rare), rare_classes=";".join(rare), modifiers=";".join(mods), confidence=conf,
                rules_fired=";".join(fired), therapy_relevance=" | ".join(therapy), flags=" | ".join(flags)), trace


def main():
    results, traces = [], {}
    for pid, f in feats.items():
        r, t = reason(pid, f)
        results.append(r); traces[pid] = t
    res = pd.DataFrame(results)
    res.to_csv(OUT / "patient_diagnoses.csv", index=False)
    json.dump(traces, open(OUT / "evidence_trace.json", "w"), default=str)

    # ------------------------------------------------------------------ patient knowledge graph
    PT = Namespace("http://example.org/oncokg/patient/")
    K = Graph()
    K.bind("okg", OKG); K.bind("pt", PT)
    kg_onto = URIRef("http://example.org/oncokg/patients")
    K.add((kg_onto, RDF.type, OWL.Ontology)); K.add((kg_onto, OWL.imports, URIRef("http://example.org/oncokg")))
    FK2BM = {str(k): b for b, k in q("SELECT ?b ?k WHERE { ?b a ?t ; okg:featureKey ?k . ?t rdfs:subClassOf* okg:Biomarker }")}
    USED = {c["feature"] for g in GROUPS for r in g["rules"] for c in
            [x for cr in r["criteria"] for x in ([cr] + cr.get("any", []))] if "feature" in c}
    for r in results:
        pid, f = r["patient_id"], feats[r["patient_id"]]
        p, s, t, a = PT[pid], PT[pid + "_specimen"], PT[pid + "_tumor"], PT[pid + "_dx"]
        K.add((p, RDF.type, OKG.Patient)); K.add((p, OKG.hasSpecimen, s)); K.add((s, RDF.type, OKG.Specimen))
        K.add((s, OKG.sourceSite, Literal(f.get("site", "")))); K.add((s, RDFS.label, Literal(f.get("specimen_site", ""))))
        K.add((p, OKG.hasTumor, t)); K.add((t, RDF.type, OKG.Tumor))
        for c in r["classes"].split(";"):
            if c: K.add((t, RDF.type, OKG[c]))
        for m in filter(None, r["modifiers"].split(";")):
            K.add((t, OKG.hasModifier, OKG[m]))
        if r["lineage"]:
            K.add((p, OKG.assignedLineage, OKG[r["lineage"]]))
        K.add((p, OKG.hasDiagnosis, a)); K.add((a, RDF.type, OKG.DiagnosisAssertion))
        K.add((a, OKG.confidence, Literal(r["confidence"]))); K.add((a, RDFS.label, Literal(r["diagnosis_label"])))
        for c in r["classes"].split(";"):
            if c: K.add((a, OKG.assertsClass, OKG[c]))
        for rid in filter(None, r["rules_fired"].split(";")):
            K.add((a, OKG.derivedByRule, OKG[rid]))
        fired_feats = {c["feature"] for rid in filter(None, r["rules_fired"].split(";")) for cr in RULE_BY_ID[rid]["criteria"]
                       for c in [cr] + cr.get("any", []) if "feature" in c}
        for fl in filter(None, r["flags"].split(" | ")):
            K.add((a, OKG.clinicalNote, Literal(fl)))
        for k, v in f.items():
            if k == "flags" or k.startswith("rna.") or k.endswith("__lo") or k.endswith("__hi") or v is None:
                continue
            if not (k.startswith(("status.", "alt.", "serum.")) or k in USED or k in ("HER2_status", "Ki67_pct", "PDL1_TPS", "purity", "tmb", "fga")):
                continue
            if k.startswith("alt.") and v is False:
                continue
            fn = PT[f"{pid}_f_{k.replace('.', '_').replace('-', '_')}"]
            K.add((s, OKG.hasFinding, fn)); K.add((fn, RDF.type, OKG.NormalizedFinding))
            K.add((fn, OKG.featureKey, Literal(k))); K.add((fn, OKG.normalizedValue, Literal(str(v))))
            if k in FK2BM:
                K.add((fn, OKG.aboutBiomarker, FK2BM[k]))
            if k in fired_feats:
                K.add((a, OKG.supportedByFinding, fn))
    K.serialize(OUT / "patient_kg.ttl", format="turtle")

    print(f"patients={len(res)}  KG triples={len(K)}")
    print(res.histology.value_counts().to_string())
    print("\nmolecular:\n" + res.molecular.value_counts().to_string())
    print("\nconfidence:\n" + res.confidence.value_counts().to_string())


if __name__ == "__main__":
    main()
