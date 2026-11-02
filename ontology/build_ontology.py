"""
Build OncoKG: an OWL ontology + knowledge graph for identifying lung, breast and
prostate cancer (incl. rare molecular subtypes such as triple-negative breast
cancer) from molecular biomarker data.

The ontology is the single source of truth for the pipeline:
  * biomarker vocabulary + synonyms (used by the normalizer to map raw names)
  * units, conversion factors and upper limits of normal (serum markers)
  * RNA panel definition (housekeeping genes, effective lengths, call thresholds)
  * guideline thresholds (ASCO/CAP ER/PR/HER2, ALK FISH, PD-L1, QC cut-offs)
  * genomic alteration definitions (HGVS / fusion / copy-number matchers)
  * lineage evidence weights (tissue-of-origin scoring)
  * malignancy evidence weights
  * diagnostic rules (reified criteria) and molecular modifiers
  * disease class hierarchy with prevalence + rarity annotations

Output: oncokg.ttl (Turtle)
"""
from pathlib import Path
from rdflib import Graph, Namespace, Literal, URIRef, BNode, RDF, RDFS, OWL, XSD
from rdflib.namespace import SKOS, DCTERMS

OKG = Namespace("http://example.org/oncokg#")
g = Graph()
g.bind("okg", OKG); g.bind("skos", SKOS); g.bind("owl", OWL); g.bind("dcterms", DCTERMS)

ONTO = URIRef("http://example.org/oncokg")
g.add((ONTO, RDF.type, OWL.Ontology))
g.add((ONTO, DCTERMS.title, Literal("OncoKG - Molecular Biomarker Ontology for Lung, Breast and Prostate Cancer Identification")))
g.add((ONTO, OWL.versionInfo, Literal("1.0.0")))
g.add((ONTO, DCTERMS.created, Literal("2026-11-02", datatype=XSD.date)))
g.add((ONTO, RDFS.comment, Literal(
    "Research/educational knowledge model. Thresholds follow published guidelines where cited, "
    "but this is NOT a validated clinical decision tool.")))


def lit(v):
    if isinstance(v, bool):
        return Literal(v, datatype=XSD.boolean)
    if isinstance(v, int):
        return Literal(v, datatype=XSD.integer)
    if isinstance(v, float):
        return Literal(v, datatype=XSD.decimal)
    return Literal(v)


def cls(name, parent=None, label=None, comment=None, **ann):
    u = OKG[name]
    g.add((u, RDF.type, OWL.Class))
    if parent:
        for p in (parent if isinstance(parent, (list, tuple)) else [parent]):
            g.add((u, RDFS.subClassOf, OKG[p]))
    g.add((u, RDFS.label, Literal(label or name)))
    if comment:
        g.add((u, RDFS.comment, Literal(comment)))
    for k, v in ann.items():
        if v is None:
            continue
        g.add((u, OKG[k], lit(v)))
    return u


def oprop(name, domain=None, rng=None, label=None):
    u = OKG[name]
    g.add((u, RDF.type, OWL.ObjectProperty))
    g.add((u, RDFS.label, Literal(label or name)))
    if domain: g.add((u, RDFS.domain, OKG[domain]))
    if rng: g.add((u, RDFS.range, OKG[rng]))


def dprop(name, domain=None, rng=XSD.string, label=None):
    u = OKG[name]
    g.add((u, RDF.type, OWL.DatatypeProperty))
    g.add((u, RDFS.label, Literal(label or name)))
    if domain: g.add((u, RDFS.domain, OKG[domain]))
    if rng: g.add((u, RDFS.range, rng))


def aprop(name):
    g.add((OKG[name], RDF.type, OWL.AnnotationProperty))


def ind(name, typ, label=None, **props):
    u = OKG[name]
    g.add((u, RDF.type, OWL.NamedIndividual))
    for t in (typ if isinstance(typ, (list, tuple)) else [typ]):
        g.add((u, RDF.type, OKG[t]))
    if label:
        g.add((u, RDFS.label, Literal(label)))
        g.add((u, SKOS.prefLabel, Literal(label)))
    for k, v in props.items():
        if v is None:
            continue
        vals = v if isinstance(v, (list, tuple)) else [v]
        for x in vals:
            if isinstance(x, URIRef):
                g.add((u, OKG[k], x))
            else:
                g.add((u, OKG[k], lit(x)))
    return u


# ---------------------------------------------------------------------------
# 1. Upper-level classes
# ---------------------------------------------------------------------------
for n, p, c in [
    ("Patient", None, "A person whose specimens are profiled."),
    ("Specimen", None, "Tissue, cytology or blood sample."),
    ("Tumor", None, "The neoplasm (or tissue) characterised in a patient; typed with Disease classes."),
    ("Assay", None, "Laboratory method producing observations."),
    ("RawObservation", None, "A value as exported by a source laboratory, before normalization."),
    ("NormalizedFinding", None, "A canonical, unit-harmonised finding derived from raw observations."),
    ("Biomarker", None, "Measurable molecular characteristic."),
    ("Gene", None, "HGNC gene."),
    ("Lineage", None, "Tissue-of-origin lineage program."),
    ("Disease", None, "Diagnostic category."),
    ("MolecularModifier", None, "Non-exclusive molecular annotation of a tumor (e.g. HER2-low, BRCA1-altered)."),
    ("RuleGroup", None, "Set of mutually exclusive diagnostic rules evaluated in priority order."),
    ("DiagnosticRule", None, "Rule assigning a Disease class when all criteria are satisfied."),
    ("Criterion", None, "Atomic test on a normalized feature."),
    ("LineageEvidence", None, "Weighted link from a feature state to a lineage (tissue-of-origin scoring)."),
    ("MalignancyEvidence", None, "Weighted feature state supporting malignancy."),
    ("SpecimenSitePrior", None, "Prior weight from anatomic specimen site to a lineage."),
    ("Unit", None, "Unit of measure."),
    ("UnitConversion", None, "Factor to convert a source unit to the canonical unit."),
    ("Parameter", None, "Named numeric threshold used by normalization or reasoning."),
    ("Guideline", None, "Published source of a rule or threshold."),
    ("DiagnosisAssertion", None, "Provenance record for an inferred diagnosis."),
]:
    cls(n, p, comment=c)

for n, p in [("IHCAssay", "Assay"), ("ISHAssay", "Assay"), ("DNASequencingAssay", "Assay"),
             ("RNAFusionAssay", "Assay"), ("RNAExpressionAssay", "Assay"),
             ("CopyNumberAssay", "Assay"), ("SerumImmunoassay", "Assay")]:
    cls(n, p)

for n, p, c in [
    ("ProteinExpressionBiomarker", "Biomarker", "Assessed by IHC and/or RNA expression of the encoding gene."),
    ("SerumBiomarker", "Biomarker", "Circulating protein measured by immunoassay."),
    ("CompositeBiomarker", "Biomarker", "Status integrated from several assays by a guideline algorithm (e.g. HER2 status)."),
    ("GenomicAlteration", "Biomarker", "Defined DNA/RNA alteration."),
    ("SequenceVariantAlteration", "GenomicAlteration", None),
    ("GeneFusionAlteration", "GenomicAlteration", None),
    ("CopyNumberAlteration", "GenomicAlteration", None),
    ("HousekeepingGene", "Gene", "Reference gene used for RNA normalization."),
]:
    cls(n, p, comment=c)

# Properties
for n, d, r in [
    ("hasSpecimen", "Patient", "Specimen"), ("hasTumor", "Patient", "Tumor"),
    ("hasFinding", "Specimen", "NormalizedFinding"), ("aboutBiomarker", "NormalizedFinding", "Biomarker"),
    ("measuredBy", "NormalizedFinding", "Assay"), ("derivedFrom", "NormalizedFinding", "RawObservation"),
    ("encodedBy", "Biomarker", "Gene"), ("hasLineageEvidence", "Biomarker", "LineageEvidence"),
    ("forLineage", None, "Lineage"), ("lineageRootDisease", "Lineage", "Disease"),
    ("ruleAssigns", "DiagnosticRule", "Disease"), ("inRuleGroup", "DiagnosticRule", "RuleGroup"),
    ("groupAppliesTo", "RuleGroup", "Disease"), ("hasCriterion", None, "Criterion"),
    ("anyOf", "Criterion", "Criterion"), ("modifierAppliesTo", "MolecularModifier", "Disease"),
    ("hasGuideline", None, "Guideline"), ("canonicalUnit", None, "Unit"),
    ("fromUnit", "UnitConversion", "Unit"), ("toUnit", "UnitConversion", "Unit"),
    ("hasDiagnosis", "Patient", "DiagnosisAssertion"), ("assertsClass", "DiagnosisAssertion", None),
    ("derivedByRule", "DiagnosisAssertion", "DiagnosticRule"), ("supportedByFinding", "DiagnosisAssertion", "NormalizedFinding"),
    ("hasModifier", "Tumor", "MolecularModifier"), ("assignedLineage", "Patient", "Lineage"),
]:
    oprop(n, d, r)

for n, d, r in [
    ("featureKey", None, XSD.string), ("requiresFeature", "Disease", XSD.string), ("hgncSymbol", None, XSD.string),
    ("effectiveLengthBp", "Gene", XSD.integer), ("rnaHighThreshold", "Gene", XSD.decimal),
    ("ihcName", None, XSD.string), ("operator", "Criterion", XSD.string),
    ("criterionValue", "Criterion", XSD.string), ("missingIs", "Criterion", XSD.string),
    ("weight", None, XSD.decimal), ("priority", None, XSD.integer), ("groupOrder", "RuleGroup", XSD.integer),
    ("outputSlot", "RuleGroup", XSD.string), ("approxPrevalence", "Disease", XSD.string),
    ("isRareSubtype", "Disease", XSD.boolean), ("rarityTier", "Disease", XSD.string),
    ("requiredSex", "Lineage", XSD.string), ("conversionFactor", "UnitConversion", XSD.decimal),
    ("upperLimitNormal", "SerumBiomarker", XSD.decimal), ("analyteKey", "UnitConversion", XSD.string),
    ("paramValue", "Parameter", XSD.decimal), ("siteKeyword", "SpecimenSitePrior", XSD.string),
    ("matchGene", "GenomicAlteration", XSD.string), ("matchProteinRegex", "SequenceVariantAlteration", XSD.string),
    ("matchCdnaRegex", "SequenceVariantAlteration", XSD.string), ("matchVariantType", "SequenceVariantAlteration", XSD.string),
    ("requiresPathogenic", "GenomicAlteration", XSD.boolean), ("fusion3Prime", "GeneFusionAlteration", XSD.string),
    ("fusion5Prime", "GeneFusionAlteration", XSD.string), ("cnaType", "CopyNumberAlteration", XSD.string),
    ("isDriver", "GenomicAlteration", XSD.boolean), ("therapyRelevance", None, XSD.string),
    ("confidence", "DiagnosisAssertion", XSD.string), ("normalizedValue", "NormalizedFinding", XSD.string),
    ("rawValue", "RawObservation", XSD.string), ("sourceSite", None, XSD.string),
    ("ucumCode", "Unit", XSD.string), ("clinicalNote", None, XSD.string),
]:
    dprop(n, d, r)

# ---------------------------------------------------------------------------
# 2. Guidelines / sources
# ---------------------------------------------------------------------------
GL = {
    "ASCO_CAP_HER2_2018": ("ASCO/CAP HER2 testing in breast cancer 2018 (affirmed 2023)", "https://doi.org/10.1200/JCO.2018.77.8738"),
    "ASCO_CAP_ER_PR_2020": ("ASCO/CAP ER and PgR testing in breast cancer 2020", "https://doi.org/10.1200/JCO.19.02309"),
    "NCCN_BREAST_2026": ("NCCN Breast Cancer 2026: HER2-null / ultralow / low categories", "https://www.nccn.org"),
    "ST_GALLEN_2013": ("St Gallen 2013 surrogate intrinsic subtypes (Ki-67 / PR thresholds)", "https://doi.org/10.1093/annonc/mdt303"),
    "CAP_IASLC_AMP_2018": ("CAP/IASLC/AMP molecular testing guideline for lung cancer 2018", "https://doi.org/10.5858/arpa.2017-0388-CP"),
    "WHO_THORACIC_2021": ("WHO Classification of Thoracic Tumours, 5th ed. (TTF-1/p40 IHC)", "https://publications.iarc.fr"),
    "WHO_URINARY_MALE_GENITAL_2022": ("WHO Classification of Urinary and Male Genital Tumours, 5th ed.", "https://publications.iarc.fr"),
    "LEHMANN_TNBC": ("Lehmann et al. TNBC molecular subtypes incl. LAR", "https://doi.org/10.1172/JCI45014"),
    "AMP_ASCO_CAP_2017": ("AMP/ASCO/CAP somatic variant interpretation 2017", "https://doi.org/10.1016/j.jmoldx.2016.10.002"),
    "SYNTHETIC_DESIGN": ("OncoKG synthetic design choice (not a clinical guideline)", None),
}
for k, (t, url) in GL.items():
    ind(k, "Guideline", t, clinicalNote=url)

# ---------------------------------------------------------------------------
# 3. Parameters (thresholds) used by normalization and reasoning
# ---------------------------------------------------------------------------
PARAMS = [
    # name, value, guideline, note
    ("ER_LOW_MIN_PCT", 1, "ASCO_CAP_ER_PR_2020", "ER <1% = negative"),
    ("ER_LOW_MAX_PCT", 10, "ASCO_CAP_ER_PR_2020", "ER 1-10% = ER-low positive"),
    ("PR_POS_MIN_PCT", 1, "ASCO_CAP_ER_PR_2020", "PR >=1% positive"),
    ("PR_LUMA_MIN_PCT", 20, "ST_GALLEN_2013", "PR >=20% for luminal A-like"),
    ("KI67_LUMB_MIN_PCT", 20, "ST_GALLEN_2013", "Ki-67 >=20% luminal B-like"),
    ("HER2_ISH_RATIO_CUT", 2.0, "ASCO_CAP_HER2_2018", "HER2/CEP17 ratio"),
    ("HER2_ISH_COPY_LOW", 4.0, "ASCO_CAP_HER2_2018", "average HER2 copies/cell lower cut"),
    ("HER2_ISH_COPY_HIGH", 6.0, "ASCO_CAP_HER2_2018", "average HER2 copies/cell upper cut"),
    ("HER2_ULTRALOW_MAX_PCT", 10, "NCCN_BREAST_2026", "IHC 0 with faint/incomplete staining in <=10% cells"),
    ("AR_POS_MIN_PCT", 10, "LEHMANN_TNBC", "AR IHC >=10% used for LAR-like TNBC"),
    ("ALK_FISH_POS_PCT", 15, "CAP_IASLC_AMP_2018", ">=15% cells with split signals"),
    ("PDL1_HIGH_TPS", 50, "CAP_IASLC_AMP_2018", "TPS >=50% high"),
    ("PDL1_POS_TPS", 1, "CAP_IASLC_AMP_2018", "TPS >=1% positive"),
    ("MIN_VARIANT_DEPTH", 100, "AMP_ASCO_CAP_2017", "minimum read depth"),
    ("MIN_VAF", 0.05, "SYNTHETIC_DESIGN", "VAF reporting floor (filters CHIP-like noise)"),
    ("MIN_FUSION_READS", 5, "SYNTHETIC_DESIGN", "minimum supporting reads"),
    ("CN_AMP_MIN", 6.0, "SYNTHETIC_DESIGN", "purity-adjusted copy number for amplification"),
    ("CN_HOMDEL_MAX", 0.5, "SYNTHETIC_DESIGN", "purity-adjusted copy number for homozygous deletion"),
    ("MIN_PURITY_FOR_CN", 0.2, "SYNTHETIC_DESIGN", "below this purity, CN calls are not made"),
    ("RNA_NORM_OFFSET", 10.0, "SYNTHETIC_DESIGN", "housekeeping-centred log2 scale offset"),
    ("PROLIF_HIGH_NLOG2", 7.5, "SYNTHETIC_DESIGN", "mean proliferation gene nlog2 for 'high'"),
    ("LINEAGE_MIN_SCORE", 3.0, "SYNTHETIC_DESIGN", "minimum lineage evidence to assign tissue of origin"),
    ("LINEAGE_MIN_MARGIN", 2.0, "SYNTHETIC_DESIGN", "required margin over runner-up lineage"),
    ("MALIGNANCY_MIN_SCORE", 3.0, "SYNTHETIC_DESIGN", "score >= this: malignant"),
    ("BENIGN_MAX_SCORE", 0.5, "SYNTHETIC_DESIGN", "score <= this (and genomic profiling done): no malignancy detected"),
    ("RNA_STRONG_DELTA", 2.0, "SYNTHETIC_DESIGN", "nlog2 above gene threshold for an RNA call to count as Positive (else Weak)"),
    ("LOW_PURITY_FLAG", 0.2, "SYNTHETIC_DESIGN", "flag specimens with low tumor cellularity"),
]
for n, v, gl, note in PARAMS:
    ind(n, "Parameter", n, paramValue=float(v), hasGuideline=OKG[gl], clinicalNote=note)

# ---------------------------------------------------------------------------
# 4. Lineages and disease hierarchy
# ---------------------------------------------------------------------------
cls("Cancer", "Disease")
nm = cls("NoMalignancyDetected", "Disease", "No malignancy detected",
    "Lineage identified but molecular evidence does not support malignancy (e.g. BPH, benign nodule).")
g.add((nm, OKG.requiresFeature, Literal("fga")))  # benign call only allowed when genomic profiling was performed
cls("IndeterminateDiagnosis", "Disease", "Indeterminate",
    "Insufficient or conflicting evidence; human review required.")

D = [
    # name, parent, label, prevalence (approx, within organ), rare, tier, guideline, therapy
    ("LungCancer", "Cancer", "Lung cancer", None, False, None, None, None),
    ("NonSmallCellLungCancer", "LungCancer", "Non-small cell lung cancer (NSCLC)", "~85% of lung cancer", False, "common", "WHO_THORACIC_2021", None),
    ("LungAdenocarcinoma", "NonSmallCellLungCancer", "Lung adenocarcinoma", "~50% of lung cancer", False, "common", "WHO_THORACIC_2021", None),
    ("LungSquamousCellCarcinoma", "NonSmallCellLungCancer", "Lung squamous cell carcinoma", "~25-30% of lung cancer", False, "common", "WHO_THORACIC_2021", None),
    ("SmallCellLungCancer", "LungCancer", "Small cell lung cancer (SCLC)", "~13-15% of lung cancer", True, "uncommon", "WHO_THORACIC_2021", "platinum-etoposide + IO; DLL3-targeted (tarlatamab)"),
    ("EGFRMutantNSCLC", "NonSmallCellLungCancer", "EGFR-mutant NSCLC", "~15% Western / ~40-50% East Asian LUAD", False, "common", "CAP_IASLC_AMP_2018", "EGFR TKI (osimertinib)"),
    ("KRASG12CMutantNSCLC", "NonSmallCellLungCancer", "KRAS G12C-mutant NSCLC", "~12-13% of LUAD", True, "uncommon", "CAP_IASLC_AMP_2018", "sotorasib / adagrasib"),
    ("ALKRearrangedNSCLC", "NonSmallCellLungCancer", "ALK-rearranged NSCLC", "~3-5% of LUAD", True, "rare", "CAP_IASLC_AMP_2018", "ALK TKI (alectinib, lorlatinib)"),
    ("ROS1RearrangedNSCLC", "NonSmallCellLungCancer", "ROS1-rearranged NSCLC", "~1-2% of LUAD", True, "rare", "CAP_IASLC_AMP_2018", "crizotinib / entrectinib / repotrectinib"),
    ("BRAFV600EMutantNSCLC", "NonSmallCellLungCancer", "BRAF V600E-mutant NSCLC", "~1-2% of LUAD", True, "rare", "CAP_IASLC_AMP_2018", "dabrafenib + trametinib"),
    ("METEx14SkippingNSCLC", "NonSmallCellLungCancer", "MET exon 14 skipping NSCLC", "~3% of LUAD", True, "rare", "CAP_IASLC_AMP_2018", "capmatinib / tepotinib"),
    ("RETRearrangedNSCLC", "NonSmallCellLungCancer", "RET-rearranged NSCLC", "~1-2% of LUAD", True, "rare", "CAP_IASLC_AMP_2018", "selpercatinib / pralsetinib"),
    ("NTRKFusionNSCLC", "NonSmallCellLungCancer", "NTRK fusion-positive NSCLC", "<1% of LUAD", True, "rare", "CAP_IASLC_AMP_2018", "larotrectinib / entrectinib"),

    ("BreastCancer", "Cancer", "Breast cancer", None, False, None, None, None),
    ("HRPositiveHER2NegativeBreastCancer", "BreastCancer", "HR+/HER2- breast cancer", "~70% of breast cancer", False, "common", "ASCO_CAP_ER_PR_2020", "endocrine therapy +/- CDK4/6i"),
    ("LuminalALike", "HRPositiveHER2NegativeBreastCancer", "Luminal A-like", "~40-50%", False, "common", "ST_GALLEN_2013", None),
    ("LuminalBLikeHER2Negative", "HRPositiveHER2NegativeBreastCancer", "Luminal B-like (HER2-negative)", "~20%", False, "common", "ST_GALLEN_2013", None),
    ("HER2PositiveBreastCancer", "BreastCancer", "HER2-positive breast cancer", "~15-20%", False, "common", "ASCO_CAP_HER2_2018", "trastuzumab / pertuzumab / T-DXd"),
    ("LuminalBLikeHER2Positive", "HER2PositiveBreastCancer", "Luminal B-like (HER2-positive)", "~10%", False, "uncommon", "ST_GALLEN_2013", None),
    ("HER2Enriched", "HER2PositiveBreastCancer", "HER2-enriched (HR-/HER2+)", "~5-10%", True, "uncommon", "ST_GALLEN_2013", None),
    ("TripleNegativeBreastCancer", "BreastCancer", "Triple-negative breast cancer (TNBC)", "~10-15%", True, "uncommon", "ASCO_CAP_ER_PR_2020", "chemo-IO (pembrolizumab), sacituzumab govitecan, PARPi if gBRCA"),
    ("TNBCLuminalAndrogenReceptor", "TripleNegativeBreastCancer", "TNBC, luminal androgen receptor (LAR)-like", "~10-15% of TNBC (~1-2% of all BC)", True, "rare", "LEHMANN_TNBC", "AR-directed / PI3K pathway trials"),
    ("ERLowHER2NegativeBreastCancer", "BreastCancer", "ER-low (1-10%), HER2-negative breast cancer", "~2-6%", True, "rare", "ASCO_CAP_ER_PR_2020", "biologically TNBC-like; often treated as TNBC"),

    ("ProstateCancer", "Cancer", "Prostate cancer", None, False, None, None, None),
    ("ProstateAdenocarcinoma", "ProstateCancer", "Prostate adenocarcinoma (acinar)", ">95%", False, "common", "WHO_URINARY_MALE_GENITAL_2022", "ADT / ARPI"),
    ("HRRDeficientProstateCancer", "ProstateAdenocarcinoma", "HRR-deficient prostate adenocarcinoma", "~20-25% of mCRPC (BRCA2 ~10%)", True, "uncommon", "WHO_URINARY_MALE_GENITAL_2022", "PARP inhibitors (olaparib, talazoparib combos)"),
    ("NeuroendocrineProstateCancer", "ProstateCancer", "Neuroendocrine prostate cancer (NEPC)", "<2% de novo; ~10-17% treatment-emergent in mCRPC", True, "rare", "WHO_URINARY_MALE_GENITAL_2022", "platinum-based chemotherapy; trials"),
]
for name, parent, label, prev, rare, tier, gl, tx in D:
    u = cls(name, parent, label, approxPrevalence=prev, isRareSubtype=rare if prev else None,
            rarityTier=tier, therapyRelevance=tx)
    if gl:
        g.add((u, OKG.hasGuideline, OKG[gl]))

# Lineages
LINEAGES = {
    "LungLineage": ("Lung (thoracic) lineage", "LungCancer", None),
    "BreastLineage": ("Breast lineage", "BreastCancer", None),
    "ProstateLineage": ("Prostate lineage", "ProstateCancer", "M"),
}
for k, (lab, root, sex) in LINEAGES.items():
    ind(k, "Lineage", lab, lineageRootDisease=OKG[root], requiredSex=sex)

# ---------------------------------------------------------------------------
# 5. Units and serum biomarkers
# ---------------------------------------------------------------------------
UNITS = {"ng_per_mL": "ng/mL", "ug_per_L": "ug/L", "pg_per_mL": "pg/mL", "ng_per_L": "ng/L",
         "U_per_mL": "U/mL", "kU_per_L": "kU/L", "percent": "%", "fraction": "1", "mut_per_Mb": "/Mb"}
for k, code in UNITS.items():
    ind(k, "Unit", code, ucumCode=code)

UNIT_SYNONYMS = {  # raw spelling -> unit individual
    "ng/mL": "ng_per_mL", "ng/ml": "ng_per_mL", "NG/ML": "ng_per_mL",
    "ug/L": "ug_per_L", "µg/L": "ug_per_L", "mcg/L": "ug_per_L",
    "pg/mL": "pg_per_mL", "pg/ml": "pg_per_mL", "ng/L": "ng_per_L",
    "U/mL": "U_per_mL", "U/ml": "U_per_mL", "kU/L": "kU_per_L", "KU/L": "kU_per_L",
}
for raw, k in UNIT_SYNONYMS.items():
    g.add((OKG[k], SKOS.altLabel, Literal(raw)))

SERUM = [
    # key, label, canonical unit, ULN, synonyms, lineage weights, conversions(from_unit -> factor)
    ("PSA", "Prostate-specific antigen (serum)", "ng_per_mL", 4.0, ["PSA", "Total PSA", "tPSA", "PSA, total"],
     {"ProstateLineage": 1.0}, {"ng_per_mL": 1, "ug_per_L": 1}),
    ("CEA", "Carcinoembryonic antigen", "ng_per_mL", 5.0, ["CEA", "Carcinoembryonic Ag"],
     {}, {"ng_per_mL": 1, "ug_per_L": 1}),
    ("CA15_3", "Cancer antigen 15-3", "U_per_mL", 30.0, ["CA 15-3", "CA15-3", "CA15.3"],
     {"BreastLineage": 0.5}, {"U_per_mL": 1, "kU_per_L": 1}),
    ("CYFRA21_1", "Cytokeratin-19 fragment (CYFRA 21-1)", "ng_per_mL", 3.3, ["CYFRA 21-1", "CYFRA21-1", "Cyfra"],
     {"LungLineage": 0.5}, {"ng_per_mL": 1, "ug_per_L": 1}),
    ("ProGRP", "Pro-gastrin-releasing peptide", "pg_per_mL", 65.0, ["ProGRP", "Pro-GRP"],
     {"LungLineage": 1.0}, {"pg_per_mL": 1, "ng_per_L": 1, "ng_per_mL": 1000}),
    ("NSE", "Neuron-specific enolase", "ng_per_mL", 16.3, ["NSE", "Neuron specific enolase"],
     {}, {"ng_per_mL": 1, "ug_per_L": 1}),
]
for key, label, unit, uln, syns, lw, conv in SERUM:
    u = ind(f"Serum_{key}", "SerumBiomarker", label, featureKey=f"serum.{key}",
            canonicalUnit=OKG[unit], upperLimitNormal=float(uln), hasGuideline=OKG["SYNTHETIC_DESIGN"],
            clinicalNote="ULN is assay-dependent; value used here is a typical reference limit.")
    for s in syns:
        g.add((u, SKOS.altLabel, Literal(s)))
    for fu, f in conv.items():
        c = ind(f"Conv_{key}_{fu}", "UnitConversion", f"{key}: {UNITS[fu]} -> {UNITS[unit]}",
                analyteKey=key, fromUnit=OKG[fu], toUnit=OKG[unit], conversionFactor=float(f))
    for lin, w in lw.items():
        e = ind(f"LE_serum_{key}_{lin}", "LineageEvidence", f"serum {key} > ULN",
                featureKey=f"serum.{key}_xULN", operator="gt", criterionValue="1.0",
                forLineage=OKG[lin], weight=float(w))
        g.add((u, OKG.hasLineageEvidence, e))

# ---------------------------------------------------------------------------
# 6. Genes / protein-expression biomarkers (IHC + RNA)
# ---------------------------------------------------------------------------
# symbol, label, IHC synonyms (raw names), RNA aliases, length bp, rna threshold, lineage weights, role
GENES = [
    ("NKX2-1", "TTF-1 (NKX2-1)", ["TTF-1", "TTF1", "NKX2-1", "Thyroid transcription factor 1"], ["TTF1"], 2300, 6.5, {"LungLineage": 3.0}, "lineage"),
    ("NAPSA", "Napsin A", ["Napsin A", "NAPSA", "napsin-a"], ["NAPSINA"], 1700, 6.5, {"LungLineage": 2.0}, "lineage"),
    ("SFTPC", "Surfactant protein C", [], ["SP-C"], 900, 6.5, {"LungLineage": 1.0}, "lineage"),
    ("TP63", "p40 / p63 (TP63)", ["p40", "TP63", "deltaNp63", "dNp63"], ["P63"], 4600, 6.5, {"LungLineage": 1.5}, "squamous"),
    ("KRT5", "Cytokeratin 5", ["CK5/6", "KRT5"], ["CK5"], 2400, 6.5, {"LungLineage": 0.5}, "squamous"),
    ("SOX2", "SOX2", ["SOX2"], [], 2500, 6.5, {"LungLineage": 0.5}, "squamous"),
    ("GATA3", "GATA3", ["GATA3", "GATA-3"], [], 3300, 6.5, {"BreastLineage": 2.0}, "lineage"),
    ("TRPS1", "TRPS1", ["TRPS1"], [], 9800, 6.5, {"BreastLineage": 3.0}, "lineage"),
    ("SCGB2A2", "Mammaglobin (SCGB2A2)", ["Mammaglobin", "SCGB2A2", "MGB1"], ["MAMMAGLOBIN"], 500, 6.5, {"BreastLineage": 2.0}, "lineage"),
    ("ESR1", "Estrogen receptor (ESR1)", ["ER", "ESR1", "Estrogen receptor", "ER-alpha"], ["ER"], 6300, 6.5, {"BreastLineage": 1.5}, "receptor"),
    ("PGR", "Progesterone receptor (PGR)", ["PR", "PgR", "PGR", "Progesterone receptor"], [], 13000, 6.5, {"BreastLineage": 1.0}, "receptor"),
    ("ERBB2", "HER2 (ERBB2)", ["HER2", "ERBB2", "HER-2/neu", "c-erbB-2"], ["HER2"], 4600, 11.0, {}, "receptor"),
    ("AR", "Androgen receptor", ["AR", "Androgen receptor"], [], 10600, 6.5, {"ProstateLineage": 1.0}, "receptor"),
    ("KLK3", "PSA (KLK3)", ["PSA", "KLK3", "Prostate specific antigen"], ["PSA"], 1500, 6.5, {"ProstateLineage": 3.0}, "lineage"),
    ("NKX3-1", "NKX3.1", ["NKX3.1", "NKX3-1"], ["NKX31"], 3100, 6.5, {"ProstateLineage": 3.0}, "lineage"),
    ("FOLH1", "PSMA (FOLH1)", ["PSMA", "FOLH1"], ["PSMA"], 2800, 6.5, {"ProstateLineage": 1.5}, "lineage"),
    ("AMACR", "AMACR (P504S)", ["AMACR", "P504S"], ["P504S"], 2000, 6.5, {"ProstateLineage": 0.5}, "malignancy"),
    ("HOXB13", "HOXB13", ["HOXB13"], [], 2000, 6.5, {"ProstateLineage": 2.0}, "lineage"),
    ("ERG", "ERG", ["ERG"], [], 5000, 6.5, {}, "fusion-surrogate"),
    ("SYP", "Synaptophysin", ["Synaptophysin", "SYP"], [], 2300, 6.5, {}, "neuroendocrine"),
    ("CHGA", "Chromogranin A", ["Chromogranin A", "CHGA", "CgA"], [], 2000, 6.5, {}, "neuroendocrine"),
    ("INSM1", "INSM1", ["INSM1"], [], 2600, 6.5, {}, "neuroendocrine"),
    ("ASCL1", "ASCL1", [], [], 2400, 6.5, {}, "neuroendocrine"),
    ("MKI67", "Ki-67 (MKI67)", ["Ki-67", "Ki67", "MIB-1", "MKI67"], ["KI67"], 12500, 7.5, {}, "proliferation"),
    ("TOP2A", "TOP2A", [], [], 5600, 7.5, {}, "proliferation"),
    ("BIRC5", "Survivin (BIRC5)", [], [], 1600, 7.5, {}, "proliferation"),
    ("CCNB1", "Cyclin B1", [], [], 2200, 7.5, {}, "proliferation"),
    ("CD274", "PD-L1 (CD274)", ["PD-L1", "PD-L1 (22C3)", "CD274", "PDL1"], ["PDL1"], 3600, 6.5, {}, "immune"),
    ("ALK", "ALK protein", ["ALK", "ALK (D5F3)", "ALK D5F3"], [], 6200, 6.5, {}, "fusion-surrogate"),
]
HOUSEKEEPING = [("ACTB", 1800), ("GAPDH", 1400), ("RPLP0", 1100), ("PGK1", 2400), ("HPRT1", 1400), ("TBP", 1900)]

for sym, label, ihc_syn, rna_alias, length, thr, lw, role in GENES:
    gu = ind(f"Gene_{sym.replace('-', '_')}", "Gene", sym, hgncSymbol=sym, effectiveLengthBp=length,
             rnaHighThreshold=float(thr), clinicalNote=f"role: {role}")
    for a in rna_alias:
        g.add((gu, SKOS.altLabel, Literal(a)))
    bu = ind(f"BM_{sym.replace('-', '_')}", "ProteinExpressionBiomarker", label,
             featureKey=f"status.{sym}", encodedBy=gu, hgncSymbol=sym, clinicalNote=f"role: {role}")
    for s in ihc_syn:
        g.add((bu, OKG.ihcName, Literal(s)))
        g.add((bu, SKOS.altLabel, Literal(s)))
    for lin, w in lw.items():
        e = ind(f"LE_{sym.replace('-', '_')}_{lin}", "LineageEvidence", f"{sym} positive",
                featureKey=f"status.{sym}", operator="in", criterionValue="Positive|Weak",
                forLineage=OKG[lin], weight=float(w))
        g.add((bu, OKG.hasLineageEvidence, e))

for sym, length in HOUSEKEEPING:
    ind(f"Gene_{sym}", ["Gene", "HousekeepingGene"], sym, hgncSymbol=sym, effectiveLengthBp=length)

# Composite (guideline-integrated) biomarkers
for k, lab, gl in [("HER2_status", "HER2 status (IHC + ISH, ASCO/CAP; NCCN null/ultralow/low)", "ASCO_CAP_HER2_2018"),
                   ("ER_status", "ER status (negative / low / positive)", "ASCO_CAP_ER_PR_2020"),
                   ("PR_status", "PR status", "ASCO_CAP_ER_PR_2020"),
                   ("ALK_rearranged", "ALK rearrangement (NGS fusion, FISH >=15% or D5F3 IHC)", "CAP_IASLC_AMP_2018"),
                   ("PDL1_status", "PD-L1 tumor proportion score category", "CAP_IASLC_AMP_2018"),
                   ("NE_count", "Number of positive neuroendocrine markers (SYP, CHGA, INSM1, ASCL1)", "WHO_THORACIC_2021"),
                   ("prolif_high", "Proliferation RNA score high", "SYNTHETIC_DESIGN")]:
    ind(f"CB_{k}", "CompositeBiomarker", lab, featureKey=k, hasGuideline=OKG[gl])

# ---------------------------------------------------------------------------
# 7. Genomic alterations (matchers)
# ---------------------------------------------------------------------------
EX19 = r"^p\.[A-Z](74[5-9]|75[0-9])_[A-Z](74[5-9]|75[0-9])(del|delins[A-Z]+)$"
ALTS = [
    # id, class, label, gene, kwargs
    ("EGFR_SENSITIZING", "SequenceVariantAlteration", "EGFR sensitizing mutation (L858R, exon19 del, G719X, L861Q, S768I)", "EGFR",
     dict(matchProteinRegex=[r"^p\.(L858R|G719[ACSD]|L861Q|S768I)$", EX19], requiresPathogenic=False, isDriver=True)),
    ("EGFR_T790M", "SequenceVariantAlteration", "EGFR T790M", "EGFR", dict(matchProteinRegex=[r"^p\.T790M$"], requiresPathogenic=False, isDriver=False)),
    ("KRAS_G12C", "SequenceVariantAlteration", "KRAS G12C", "KRAS", dict(matchProteinRegex=[r"^p\.G12C$"], requiresPathogenic=False, isDriver=True)),
    ("KRAS_OTHER", "SequenceVariantAlteration", "KRAS non-G12C hotspot", "KRAS", dict(matchProteinRegex=[r"^p\.(G12[ADRSV]|G13[CD]|Q61[HLKR])$"], requiresPathogenic=False, isDriver=True)),
    ("BRAF_V600E", "SequenceVariantAlteration", "BRAF V600E", "BRAF", dict(matchProteinRegex=[r"^p\.V600E$"], requiresPathogenic=False, isDriver=True)),
    ("MET_EX14", "SequenceVariantAlteration", "MET exon 14 skipping", "MET",
     dict(matchProteinRegex=[r"^exon14skip$"], matchCdnaRegex=[r"^c\.(3028|3082)\+[1-3]"], requiresPathogenic=False, isDriver=True)),
    ("PIK3CA_HOTSPOT", "SequenceVariantAlteration", "PIK3CA hotspot (H1047R/L, E545K, E542K)", "PIK3CA",
     dict(matchProteinRegex=[r"^p\.(H1047[RL]|E545K|E542K)$"], requiresPathogenic=False, isDriver=True)),
    ("ESR1_LBD", "SequenceVariantAlteration", "ESR1 ligand-binding-domain mutation", "ESR1",
     dict(matchProteinRegex=[r"^p\.(Y537[SNC]|D538G|E380Q)$"], requiresPathogenic=False, isDriver=False)),
    ("SPOP_MUT", "SequenceVariantAlteration", "SPOP MATH-domain mutation", "SPOP",
     dict(matchProteinRegex=[r"^p\.(F133[VLC]|W131[GC]|Y87[CN]|F102C)$"], requiresPathogenic=False, isDriver=True)),
]
for gene in ["TP53", "BRCA1", "BRCA2", "ATM", "CDK12", "PALB2", "STK11", "KEAP1"]:
    ALTS.append((f"{gene}_PATHOGENIC", "SequenceVariantAlteration", f"{gene} pathogenic / likely pathogenic variant", gene,
                 dict(requiresPathogenic=True, isDriver=gene in ("TP53",))))
for gene in ["RB1", "PTEN"]:
    ALTS.append((f"{gene}_TRUNC", "SequenceVariantAlteration", f"{gene} pathogenic variant", gene, dict(requiresPathogenic=True, isDriver=True)))
    ALTS.append((f"{gene}_HOMDEL", "CopyNumberAlteration", f"{gene} homozygous deletion", gene, dict(cnaType="homdel", isDriver=True)))
for gene in ["ERBB2", "MYCN", "AURKA"]:
    ALTS.append((f"{gene}_AMP", "CopyNumberAlteration", f"{gene} amplification", gene, dict(cnaType="amp", isDriver=True)))
for alt_id, three, five, lab in [("ALK_FUSION", "ALK", None, "ALK fusion (any 5' partner)"),
                                 ("ROS1_FUSION", "ROS1", None, "ROS1 fusion"), ("RET_FUSION", "RET", None, "RET fusion"),
                                 ("NTRK1_FUSION", "NTRK1", None, "NTRK1 fusion"), ("NTRK2_FUSION", "NTRK2", None, "NTRK2 fusion"),
                                 ("NTRK3_FUSION", "NTRK3", None, "NTRK3 fusion"),
                                 ("TMPRSS2_ERG", "ERG", "TMPRSS2", "TMPRSS2::ERG fusion")]:
    ALTS.append((alt_id, "GeneFusionAlteration", lab, three, dict(fusion3Prime=three, fusion5Prime=five, isDriver=True)))

for alt_id, klass, label, gene, kw in ALTS:
    ind(f"ALT_{alt_id}", klass, label, featureKey=f"alt.{alt_id}", matchGene=gene, **kw)

# Gene symbol aliases appearing in raw genomic files
GENE_ALIASES = {"ERBB2": ["HER2", "HER-2", "NEU"], "MET": ["c-MET", "HGFR"], "NKX2-1": ["TTF1"],
                "CDKN2A": ["p16"], "TMPRSS2": ["TMPRSS-2"]}
for sym, al in GENE_ALIASES.items():
    gu = OKG[f"Gene_{sym.replace('-', '_')}"]
    g.add((gu, RDF.type, OKG.Gene)); g.add((gu, OKG.hgncSymbol, Literal(sym)))
    for a in al:
        g.add((gu, SKOS.altLabel, Literal(a)))

# ---------------------------------------------------------------------------
# 8. Clinical lineage evidence + specimen-site priors
# ---------------------------------------------------------------------------
for name, fk, op, val, lin, w, lab in [
    ("LE_tmprss2erg_prostate", "alt.TMPRSS2_ERG", "eq", "true", "ProstateLineage", 3.0, "TMPRSS2::ERG fusion (prostate-specific)"),
    ("LE_history_prostate", "history.prostate_cancer", "eq", "true", "ProstateLineage", 3.0, "history of prostate cancer"),
    ("LE_history_breast", "history.breast_cancer", "eq", "true", "BreastLineage", 3.0, "history of breast cancer"),
    ("LE_smoking_lung", "clin.pack_years", "ge", "20", "LungLineage", 0.5, "smoking >= 20 pack-years"),
    ("LE_female_breast", "clin.sex", "eq", "F", "BreastLineage", 0.5, "female sex"),
]:
    ind(name, "LineageEvidence", lab, featureKey=fk, operator=op, criterionValue=val, forLineage=OKG[lin], weight=w)

for kw, lin in [("lung", "LungLineage"), ("bronch", "LungLineage"), ("breast", "BreastLineage"),
                ("prostat", "ProstateLineage")]:
    ind(f"SSP_{kw}", "SpecimenSitePrior", f"specimen site '{kw}*'", siteKeyword=kw, forLineage=OKG[lin], weight=1.5)

# ---------------------------------------------------------------------------
# 9. Malignancy evidence
# ---------------------------------------------------------------------------
MAL = [
    ("ME_driver1", "n_somatic_drivers", "ge", "1", 2.0, "somatic pathogenic/hotspot alteration (VAF >= floor)"),
    ("ME_driver2", "n_somatic_drivers", "ge", "2", 1.0, "two or more somatic alterations"),
    ("ME_fusion", "n_driver_fusions", "ge", "1", 3.0, "driver gene fusion"),
    ("ME_cna", "n_focal_cna", "ge", "1", 2.0, "focal amplification / homozygous deletion"),
    ("ME_fga", "fga", "ge", "0.1", 2.0, "fraction genome altered >= 10% (aneuploidy)"),
    ("ME_ki67_10", "Ki67_pct", "ge", "10", 1.0, "Ki-67 >= 10%"),
    ("ME_ki67_20", "Ki67_pct", "ge", "20", 1.0, "Ki-67 >= 20%"),
    ("ME_prolif", "prolif_high", "eq", "true", 1.0, "high proliferation RNA score"),
    ("ME_amacr", "status.AMACR", "in", "Positive", 1.0, "AMACR (P504S) positive"),
    ("ME_ne", "NE_count", "ge", "2", 1.0, "neuroendocrine program"),
    ("ME_tmb", "tmb", "ge", "5", 1.0, "TMB >= 5 mut/Mb"),
    ("ME_tmb_low", "tmb", "ge", "1.5", 1.0, "TMB >= 1.5 mut/Mb (above normal-tissue background)"),
    ("ME_fga_low", "fga", "ge", "0.05", 1.0, "fraction genome altered >= 5%"),
]
for n, fk, op, v, w, lab in MAL:
    ind(n, "MalignancyEvidence", lab, featureKey=fk, operator=op, criterionValue=v, weight=w)

# ---------------------------------------------------------------------------
# 10. Diagnostic rules (reified criteria)
# ---------------------------------------------------------------------------
_crit_n = [0]


def crit(fk, op, val=None, missing="fail", any_of=None, label=None):
    _crit_n[0] += 1
    u = OKG[f"Crit_{_crit_n[0]:03d}"]
    g.add((u, RDF.type, OWL.NamedIndividual)); g.add((u, RDF.type, OKG.Criterion))
    if any_of:
        g.add((u, OKG.operator, Literal("anyOf")))
        for c in any_of:
            g.add((u, OKG.anyOf, c))
        g.add((u, RDFS.label, Literal(label or "any of: " + " | ".join(str(g.value(c, RDFS.label)) for c in any_of))))
        return u
    g.add((u, OKG.featureKey, Literal(fk))); g.add((u, OKG.operator, Literal(op)))
    if val is not None:
        g.add((u, OKG.criterionValue, Literal(str(val))))
    g.add((u, OKG.missingIs, Literal(missing)))
    g.add((u, RDFS.label, Literal(label or f"{fk} {op} {val}")))
    return u


HER2_NEG = "Null|Ultralow|Low|Negative"
GROUPS = [
    # id, applies to, order, slot, label
    ("RG_lung_histology", "LungCancer", 1, "histology", "Lung histologic type (IHC/NE program)"),
    ("RG_nsclc_driver", "NonSmallCellLungCancer", 2, "molecular", "NSCLC oncogenic driver (mutually exclusive, rarest first)"),
    ("RG_breast_receptor", "BreastCancer", 1, "histology", "Breast receptor-defined surrogate subtype"),
    ("RG_tnbc_subtype", "TripleNegativeBreastCancer", 2, "molecular", "TNBC molecular subtype"),
    ("RG_prostate_histology", "ProstateCancer", 1, "histology", "Prostate histologic program"),
    ("RG_prostate_molecular", "ProstateAdenocarcinoma", 2, "molecular", "Prostate actionable molecular class"),
]
for gid, applies, order, slot, lab in GROUPS:
    ind(gid, "RuleGroup", lab, groupAppliesTo=OKG[applies], groupOrder=order, outputSlot=slot)

RULES = [
    # id, group, priority, assigns, guideline, criteria builder
    ("R_SCLC", "RG_lung_histology", 1, "SmallCellLungCancer", "WHO_THORACIC_2021",
     lambda: [crit("NE_count", "ge", 2, label=">=2 neuroendocrine markers positive"),
              crit(None, None, any_of=[crit("alt.RB1_ANY", "eq", "true", label="RB1 loss"),
                                       crit("Ki67_pct", "ge", 50, label="Ki-67 >= 50%"),
                                       crit("prolif_high", "eq", "true", label="proliferation RNA high")])]),
    ("R_LUSC", "RG_lung_histology", 2, "LungSquamousCellCarcinoma", "WHO_THORACIC_2021",
     lambda: [crit("status.TP63", "in", "Positive", label="p40 (TP63) positive"),
              crit("status.NKX2-1", "nin", "Positive|Weak", missing="pass", label="TTF-1 not positive")]),
    ("R_LUAD", "RG_lung_histology", 3, "LungAdenocarcinoma", "WHO_THORACIC_2021",
     lambda: [crit(None, None, any_of=[crit("status.NKX2-1", "in", "Positive|Weak", label="TTF-1 positive"),
                                       crit("status.NAPSA", "in", "Positive|Weak", label="Napsin A positive")]),
              crit("status.TP63", "nin", "Positive", missing="pass", label="p40 not diffusely positive"),
              crit("NE_count", "lt", 2, missing="pass", label="no neuroendocrine program")]),
    ("R_NSCLC_NOS", "RG_lung_histology", 9, "NonSmallCellLungCancer", "WHO_THORACIC_2021",
     lambda: [crit("NE_count", "lt", 2, missing="pass", label="no neuroendocrine program")]),

    ("R_ALK", "RG_nsclc_driver", 1, "ALKRearrangedNSCLC", "CAP_IASLC_AMP_2018", lambda: [crit("ALK_rearranged", "eq", "true", label="ALK rearranged")]),
    ("R_ROS1", "RG_nsclc_driver", 2, "ROS1RearrangedNSCLC", "CAP_IASLC_AMP_2018", lambda: [crit("alt.ROS1_FUSION", "eq", "true")]),
    ("R_NTRK", "RG_nsclc_driver", 3, "NTRKFusionNSCLC", "CAP_IASLC_AMP_2018", lambda: [crit("alt.NTRK_ANY", "eq", "true", label="NTRK1/2/3 fusion")]),
    ("R_RET", "RG_nsclc_driver", 4, "RETRearrangedNSCLC", "CAP_IASLC_AMP_2018", lambda: [crit("alt.RET_FUSION", "eq", "true")]),
    ("R_MET", "RG_nsclc_driver", 5, "METEx14SkippingNSCLC", "CAP_IASLC_AMP_2018", lambda: [crit("alt.MET_EX14", "eq", "true")]),
    ("R_BRAF", "RG_nsclc_driver", 6, "BRAFV600EMutantNSCLC", "CAP_IASLC_AMP_2018", lambda: [crit("alt.BRAF_V600E", "eq", "true")]),
    ("R_EGFR", "RG_nsclc_driver", 7, "EGFRMutantNSCLC", "CAP_IASLC_AMP_2018", lambda: [crit("alt.EGFR_SENSITIZING", "eq", "true")]),
    ("R_KRASG12C", "RG_nsclc_driver", 8, "KRASG12CMutantNSCLC", "CAP_IASLC_AMP_2018", lambda: [crit("alt.KRAS_G12C", "eq", "true")]),

    ("R_HER2E", "RG_breast_receptor", 1, "HER2Enriched", "ASCO_CAP_HER2_2018",
     lambda: [crit("HER2_status", "in", "Positive"), crit("ER_status", "in", "Negative"), crit("PR_status", "in", "Negative")]),
    ("R_LUMB_HER2P", "RG_breast_receptor", 2, "LuminalBLikeHER2Positive", "ASCO_CAP_HER2_2018",
     lambda: [crit("HER2_status", "in", "Positive"),
              crit(None, None, any_of=[crit("ER_status", "in", "Positive|Low"), crit("PR_status", "in", "Positive")])]),
    ("R_TNBC", "RG_breast_receptor", 3, "TripleNegativeBreastCancer", "ASCO_CAP_ER_PR_2020",
     lambda: [crit("ER_status", "in", "Negative", label="ER < 1%"), crit("PR_status", "in", "Negative", label="PR < 1%"),
              crit("HER2_status", "in", HER2_NEG, label="HER2-negative (IHC 0/1+ or 2+/ISH-)")]),
    ("R_ERLOW", "RG_breast_receptor", 4, "ERLowHER2NegativeBreastCancer", "ASCO_CAP_ER_PR_2020",
     lambda: [crit("ER_status", "in", "Low", label="ER 1-10%"), crit("HER2_status", "in", HER2_NEG)]),
    ("R_LUMB_HER2N", "RG_breast_receptor", 5, "LuminalBLikeHER2Negative", "ST_GALLEN_2013",
     lambda: [crit("ER_status", "in", "Positive", label="ER > 10%"), crit("HER2_status", "in", HER2_NEG),
              crit(None, None, any_of=[crit("Ki67_pct", "ge", 20, label="Ki-67 >= 20%"),
                                       crit("PR_pct", "lt", 20, label="PR < 20%")])]),
    ("R_LUMA", "RG_breast_receptor", 6, "LuminalALike", "ST_GALLEN_2013",
     lambda: [crit("ER_status", "in", "Positive"), crit("HER2_status", "in", HER2_NEG),
              crit("Ki67_pct", "lt", 20, label="Ki-67 < 20%"), crit("PR_pct", "ge", 20, label="PR >= 20%")]),
    ("R_HRPOS_NOS", "RG_breast_receptor", 7, "HRPositiveHER2NegativeBreastCancer", "ASCO_CAP_ER_PR_2020",
     lambda: [crit("ER_status", "in", "Positive"), crit("HER2_status", "in", HER2_NEG)]),

    ("R_TNBC_LAR", "RG_tnbc_subtype", 1, "TNBCLuminalAndrogenReceptor", "LEHMANN_TNBC",
     lambda: [crit("AR_pct", "ge", 10, label="AR IHC >= 10%")]),

    ("R_NEPC", "RG_prostate_histology", 1, "NeuroendocrineProstateCancer", "WHO_URINARY_MALE_GENITAL_2022",
     lambda: [crit("NE_count", "ge", 2, label=">=2 neuroendocrine markers positive"),
              crit("status.AR", "nin", "Positive", missing="pass", label="AR not positive"),
              crit("status.KLK3", "nin", "Positive", missing="pass", label="PSA (tissue) not positive")]),
    ("R_PRAD", "RG_prostate_histology", 2, "ProstateAdenocarcinoma", "WHO_URINARY_MALE_GENITAL_2022", lambda: []),
    ("R_HRR", "RG_prostate_molecular", 1, "HRRDeficientProstateCancer", "WHO_URINARY_MALE_GENITAL_2022",
     lambda: [crit("alt.HRR_ANY", "eq", "true", label="BRCA1/2, ATM, CDK12 or PALB2 pathogenic")]),
]
for rid, grp, prio, assigns, gl, builder in RULES:
    u = ind(rid, "DiagnosticRule", f"{rid} -> {assigns}", inRuleGroup=OKG[grp], priority=prio,
            ruleAssigns=OKG[assigns], hasGuideline=OKG[gl])
    for c in builder():
        g.add((u, OKG.hasCriterion, c))

# ---------------------------------------------------------------------------
# 11. Molecular modifiers (non-exclusive annotations)
# ---------------------------------------------------------------------------
MODS = [
    ("MOD_HER2_LOW", "HER2-low (IHC 1+ or 2+/ISH-)", "BreastCancer", [("HER2_status", "in", "Low")], "T-DXd eligible (metastatic)"),
    ("MOD_HER2_ULTRALOW", "HER2-ultralow (IHC 0, faint staining <=10%)", "BreastCancer", [("HER2_status", "in", "Ultralow")], "T-DXd eligible (HR+, post-endocrine)"),
    ("MOD_HER2_NULL", "HER2-null (IHC 0, no staining)", "BreastCancer", [("HER2_status", "in", "Null")], None),
    ("MOD_ER_LOW", "ER-low positive (1-10%)", "BreastCancer", [("ER_status", "in", "Low")], None),
    ("MOD_BRCA1", "BRCA1 pathogenic variant", "Cancer", [("alt.BRCA1_PATHOGENIC", "eq", "true")], "PARP inhibitor; germline testing"),
    ("MOD_BRCA2", "BRCA2 pathogenic variant", "Cancer", [("alt.BRCA2_PATHOGENIC", "eq", "true")], "PARP inhibitor; germline testing"),
    ("MOD_PIK3CA", "PIK3CA hotspot mutation", "BreastCancer", [("alt.PIK3CA_HOTSPOT", "eq", "true")], "alpelisib / capivasertib / inavolisib"),
    ("MOD_ESR1", "ESR1 LBD mutation (endocrine resistance)", "BreastCancer", [("alt.ESR1_LBD", "eq", "true")], "elacestrant"),
    ("MOD_TP53", "TP53 pathogenic variant", "Cancer", [("alt.TP53_PATHOGENIC", "eq", "true")], None),
    ("MOD_PDL1_HIGH", "PD-L1 TPS >= 50%", "LungCancer", [("PDL1_status", "in", "High")], "IO monotherapy eligible (no driver)"),
    ("MOD_PDL1_LOW", "PD-L1 TPS 1-49%", "LungCancer", [("PDL1_status", "in", "Low")], None),
    ("MOD_PDL1_NEG", "PD-L1 TPS < 1%", "LungCancer", [("PDL1_status", "in", "Negative")], None),
    ("MOD_EGFR_T790M", "EGFR T790M", "LungCancer", [("alt.EGFR_T790M", "eq", "true")], "osimertinib"),
    ("MOD_KRAS_OTHER", "KRAS non-G12C mutation", "LungCancer", [("alt.KRAS_OTHER", "eq", "true")], None),
    ("MOD_STK11", "STK11 loss", "LungCancer", [("alt.STK11_PATHOGENIC", "eq", "true")], "may blunt IO response"),
    ("MOD_KEAP1", "KEAP1 loss", "LungCancer", [("alt.KEAP1_PATHOGENIC", "eq", "true")], None),
    ("MOD_RB1", "RB1 loss", "Cancer", [("alt.RB1_ANY", "eq", "true")], None),
    ("MOD_TMPRSS2_ERG", "TMPRSS2::ERG fusion", "ProstateCancer", [("alt.TMPRSS2_ERG", "eq", "true")], None),
    ("MOD_PTEN", "PTEN loss", "ProstateCancer", [("alt.PTEN_ANY", "eq", "true")], "capivasertib (PTEN-deficient mHSPC)"),
    ("MOD_SPOP", "SPOP mutation", "ProstateCancer", [("alt.SPOP_MUT", "eq", "true")], None),
    ("MOD_MYCN", "MYCN amplification", "ProstateCancer", [("alt.MYCN_AMP", "eq", "true")], None),
    ("MOD_AURKA", "AURKA amplification", "ProstateCancer", [("alt.AURKA_AMP", "eq", "true")], None),
    ("MOD_ATM", "ATM pathogenic variant", "ProstateCancer", [("alt.ATM_PATHOGENIC", "eq", "true")], "PARP inhibitor (lower benefit than BRCA)"),
    ("MOD_CDK12", "CDK12 biallelic loss", "ProstateCancer", [("alt.CDK12_PATHOGENIC", "eq", "true")], None),
]
for mid, lab, applies, crits, tx in MODS:
    u = ind(mid, "MolecularModifier", lab, modifierAppliesTo=OKG[applies], therapyRelevance=tx)
    for fk, op, v in crits:
        g.add((u, OKG.hasCriterion, crit(fk, op, v)))

# Aggregate alteration features computed by the normalizer (documented here)
for k, lab in [("alt.RB1_ANY", "RB1 pathogenic variant OR homozygous deletion"),
               ("alt.PTEN_ANY", "PTEN pathogenic variant OR homozygous deletion"),
               ("alt.NTRK_ANY", "any NTRK1/2/3 fusion"),
               ("alt.HRR_ANY", "BRCA1/BRCA2/ATM/CDK12/PALB2 pathogenic variant")]:
    ind("AGG_" + k.split(".")[1], "CompositeBiomarker", lab, featureKey=k)

out = Path(__file__).with_name("oncokg.ttl")
g.serialize(out, format="turtle")
print(f"wrote {out} with {len(g)} triples")
n_cls = len(set(g.subjects(RDF.type, OWL.Class)))
print(f"classes={n_cls} rules={len(RULES)} modifiers={len(MODS)} alterations={len(ALTS)}")
