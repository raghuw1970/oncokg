# OncoKG — Explainable Oncology Reasoning Demo

**[▶ Live demo](https://raghuw1970.github.io/oncokg/)** · synthetic data · not for clinical use

**Problem.** Clinical and molecular test data arrives in inconsistent formats, while oncology classification depends on rules, thresholds, synonyms, and subtype logic that change over time.

**What I built.** I originally built this capability within a commercial oncology data platform in partnership with bioinformatics SMEs. Because the production system is proprietary, I rebuilt OncoKG as a live, non-proprietary demonstration using synthetic data.

**How it works.** Ontology → normalization → reasoning → evidence trace → human review.

**My role.** I translated SME knowledge into the product model, ontology structure, normalization logic, reasoning flow, explainability model, and interactive explorer.

**Why it matters.** The architecture separates domain knowledge from application code, making rules easier to extend, test, audit, and govern.

---

## Technical overview

A working, end-to-end prototype:

1. an **OWL ontology** (Turtle) that encodes cancer classes, biomarkers, thresholds and diagnostic rules,
2. a **synthetic raw molecular dataset** from three labs with deliberately inconsistent formats,
3. a **normalization pipeline** that turns raw exports into canonical, guideline-interpreted findings, and
4. an **ontology-driven reasoner** that identifies tissue of origin, malignancy, subtype (including rare types such as triple-negative breast cancer) and actionable modifiers, with a full evidence trail and a patient knowledge graph.

> Research and education only. Thresholds follow published guidelines where cited, but this is not a validated clinical decision tool, and all patients are synthetic.

## Run it

```bash
pip install -r requirements.txt
./run_all.sh                       # ontology -> data -> normalize -> reason -> evaluate -> explorer
python3 pipeline/example_queries.py  # SPARQL over the patient knowledge graph
```

## Layout

| Path | What it is |
|---|---|
| `ontology/build_ontology.py` → `oncokg.ttl` | Ontology source and output (about 3,350 triples) |
| `data/generate_synthetic.py` | Synthetic cohort generator (598 patients, 3 labs) |
| `data/raw/*.csv` | Raw lab exports (input to the pipeline) |
| `data/ground_truth.csv` | Hidden labels, used **only** by `evaluate.py` |
| `pipeline/normalize.py` | Raw → normalized findings (all mappings and thresholds read from the ontology) |
| `pipeline/reason.py` | Rules from the ontology via SPARQL → diagnoses, traces, patient KG |
| `pipeline/evaluate.py` | Accuracy against ground truth |
| `pipeline/build_explorer.py` | Builds the interactive HTML explorer |
| `docs/index.html` | Published demo (GitHub Pages serves `/docs`) |
| `output/` | `normalized_findings.csv` (long, with provenance), `normalized_features.csv` (wide), `patient_diagnoses.csv`, `evidence_trace.json`, `patient_kg.ttl`, `evaluation.json`, `oncokg_explorer.html` |

## 1. Ontology (`okg:` = `http://example.org/oncokg#`)

**The ontology is the single source of truth.** Neither the normalizer nor the reasoner hard-codes a synonym, unit, cut-off or rule.

| Component | Contents |
|---|---|
| Disease hierarchy | `Cancer` → `LungCancer` / `BreastCancer` / `ProstateCancer`, each with subtypes annotated with `approxPrevalence`, `isRareSubtype`, `rarityTier` and `therapyRelevance`. Also `NoMalignancyDetected` and `IndeterminateDiagnosis`. |
| Lung | NSCLC → adenocarcinoma, squamous; SCLC. Driver classes: EGFR, KRAS G12C, ALK, ROS1, BRAF V600E, MET exon 14, RET, NTRK. |
| Breast | Luminal A-like, Luminal B-like (HER2−/HER2+), HER2-enriched, **TNBC**, **TNBC-LAR**, **ER-low (1–10%)**, HR+/HER2− (when Luminal A vs B can't be decided). |
| Prostate | Adenocarcinoma, **HRR-deficient** adenocarcinoma, **neuroendocrine prostate cancer (NEPC)**. |
| Biomarkers | 29 protein/RNA markers (IHC names and RNA aliases as `skos:altLabel`), 6 serum markers (canonical unit, upper limit of normal (ULN), unit conversions), 6 housekeeping genes, 31 genomic alteration definitions (HGVS regex, fusion 5′/3′ partners, amplification/deletion). |
| Composite biomarkers | HER2 status (ASCO/CAP 2018, affirmed 2023, plus the NCCN 2026 null/ultralow/low split), ER/PR status (ASCO/CAP 2020), ALK rearrangement, PD-L1 TPS category, neuroendocrine score. |
| Reasoning knowledge | `LineageEvidence` (weighted marker → tissue of origin), `SpecimenSitePrior`, `MalignancyEvidence`, `RuleGroup` / `DiagnosticRule` / `Criterion` (reified, with `anyOf`, `missingIs`), 24 `MolecularModifier`s, 27 `Parameter`s, each with a `Guideline` source. |

### How rare types are identified
* **TNBC**: ER < 1% **and** PR < 1% **and** HER2 ∈ {null, ultralow, low} (IHC 0/1+, or 2+ with negative ISH). HER2 IHC 2+ without ISH is held as *equivocal*: the reasoner refuses to call TNBC and asks for reflex ISH.
* **TNBC-LAR**: TNBC plus AR ≥ 10%.
* **ER-low**: ER 1–10%, HER2-negative. This is kept separate from TNBC because guidelines report it as ER-positive, while its biology is TNBC-like.
* **NEPC**: two or more neuroendocrine markers, AR and tissue PSA not positive, and prostate lineage. Prostate lineage can still be established after AR/PSA loss from a TMPRSS2::ERG fusion or a history of prostate cancer.
* **Rare NSCLC drivers** are evaluated rarest-first in one mutually exclusive group; any co-occurrence is flagged.

## 2. Synthetic raw data (input)

598 patients: 210 lung, 228 breast, 120 prostate and 40 benign controls (benign lung nodule, benign breast, BPH with raised PSA). They come from three labs whose conventions differ:

| Data | SITE_A | SITE_B | SITE_C |
|---|---|---|---|
| ER/PR | `85%` | `Allred PS5 IS3 (TS 8)` | `Positive, 90%, strong` |
| Ki-67 | `25%` | `MIB-1 index 0.25` | `Ki-67: 25% (hot spot)` |
| HER2 IHC | `2+` | `Score 2+ (equivocal)` | `HER2 IHC 0 - faint/… incomplete staining in 8%` |
| Variants | `p.L858R`, VAF `0.231` | `p.Leu858Arg`, VAF `23.1%` | `L858R`, VAF `23.1` |
| Fusions | `EML4::ALK` | `EML4-ALK` | `EML4/ALK (e13:e20)` |
| Copy number | log2 ratio | observed CN | fold change |
| Serum | ng/mL | µg/L, kU/L, ng/L | ProGRP in ng/mL |
| RNA counts | gene symbols | gene symbols | aliases (`TTF1`, `PSA`, `HER-2`…) |

Built-in noise: about 8% of IHC stains missing; RNA available for about 70% of patients and serum for about 60%; 5% of specimens with low tumour purity; CHIP-like low-VAF variants; VUS; low-depth calls; 1–4-read fusion artefacts; 2% IHC/RNA discordance. Designed edge cases: HER2 2+ awaiting ISH, ER−/PR+, TTF-1/Napsin A double-negative adenocarcinoma, near-diploid luminal A, male breast cancer, and metastatic specimens.

## 3. Normalization

| Step | Method (parameters from the ontology) |
|---|---|
| IHC | Stain synonym → HGNC symbol. %, Allred proportion score (→ % range), fraction and free text → %. ER: <1% neg / 1–10% low / >10% pos. PR ≥ 1% pos. AR ≥ 10% pos. PD-L1 TPS <1 / 1–49 / ≥50. Qualitative stains → Positive / Weak / Negative. |
| HER2 | ISH ratio and copies → ASCO/CAP groups 1–5, combined with the IHC score → Positive / Low / Ultralow / Null / Equivocal. |
| Variants | Gene alias → HGNC symbol; 3-letter → 1-letter HGVS; textual `METex14`; VAF scale inferred per lab; ACMG/AMP class harmonised; QC: depth ≥ 100, VAF ≥ 5%; matched against alteration definitions. |
| Fusions | Separators harmonised; ≥ 5 supporting reads; a fusion missing from a lab that never ran the panel counts as **unknown**, not negative. |
| Copy number | Converted to a ratio, then purity-adjusted CN = (2r − 2(1−p))/p. Amplification ≥ 6, homozygous deletion ≤ 0.5; no calls when purity is below 20%. |
| RNA | counts / length (kb) → log2(+1) → centred on housekeeping genes (+10). Per-gene calls: High / intermediate / Low. Proliferation score. |
| Serum | Analyte synonyms, unit conversion, values below detection limit → LOD/2, × ULN. |
| Consolidation | IHC preferred; RNA as fallback. ER/PR from RNA only when IHC is missing, and always flagged. Discordances flagged. |

## 4. Reasoning

1. **Tissue of origin**: weighted lineage evidence plus specimen-site prior. Sex constraint (prostate requires male). Assigned at score ≥ 3 with a margin ≥ 2; otherwise unresolved.
2. **Malignancy**: weighted evidence (drivers, fusions, focal CNAs, aneuploidy, TMB, Ki-67, proliferation, AMACR, neuroendocrine program). **A benign call requires genomic profiling to have been done**; otherwise the system abstains.
3. **Subtype**: rule groups in order; within a group, the first rule whose criteria all hold wins. Missing data or an Allred range straddling a cut-off is treated as *uncertain*, never as satisfied.
4. **Modifiers**, **therapy relevance**, **rarity**, **confidence** and **flags** for human review.
5. **Patient knowledge graph** (`patient_kg.ttl`): each patient is linked to a specimen, findings, a tumour typed with ontology classes, modifiers, and a `DiagnosisAssertion` with `derivedByRule` and `supportedByFinding` (only the findings the fired rules used).

## 5. Results on the synthetic cohort

| Metric | Result |
|---|---|
| Tissue of origin | 587/598 correct, **0 wrong**, 11 unresolved |
| Subtype | 527 exact, 27 correct but less specific, 43 abstained, **1 wrong** |
| Benign controls | 37/40 called benign, 3 abstained, 0 called cancer; **0 cancers called benign** |
| Rare subtypes | 197/215 found, PPV 99.5% |
| **TNBC** (incl. LAR) | 44/45 found, PPV 97.8% |
| Modifiers | precision 100% for every modifier; recall 50–100% (limited by missing tests, purity and QC) |

The single wrong call was an ER−/PR+ tumour whose PR IHC was missing. PR was inferred from RNA, and that RNA value was a simulated discordant read, so the case came out as TNBC. The call carried a "confirm PR by IHC" flag.

**Caveat:** the rules and the generator were designed together, so these numbers are optimistic. Real-world performance needs validation on real cohorts (e.g. TCGA/CPTAC with IHC annotations, or institutional data).

## Extending it
* Add a subtype: add a class, then a `DiagnosticRule` with `Criterion` nodes, in `build_ontology.py`. No pipeline code changes are needed.
* Add a lab format: add `skos:altLabel` synonyms or a `UnitConversion` individual.
* Other options: swap the reified rules for SWRL/SHACL if you need an OWL reasoner, map classes to NCIt/OncoTree/HGNC IRIs, or load `patient_kg.ttl` into GraphDB/Neo4j (n10s).

## License
MIT (see `LICENSE`). All data in this repository is synthetic.
