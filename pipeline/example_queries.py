"""Example SPARQL queries over the patient knowledge graph (ontology + inferred patient facts)."""
from rdflib import Graph
from common import ROOT

g = Graph()
g.parse(ROOT / "ontology" / "oncokg.ttl")
g.parse(ROOT / "output" / "patient_kg.ttl")
P = """PREFIX okg: <http://example.org/oncokg#>
PREFIX rdfs: <http://www.w3.org/2000/01/rdf-schema#>
"""
QUERIES = {
    "TNBC patients with a BRCA1/2 pathogenic variant (PARP-inhibitor candidates)": """
SELECT ?patient ?mod WHERE {
  ?patient okg:hasTumor ?t . ?t a okg:TripleNegativeBreastCancer ; okg:hasModifier ?m .
  ?m rdfs:label ?mod . FILTER(?m IN (okg:MOD_BRCA1, okg:MOD_BRCA2)) } ORDER BY ?patient""",
    "All rare-subtype diagnoses by class (subclass-aware via rdfs:subClassOf*)": """
SELECT ?label (COUNT(DISTINCT ?patient) AS ?n) WHERE {
  ?patient okg:hasTumor ?t . ?t a ?c . ?c okg:isRareSubtype true ; rdfs:label ?label } GROUP BY ?label ORDER BY DESC(?n)""",
    "Every lung cancer patient (any subtype) using the class hierarchy": """
SELECT (COUNT(DISTINCT ?patient) AS ?n) WHERE {
  ?patient okg:hasTumor ?t . ?t a ?c . ?c rdfs:subClassOf* okg:LungCancer }""",
    "Diagnoses needing human review (low confidence)": """
SELECT ?patient ?label WHERE { ?patient okg:hasDiagnosis ?d . ?d okg:confidence "Low" ; rdfs:label ?label } LIMIT 10""",
    "Which findings supported the ALK-rearranged calls": """
SELECT ?patient ?feature ?value WHERE {
  ?patient okg:hasDiagnosis ?d . ?d okg:derivedByRule okg:R_ALK ; okg:supportedByFinding ?f .
  ?f okg:featureKey ?feature ; okg:normalizedValue ?value . FILTER(?feature = "ALK_rearranged") } LIMIT 5""",
}
for title, sparql in QUERIES.items():
    print(f"\n## {title}")
    for row in g.query(P + sparql):
        print("  " + " | ".join(str(x).split("/")[-1].split("#")[-1] for x in row))
