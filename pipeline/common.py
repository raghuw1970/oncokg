"""Shared helpers: load the OncoKG ontology and expose its contents to the pipeline via SPARQL."""
from pathlib import Path
from functools import lru_cache
from rdflib import Graph, Namespace, RDF, RDFS
from rdflib.namespace import SKOS

ROOT = Path(__file__).resolve().parents[1]
OKG = Namespace("http://example.org/oncokg#")
PREFIX = """PREFIX okg: <http://example.org/oncokg#>
PREFIX rdfs: <http://www.w3.org/2000/01/rdf-schema#>
PREFIX skos: <http://www.w3.org/2004/02/skos/core#>
PREFIX owl: <http://www.w3.org/2002/07/owl#>
"""


@lru_cache(maxsize=1)
def ontology() -> Graph:
    return Graph().parse(ROOT / "ontology" / "oncokg.ttl")


def q(sparql):
    return list(ontology().query(PREFIX + sparql))


def local(u):
    return str(u).split("#")[-1]


@lru_cache(maxsize=1)
def params():
    return {str(n): float(v) for n, v in q("SELECT ?n ?v WHERE { ?p a okg:Parameter ; rdfs:label ?n ; okg:paramValue ?v }")}


def label(u):
    return str(ontology().value(u, RDFS.label) or local(u))
