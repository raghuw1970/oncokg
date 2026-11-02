#!/usr/bin/env bash
# End-to-end: ontology -> synthetic raw data -> normalization -> reasoning -> evaluation -> explorer
set -euo pipefail
cd "$(dirname "$0")"
python3 ontology/build_ontology.py
python3 data/generate_synthetic.py
cd pipeline
python3 normalize.py > /dev/null
python3 reason.py > /dev/null
python3 evaluate.py > /dev/null
python3 build_explorer.py
cp ../output/oncokg_explorer.html ../docs/index.html
echo "done -> output/ and docs/index.html"
