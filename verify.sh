#!/usr/bin/env bash
# Checks the C45 certificate. Needs a C++20 compiler and Python 3.10+ (standard library only).
#   ./verify.sh --stored        # check the stored histograms in certs/ (seconds)
#   ./verify.sh cert-51.json    # rebuild the 51-prime histogram and check it (~17 min, THREADS=2)
#   ./verify.sh                 # rebuild the 52-prime histogram and check it (~1.5 h, THREADS=2, ~3.6 GB)
set -euo pipefail
cd "$(dirname "$0")"
if [[ "${1:-}" == "--stored" ]]; then
    for c in certs/c45-52 certs/c45-51; do
        python3 verify_cert.py "$c/cert.json" "$c/seed_histogram.txt" --self-test
    done
    exit 0
fi
cert="${1:-cert.json}"
cxx="${CXX:-c++}"
"$cxx" -O3 -std=c++20 -pthread dp_seed2.cpp -o dp_seed2
./dp_seed2 --self-test
# shellcheck disable=SC2046
./dp_seed2 --threads "${THREADS:-2}" $(python3 -c 'import json, sys; print(" ".join(json.load(open(sys.argv[1]))["dp_args"]))' "$cert") > seed_histogram.txt
python3 verify_cert.py "$cert" seed_histogram.txt --self-test --json-out certificate_summary.json
