#!/bin/bash
VERSION=`thapbi_pict -v | sed "s/THAPBI PICT //g"`
set -euo pipefail

./build_CURATED.sh && ./build_CURATED+NCBI.sh && ./build_ITS1_DB.sh

echo "All done, used THAPBI PICT $VERSION"
