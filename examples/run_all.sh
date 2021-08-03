#!/bin/bash
set -euo pipefail

# Only look at documented examples, and run them in that order:
for example in `grep "/index" ../docs/examples/index.rst | cut -f 1 -d "/" | cut -c 4-`; do
    echo "==========================="
    echo "Running $example"
    echo "==========================="
    cd $example
    time ./run.sh
    cd ..
done

echo "==========================="
echo "Ran all documented examples"
echo "==========================="
