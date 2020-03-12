#!/bin/bash
set -eup pipeline

for example in woody_hosts recycled_water fungal_mock microalgal_mock fecal_sequel endangered_species; do
    echo "========================="
    echo "Running $example"
    echo "========================="
    cd $example
    time ./run.sh
    cd ..
done

echo "========================="
echo "Ran all examples"
echo "========================="
