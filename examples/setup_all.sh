#!/bin/bash
set -eup pipeline

for example in woody_hosts recycled_water fungal_mock microalgal_mock fecal_sequel endangered_species; do
    echo "========================="
    echo "Setting up $example"
    echo "========================="
    cd $example
    ./setup.sh
    cd ..
done

echo "========================="
echo "Setup all examples"
echo "========================="
