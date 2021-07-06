#!/bin/bash
set -euo pipefail

for example in woody_hosts recycled_water fungal_mock great_lakes fecal_sequel soil_nematodes endangered_species; do
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
