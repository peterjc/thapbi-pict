#!/bin/bash
set -eup pipeline

if [ -f raw_data.tar.gz ]; then
    echo "Already have taxdmp_2019-12-01/"
else
    echo "Downloading raw data"
    wget https://zenodo.org/record/3342957/files/raw_data.tar.gz
    # TODO - verify MD5 of tar ball is 3435aa5b567897ab23e6aac46fb822a9
    echo "Decompressing raw data"
    tar -jxvkf raw_data.tar.gz
fi
