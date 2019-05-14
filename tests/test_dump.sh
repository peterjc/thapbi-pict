#!/bin/bash
IFS=$'\n\t'
set -eux
# Note not using "set -o pipefail" until after check error message with grep

# Note all tests here (initially) using default database:

echo "Checking dump"
thapbi_pict dump -s fallax 2>&1 | grep "species requires a single genus"
set -o pipefail

if [ `thapbi_pict dump | grep -c -v "^#"` -ne 14988 ]; then echo "Wrong record count for table dump"; false; fi
if [ `thapbi_pict dump -f fasta | grep -c "^>"` -ne 14988 ]; then echo "Wrong record count for fasta dump"; false; fi

# With genus filter,
if [ `thapbi_pict dump -f fasta -g Phytophthora | grep -c "^>"` -ne 13020 ]; then echo "Wrong record count for Phytophthora genus"; false; fi

# With genus and species filter,
if [ `thapbi_pict dump -f fasta -g Phytophthora -s "fallax, andina" | grep -c "^>"` -ne 7 ]; then echo "Wrong record count for two species"; false; fi

# With genus and clades
if [ `thapbi_pict dump -f fasta -g Phytophthora -c 8a,8b | grep -c "^>"` -ne 17 ]; then echo "Wrong record count for Phytophthora clade"; false; fi

# With clades only (by DB construction will be Phytophthora only)
if [ `thapbi_pict dump -c 8a,8b | grep -c -v "^#"` -ne 17 ]; then echo "Wrong record count for clade"; false; fi

echo "$0 - test_dump.sh passed"
