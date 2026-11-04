#!/bin/sh
# Usage: run_restart_test.sh <electric_rabe> <test/restart source dir> <field file>
set -e
electric_rabe=$1
source_dir=$2
field_file=$3

rm -rf full resumed
mkdir full
cp "$source_dir/input/rabe.in" "$field_file" full/

(cd full && "$electric_rabe")
cp -r full resumed
python3 "$source_dir/truncate.py" resumed
(cd resumed && "$electric_rabe" > run.log)
grep -Eq "Read +1 +of +2 +fieldlines" resumed/run.log || {
    echo "FAIL: resumed run did not reuse fieldline 1 from the restart file"
    exit 1
}
python3 "$source_dir/compare.py" full resumed
