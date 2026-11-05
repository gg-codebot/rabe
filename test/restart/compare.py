"""Require output and restart file of two run directories to be byte-identical."""

import filecmp
import sys
from pathlib import Path

full, resumed = Path(sys.argv[1]), Path(sys.argv[2])
failed = False
for name in ("electric_rabe.nc", "saved_bounce_integrals.dat"):
    if filecmp.cmp(full / name, resumed / name, shallow=False):
        print(f"PASS: {name} identical")
    else:
        print(f"FAIL: {name} differs between {full} and {resumed}")
        failed = True
sys.exit(1 if failed else 0)
