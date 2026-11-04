"""Make a run directory look like electric_rabe was killed mid-fieldline.

Cuts the restart file in the middle of a line of the last fieldline and
removes the output file.
"""

import sys
from pathlib import Path

run_dir = Path(sys.argv[1])
restart = run_dir / "saved_bounce_integrals.dat"

lines = restart.read_text().splitlines(keepends=True)
last_id = lines[-1].split()[0]
first_of_last = next(i for i, line in enumerate(lines) if line.split()[0] == last_id)
cut = first_of_last + (len(lines) - first_of_last) // 2
restart.write_text("".join(lines[:cut]) + lines[cut][: len(lines[cut]) // 2])

(run_dir / "electric_rabe.nc").unlink()
