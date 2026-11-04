# restart

End-to-end test of the `electric_rabe` restart file `electric_rabe.restart.dat`
(label `slow`).

1. `RestartFullRun` runs `electric_rabe` from scratch in `full/`.
2. `RestartTruncate` copies its restart file to `resumed/`, cut off in the middle
   of the last fieldline as if the run had been killed (`truncate_restart.py`).
3. `RestartResumedRun` resumes from that file in `resumed/` and must report that
   it read one of the two fieldlines.
4. `RestartCompareOutput` and `RestartCompareRestartFile` require
   `electric_rabe.nc` and the restart file of both runs to be byte-identical.

The input uses the field file of `test/integration/vmec/input` with the smallest
possible number of fieldlines (two), since every fieldline costs about a minute
of bounce integrals.

Single test: `ctest --test-dir build/test -R Restart --output-on-failure`
