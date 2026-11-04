# restart

End-to-end test of resuming `electric_rabe` from its saved bounce integrals
`saved_bounce_integrals.dat` (test `RestartResume`, label `slow`).
`run_restart_test.sh` runs in `build/test/restart/`:

1. run `electric_rabe` from scratch in `full/`
2. copy `full/` to `resumed/`
3. `truncate.py` cuts the restart file in `resumed/` in the middle of the last
   fieldline, as if the run had been killed, and removes the output
4. run `electric_rabe` again in `resumed/`; it must report that it read at
   least one fieldline from the restart file
5. `compare.py` requires `electric_rabe.nc` and the restart file of both runs to
   be byte-identical

The input uses the field file of `test/integration/vmec/input` with the smallest
possible number of fieldlines (two), since every fieldline costs about a minute
of bounce integrals.

Single test: `ctest --test-dir build/test -R RestartResume --output-on-failure`
