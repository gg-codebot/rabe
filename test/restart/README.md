# restart

End-to-end test that an interrupted `electric_rabe` run resumes from its saved
bounce integrals and reproduces an uninterrupted run exactly.

`run_restart_test.sh` runs the steps of the test.

Field: `wout_LandremanPaul2021_QH_reactorScale_lowres_reference.nc` from
`test/integration/vmec/input`.

Fails in a dirty build, which does not reuse saved bounce integrals.
