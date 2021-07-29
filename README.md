# Pawtrack_Test_Particle
Test-particle simulation programs for "pawtrack-like" pitch angle distributions.

Core Programs:
main_Fgeopack.m: the main and start function calculating test-particles' unperturbed orbits
Integ_deltaW_parallel.m: the function calculating the work done by ULF electric fields, using unperturbed orbits calculated by main_Fgeopack.m

Auxiliary Programs:
MAGNETIC_FIELD_FORTRAN.m: the function calling mex-format GEOPACK package to model the magnetic field (TS05+IGRF)
plot_trajectory.m: the function plotting test-particles' unperturbed orbits
read_mms_mec.m: the function reading MEC MMS orbits data
GEOPACK_TRACE_smallstep.m: the function tracing magnetic field lines with higher accuracy, based on the GEOPACK function GEOPACK_TRACE.m

Folders:
TS05: the files needed for TS05 model inputs (solar wind conditions)
geopack_mex: mex-format GEOPACK package, in order to decrease programs' running time
