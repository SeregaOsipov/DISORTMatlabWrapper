This repository allows to call DISORT (Fortran program) as a subtouine from Matlab. This approach is advantegeous because it skips the text input/output from the extrernal files. The results of the radiative transfer (fluxes, radiances) are returned as Matlab variables.

If you use this code, cite Osipov et al., 2021 (https://doi.org/10.1038/s43247-021-00141-7).

# How to call DISORT as a subroutine from Matlab
The most important thing is to compile the Fortran code using mex corectly. The Fortran-Matlab wrapper will pass 8-byte integer and reals. Make sure that mex compiles with the -fdefault-real-8 and -fdefault-integer-8 flags. My version of Matlab supplies -fdefault-integer-8 by default and I add -fdefault-real-8 manually.

# Possible crash issues
1. Wrong data types (explained above)
2. Wrong BDREF file (there are two)
3. Finally, you might end up compiling mex file in the wrong dirrectory (current vs root) and it will ambiguous which mex file you are using

The test case (for debugging; ported DISOTEST.f) shows how to use mex and call MATLABDisortWrapper https://github.com/SeregaOsipov/DISORTMatlabWrapper/blob/master/portedUnitTestsDisortBeta2_0.m

If you run into problems, familirize yourself with the Matlab docs on MEX https://de.mathworks.com/help/matlab/call-mex-functions.html

DISORT source code in Fortran is here: http://www.rtatmocn.com/disort/
If you want coupling with the legacy DISORTv2.0beta, then a have a look at v2.0beta tag https://github.com/SeregaOsipov/DISORTMatlabWrapper/releases/tag/v2.0beta
Since that tag, code switched to disort4.0.99.
