# How to call DISORT as a subroutine from Matlab
The most important thing is to compile the Fortran code using mex corectly. The Fortran-Matlab wrapper will pass 8-byte integer and reals. Make sure that mex compiles with the -fdefault-real-8 and -fdefault-integer-8 flags. My version of Matlab supplies -fdefault-integer-8 by default and I add -fdefault-real-8 manually.

# Possible crash issues
1. Wrong data types (explained above)
2. Wrong BDREF file (there are two)
3. Finally, you might end up compiling mex file in the wrong dirrectory (current vs root) and it will ambiguous which mex file you are using

Here is the test case for debugging (ported DISOTEST.f), which shows how to use mex and call MATLABDisortWrapper https://github.com/SeregaOsipov/DISORTMatlabWrapper/blob/master/portedUnitTestsDisortBeta2_0.m

If you run into problems, familirize yourself with the Matlab docs on MEX https://de.mathworks.com/help/matlab/call-mex-functions.html
