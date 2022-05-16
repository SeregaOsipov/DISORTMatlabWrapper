# How to call DISORT as a subroutine from Matlab
Matlab help on MEX https://de.mathworks.com/help/matlab/call-mex-functions.html

Compile the wrapper and ALL other F files using MEX functions https://github.com/SeregaOsipov/DISORTMatlabWrapper/blob/master/DISORT2.0beta/DISORTMatlabWrapper.F  
go to the root dir and execute:  

mex /work/mm0062/b302074/workspace/fortran/DISORTMatlabWrapper/DISORT2.0beta/DISORTMatlabWrapper.F /work/mm0062/b302074/workspace/fortran/DISORTMatlabWrapper/DISORT2.0beta/DISORT.f /work/mm0062/b302074/workspace/fortran/DISORTMatlabWrapper/DISORT2.0beta/RDI1MACH.f /work/mm0062/b302074/workspace/fortran/DISORTMatlabWrapper/DISORT2.0beta/BDREFCopyFromTest.f /work/mm0062/b302074/workspace/fortran/DISORTMatlabWrapper/DISORT2.0beta/ErrPackMexCompatible.f /work/mm0062/b302074/workspace/fortran/DISORTMatlabWrapper/DISORT2.0beta/LINPAK.f

