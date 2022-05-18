%the following line will compile the fortran code adopted for matlab
clc
clear all
close all

format short e
ACCUR = double(0.0);
PRNT = ones(5,1);%true(5,1);
HEADER = 'header sample string';

%%
% this one compiles DISORTMatlab. Double precision (default)
% mex /work/mm0062/b302074/workspace/fortran/DISORTMatlabWrapper/DISORT2.0beta/DISORTMatlabWrapper.F /work/mm0062/b302074/workspace/fortran/DISORTMatlabWrapper/DISORT2.0beta/DISORT.f /work/mm0062/b302074/workspace/fortran/DISORTMatlabWrapper/DISORT2.0beta/RDI1MACH.f /work/mm0062/b302074/workspace/fortran/DISORTMatlabWrapper/DISORT2.0beta/BDREFCopyFromTest.f /work/mm0062/b302074/workspace/fortran/DISORTMatlabWrapper/DISORT2.0beta/ErrPackMexCompatible.f /work/mm0062/b302074/workspace/fortran/DISORTMatlabWrapper/DISORT2.0beta/LINPAK.f
mex FFLAGS='$FFLAGS -fdefault-real-8'  /work/mm0062/b302074/workspace/fortran/DISORTMatlabWrapper/DISORT2.0beta/DISORTMatlabWrapper.F /work/mm0062/b302074/workspace/fortran/DISORTMatlabWrapper/DISORT2.0beta/DISORT.f /work/mm0062/b302074/workspace/fortran/DISORTMatlabWrapper/DISORT2.0beta/RDI1MACH.f /work/mm0062/b302074/workspace/fortran/DISORTMatlabWrapper/DISORT2.0beta/BDREFCopyFromTest.f /work/mm0062/b302074/workspace/fortran/DISORTMatlabWrapper/DISORT2.0beta/ErrPackMexCompatible.f /work/mm0062/b302074/workspace/fortran/DISORTMatlabWrapper/DISORT2.0beta/LINPAK.f
% Make sure to choose correct BDREF.f (WILL CASE FATAL ERROR) or BDREFCopyFromTest.f
% on Levante/2022, if loading gcc module, then add gfortran flags: mex FFLAGS='$FFLAGS -fallow-argument-mismatch -fallow-invalid-boz' ...
% but in cant load correct libgfortran
% add mex -v fo verbose

%%%%% below these are from the disotest.f
%% Isotropic scattering
NSTR = int64(16);
NLYR = int64(1);
NMOM = NSTR;
PMOM = double(getmom( 1, 0.0, NMOM));
% PMOM = double(zeros(NMOM+1, NLYR));
% PMOM(1) = 1;
USRTAU    = true;
NTAU      = int64(2);
UTAU = double(zeros(NTAU,1));
UTAU( 1 ) = 0.0;
USRANG    = true;
% USRANG    = false;
NUMU      = int64(6);
UMU = double(zeros(NUMU,1));
UMU( 1 )  = -1.0;
UMU( 2 )  = -0.5;
UMU( 3 )  = -0.1;
UMU( 4 )  =  0.1;
UMU( 5 )  =  0.5;
UMU( 6 )  =  1.0;

%     NUMU      = int64(8);
%     UMU = double(linspace(-1,1, NUMU));    
    if ( any(UMU==0))
        ME = MException('DISORT:UMU', 'has 0 values!');
        throw(ME) 
    end

NPHI      = int64(1);
PHI = double(zeros(NPHI));
IBCND     = int64(0);
UMU0      = double(0.1);
PHI0      = double(0.0);
LAMBER    = true;
ALBEDO    = double(0.0);
PLANK     = false;
ONLYFL    = false;

UTAU( 2 )  = 0.03125;
SSALB = double(0*ones(NLYR,1));
SSALB( 1 ) = 0.2;
FBEAM      = double(pi / UMU0);
FISOT      = double(0.0);

DTAUC = double(zeros(NLYR,1));
DTAUC( 1 ) = UTAU( 2 );

TEMPER = double(zeros(NLYR+1,1));
WVNMLO = double(0);
WVNMHI = double(0);
BTEMP = double(0);
TTEMP = double(0);
TEMIS = double(0);
% PLANK = 1;%true;
% ONLYFL = 1;%true;
ACCUR = double(0.0);
PRNT = ones(5,1);%true(5,1);
HEADER = 'header sample string';

% MAXCLY = int64(6);
% MAXMOM = int64(299);
% MAXPHI = int64(3);
% MAXULV = int64(5);
% MAXUMU = int64(10);

MAXCLY = NLYR;
MAXMOM = NMOM;
MAXULV = NLYR+1;
if ( USRTAU )
    MAXULV = NTAU;
end
MAXUMU = NSTR;
if ( USRANG )
    MAXUMU = NUMU;    
end
MAXUMU = int64(256);
MAXPHI = NPHI;
%% run the disort
for ICAS = 6:6%6
	if ( ICAS == 1 ) 
        UTAU( 2 )  = 0.03125
        SSALB( 1 ) = 0.2
        FBEAM      = double(pi / UMU0);
        FISOT      = double(0.0);
    elseif ( ICAS == 2 )
        UTAU( 2 )  = 0.03125
        SSALB( 1 ) = 1.0
        FBEAM      = double(pi / UMU0);
        FISOT      = double(0.0);
    elseif ( ICAS == 3 )
        UTAU( 2 )  = 0.03125
        SSALB( 1 ) = 0.99
        FBEAM      = double(0.0);
        FISOT      = double(1.0);
    elseif( ICAS == 4 )
        UTAU( 2 )  = 32.0
        SSALB( 1 ) = 0.2
        FBEAM      = double(pi / UMU0);
        FISOT      = double(0.0);
    elseif( ICAS == 5 )
        UTAU( 2 )  = 32.0
        SSALB( 1 ) = 1.0
        FBEAM      = double(pi / UMU0);
        FISOT      = double(0.0);
    elseif( ICAS == 6 )
        UTAU( 2 )  = 32.0
        SSALB( 1 ) = 0.99
        FBEAM      = double(0.0);
        FISOT      = double(1.0);
    end

    DTAUC( 1 ) = UTAU( 2 )

    display(['Test Case No. 1_' num2str(ICAS) ':  Isotropic Scattering, Ref. VH1, Table 12:  b =' num2str(UTAU(2)) ', a =' num2str(SSALB(1))])
    tic
    [RFLDIR, RFLDN, FLUP, DFDT, UAVG, UU, ALBMED, TRNMED] = DISORTMatlabWrapper(NLYR, DTAUC, SSALB, NMOM, PMOM, TEMPER, WVNMLO, ... %1-7
                        WVNMHI, USRTAU, NTAU, UTAU, NSTR, USRANG, NUMU, ...%8-14
                        UMU, NPHI, PHI, IBCND, FBEAM, UMU0, PHI0, ...%15-21
                        FISOT, LAMBER, ALBEDO, BTEMP, TTEMP, TEMIS, ...%22-27
                        PLANK, ONLYFL, ACCUR, PRNT, HEADER, MAXCLY, ...%28-33
                        MAXULV, MAXUMU, MAXPHI, MAXMOM)%34-37
%internally DISORT stores the UU variable in MAX possible size of the dimensions
%so we have to cut it done here manually
% UU = UU(1:NUMU, 1:NTAU, 1:NPHI)
    toc
end


%% Rayleigh Scattering, Beam Source
clc
NSTR = int64(16);
NLYR = int64(1);
NMOM = NSTR;
PMOM = double(getmom( 2, 0.0, NMOM));
USRTAU    = true;
NTAU      = int64(2);
UTAU( 1 ) = 0.0;
USRANG    = true;
NUMU      = int64(6);
UMU = double(zeros(NUMU,1));
UMU( 1 )  = -0.981986;
UMU( 2 )  = -0.538263;
UMU( 3 )  = -0.018014;
UMU( 4 )  =  0.018014;
UMU( 5 )  =  0.538263;
UMU( 6 )  =  0.981986;
NPHI      = int64(1);
PHI(1) = 0;
IBCND     = int64(0);
FBEAM      = double(pi);
UMU0      = double(0.080442);
PHI0      = double(0.0);
FISOT      = double(0.0);
LAMBER    = true;
ALBEDO    = double(0.0);
PLANK     = false;
ONLYFL    = false;

% SSALB = double(zeros(NLYR,1));
% DTAUC = double(zeros(NLYR,1));

% TEMPER = double(zeros(NLYR+1,1));
% WVNMLO = double(0);
% WVNMHI = double(0);
% BTEMP = double(0);
% TTEMP = double(0);
% TEMIS = double(0);
% ACCUR = double(0.0);
% PRNT = ones(5,1);%true(5,1);
% HEADER = 'header sample string';


% MAXCLY = NLYR;
% MAXMOM = NMOM;
% MAXULV = NLYR+1;
% if ( USRTAU )
%     MAXULV = NTAU;
% end
% MAXUMU = NSTR;
% if ( USRANG )
%     MAXUMU = NUMU;
% end
% MAXPHI = NPHI;

ICAS = 0;
for IOD = 1:2
    if ( IOD == 1 )
        UTAU( 2 ) = 0.2
    end
    if ( IOD == 2 )
        UTAU( 2 ) = 5.0
    end

    DTAUC( 1 ) = UTAU( 2 )

    for ISS = 1:2
        if( ISS == 1 ) 
            SSALB( 1 ) = 0.5
        end
        if( ISS == 2 )
            SSALB( 1 ) = 1.0
        end
        ICAS = ICAS + 1

        display(['Test Case No. 2_' num2str(ICAS) ', Rayleigh Scattering, Ref. SW, Table 1:  tau =' num2str(UTAU(2)) ', mu0 =', num2str(UMU0) ', ss-albedo =' num2str(SSALB(1))])
        tic
        [RFLDIR, RFLDN, FLUP, DFDT, UAVG, UU, ALBMED, TRNMED] = DISORTMatlabWrapper(NLYR, DTAUC, SSALB, NMOM, PMOM, TEMPER, WVNMLO, ... %1-7
                            WVNMHI, USRTAU, NTAU, UTAU, NSTR, USRANG, NUMU, ...%8-14
                            UMU, NPHI, PHI, IBCND, FBEAM, UMU0, PHI0, ...%15-21
                            FISOT, LAMBER, ALBEDO, BTEMP, TTEMP, TEMIS, ...%22-27
                            PLANK, ONLYFL, ACCUR, PRNT, HEADER, MAXCLY, ...%28-33
                            MAXULV, MAXUMU, MAXPHI, MAXMOM)%34-37
    %internally DISORT stores the UU variable in MAX possible size of the dimensions
    %so we have to cut it done here manually
    %UU = UU(1:NUMU, 1:NTAU, 1:NPHI)
        toc
    end
end

%% Test Problem 3:  Henyey-Greenstein Scattering Compare To Ref. VH2, Table 37
clc
NSTR = int64(16);
NLYR = int64(1);
SSALB( 1 ) = 1.0;
NMOM = int64(32);
PMOM = double(getmom( 3, 0.75, NMOM));
USRTAU    = true;
NTAU      = int64(2);
UTAU( 1 ) = 0.0;
USRANG    = true;
NUMU      = int64(6);
UMU = double(zeros(NUMU,1));
UMU( 1 )  = -1.0;
UMU( 2 )  = -0.5;
UMU( 3 )  = -0.1;
UMU( 4 )  =  0.1;
UMU( 5 )  =  0.5;
UMU( 6 )  =  1.0;
NPHI      = int64(1);
PHI = double(zeros(NPHI));
PHI( 1 )  = 0.0;
IBCND     = int64(0);
UMU0      = double(1.0);
PHI0      = double(0.0);
FBEAM      = double(pi / UMU0);
FISOT      = double(0.0);
LAMBER    = true;
ONLYFL    = false;
ALBEDO    = double(0.0);
PLANK     = false;

for ICAS = 1:2
    if ( ICAS == 1 )
        UTAU( 2 )  = 1.0;
    elseif ( ICAS == 2 )
        UTAU( 2 )  = 8.0;
    end

    DTAUC( 1 ) = UTAU( 2 );
    display( ['Test Case No. 3_' num2str(ICAS) ', Henyey-Greenstein Scattering, Ref. VH2, Table 37, g = 0.75, b =' num2str(UTAU(2)) ', a ='  num2str(SSALB(1))])

    
    %extra from me
    MAXMOM = NMOM;
    [RFLDIR, RFLDN, FLUP, DFDT, UAVG, UU, ALBMED, TRNMED] = DISORTMatlabWrapper(NLYR, DTAUC, SSALB, NMOM, PMOM, TEMPER, WVNMLO, ... %1-7
                    WVNMHI, USRTAU, NTAU, UTAU, NSTR, USRANG, NUMU, ...%8-14
                    UMU, NPHI, PHI, IBCND, FBEAM, UMU0, PHI0, ...%15-21
                    FISOT, LAMBER, ALBEDO, BTEMP, TTEMP, TEMIS, ...%22-27
                    PLANK, ONLYFL, ACCUR, PRNT, HEADER, MAXCLY, ...%28-33
                    MAXULV, MAXUMU, MAXPHI, MAXMOM)%34-37

end

%% Test Problem 7

UMU0 = double(0.5);
PHI0 = double(0.0);

TEMPER = double(zeros(NLYR+1,1));
TEMPER( 1 ) = 300.0;
TEMPER( 2 ) = 200.0;
WVNMLO      = double(0.0);
WVNMHI      = double(50000);
BTEMP       = double(320.0);
TTEMP       = double(100.0);
TEMIS       = double(1.0);



%% Test Problem 8:  Absorbing/Isotropic-Scattering Medium      With Two Computational Layers
clc
NSTR = int64(8);
NLYR = int64(2);
NMOM = NSTR;
PMOM = double(zeros(NMOM+1, NLYR));
PMOM(:,1) = double(getmom(1, 0.0, NMOM));
PMOM(:,2) = double(getmom(1, 0.0, NMOM));
% CALL  GETMOM( 1, 0.0, NMOM, PMOM(0,1) )
% CALL  GETMOM( 1, 0.0, NMOM, PMOM(0,2) )
USRTAU   = true;
USRANG   = true;
NUMU     = int64(4);
UMU( 1 ) = -1.0;
UMU( 2 ) = -0.2;
UMU( 3 ) =  0.2;
UMU( 4 ) =  1.0;
NPHI     = int64(1);
PHI( 1 ) = 60.0;
IBCND    = int64(0);
FBEAM    = double(0.0);
FISOT    = double(1.0 / pi);
LAMBER   = true;
ALBEDO   = double(0.0);
PLANK    = false;
ONLYFL   = false;

for ICAS = 1:3
    if ( ICAS == 1 )
        DTAUC( 1 ) = 0.25
        DTAUC( 2 ) = 0.25
        SSALB( 1 ) = 0.5
        SSALB( 2 ) = 0.3
        NTAU       = int64(3);
        UTAU( 1 )  = 0.0
        UTAU( 2 )  = 0.25
        UTAU( 3 )  = 0.5
        HEADER = 'Test Case No. 8a:  Ref. OS, Table 1, Line 4 (Two Inhomogeneous Layers)';
    elseif ( ICAS == 2 ) 
        DTAUC( 1 ) = 0.25
        DTAUC( 2 ) = 0.25
        SSALB( 1 ) = 0.8
        SSALB( 2 ) = 0.95
        NTAU       = int64(3);
        UTAU( 1 )  = 0.0
        UTAU( 2 )  = 0.25
        UTAU( 3 )  = 0.5
        HEADER = 'Test Case No. 8b:  Ref. OS, Table 1, Line 1 (Two Inhomogeneous Layers)';
    elseif ( ICAS == 3 ) 
        DTAUC( 1 ) = 1.0
        DTAUC( 2 ) = 2.0
        SSALB( 1 ) = 0.8
        SSALB( 2 ) = 0.95
        NTAU       = int64(3);
        UTAU( 1 )  = 0.0
        UTAU( 2 )  = 1.0
        UTAU( 3 )  = 3.0
        HEADER = 'Test Case No. 8c:  Ref. OS, Table 1, Line 13 (Two Inhomogeneous Layers)'
    end

    MAXCLY = NLYR;
    MAXMOM = NMOM;
    MAXULV = NLYR+1;
    if ( USRTAU )
        MAXULV = NTAU;
    end
    MAXUMU = NSTR;
    if ( USRANG )
        MAXUMU = NUMU;
    end
    MAXPHI = NPHI;

    display(HEADER)
    tic
    DTAUC = double(DTAUC);
    SSALB = double(SSALB);
    UTAU = double(UTAU);
    UMU = double(UMU);
    [RFLDIR, RFLDN, FLUP, DFDT, UAVG, UU, ALBMED, TRNMED] = DISORTMatlabWrapper(NLYR, DTAUC, SSALB, NMOM, PMOM, TEMPER, WVNMLO, ... %1-7
                        WVNMHI, USRTAU, NTAU, UTAU, NSTR, USRANG, NUMU, ...%8-14
                        UMU, NPHI, PHI, IBCND, FBEAM, UMU0, PHI0, ...%15-21
                        FISOT, LAMBER, ALBEDO, BTEMP, TTEMP, TEMIS, ...%22-27
                        PLANK, ONLYFL, ACCUR, PRNT, HEADER, MAXCLY, ...%28-33
                        MAXULV, MAXUMU, MAXPHI, MAXMOM)%34-37
%internally DISORT stores the UU variable in MAX possible size of the dimensions
%so we have to cut it done here manually
%UU = UU(1:NUMU, 1:NTAU, 1:NPHI)
    toc

end


return



                    
%% these are the required parameters for the DISORT
NLYR = int64(10);
NMOM = int64(11);
MAXCLY = int64(3);
MAXULV = int64(12);
MAXUMU = int64(13);
MAXPHI = int64(2);
MAXMOM = int64(2);

DTAUC = double(0.21*ones(NLYR,1));
SSALB = double(0.22*ones(NLYR,1));

PMOM = double(0.23*ones(NMOM+1, NLYR));
TEMPER = double(0.24*ones(NLYR+1,1));
WVNMLO = double(12.5);
WVNMHI = double(12.6);
USRTAU = 1;%true;
NSTR = int64(16);
USRANG = 1;%true;

NPHI = int64(5);
PHI = double(0.275*ones(NPHI));
IBCND = int64(2);
FBEAM = double(0.5);
UMU0 = double(0.4);
PHI0 = double(0.3);
FISOT = double(0.2);
LAMBER = 1;%true;
ALBEDO = double(0.1);
BTEMP = double(0.6);
TTEMP = double(0.7);
TEMIS = double(0.8);
PLANK = 1;%true;
ONLYFL = 1;%true;
ACCUR = double(0.1);
PRNT = ones(5,1);%true(5,1);
HEADER = 'header sample string';

%these are the output variables
NTAU = int64(10);
UTAU = 0.25*ones(NTAU,1);

NUMU = int64(10);
UMU = (0.26*ones(NUMU,1));


return

%%
%this one compile DISOTESTMatlab
%mex /work/mm0062/b302074/workspace/fortran/DISORTMatlabWrapper/DISORT2.0beta/DISOTESTMatlab.F /work/mm0062/b302074/workspace/fortran/DISORTMatlabWrapper/DISORT2.0beta/DISORT.f /work/mm0062/b302074/workspace/fortran/DISORTMatlabWrapper/DISORT2.0beta/RDI1MACH.f /work/mm0062/b302074/workspace/fortran/DISORTMatlabWrapper/DISORT2.0beta/ErrPackMexCompatible.f /work/mm0062/b302074/workspace/fortran/DISORTMatlabWrapper/DISORT2.0beta/LINPAK.f
% On Levante with module load gcc, add flags:
% FFLAGS='$FFLAGS -fallow-argument-mismatch -fallow-invalid-boz'
%%
try
    asd = DISOTESTMatlab()
catch error
    error.message
end
