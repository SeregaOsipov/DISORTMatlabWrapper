%the following line will compile the fortran code addopted for matlab
clc
clear all
close all

format short e
ACCUR = single(0.0);
PRNT = ones(5,1);%true(5,1);
HEADER = 'header sample string';
%%
%this one compiles DISORTMatlab
mex /home/osipovs/workspace/Fortran/DISORT/DISORT2.0beta/DISORTMatlabWrapper.F /home/osipovs/workspace/Fortran/DISORT/DISORT2.0beta/DISORT.f /home/osipovs/workspace/Fortran/DISORT/DISORT2.0beta/RDI1MACH.f /home/osipovs/workspace/Fortran/DISORT/DISORT2.0beta/BDREFCopyFromTest.f /home/osipovs/workspace/Fortran/DISORT/DISORT2.0beta/ErrPackMexCompatible.f /home/osipovs/workspace/Fortran/DISORT/DISORT2.0beta/LINPAK.f

%%%%% below these are from the disotest.f

%% Isotropic scattering
NSTR = int32(16);
NLYR = int32(1);
NMOM = NSTR;
PMOM = single(getmom( 1, 0.0, NMOM));
% PMOM = single(zeros(NMOM+1, NLYR));
% PMOM(1) = 1;
USRTAU    = true;
NTAU      = int32(2);
UTAU = single(zeros(NTAU,1));
UTAU( 1 ) = 0.0;
USRANG    = true;
% USRANG    = false;
NUMU      = int32(6);
UMU = single(zeros(NUMU,1));
UMU( 1 )  = -1.0;
UMU( 2 )  = -0.5;
UMU( 3 )  = -0.1;
UMU( 4 )  =  0.1;
UMU( 5 )  =  0.5;
UMU( 6 )  =  1.0;

    NUMU      = int32(16);
    UMU = single(linspace(-1,1, NUMU));    
    if ( any(UMU==0))
        ME = MException('DISORT:UMU', 'has 0 values!');
        throw(ME) 
    end

NPHI      = int32(1);
PHI = single(zeros(NPHI));
IBCND     = int32(0);
UMU0      = single(0.1);
PHI0      = single(0.0);
LAMBER    = true;
ALBEDO    = single(0.0);
PLANK     = false;
ONLYFL    = false;

UTAU( 2 )  = 0.03125;
SSALB = single(0*ones(NLYR,1));
SSALB( 1 ) = 0.2;
FBEAM      = single(pi / UMU0);
FISOT      = single(0.0);

DTAUC = single(zeros(NLYR,1));
DTAUC( 1 ) = UTAU( 2 );

TEMPER = single(zeros(NLYR+1,1));
WVNMLO = single(0);
WVNMHI = single(0);
BTEMP = single(0);
TTEMP = single(0);
TEMIS = single(0);
% PLANK = 1;%true;
% ONLYFL = 1;%true;
ACCUR = single(0.0);
PRNT = ones(5,1);%true(5,1);
HEADER = 'header sample string';

% MAXCLY = int32(6);
% MAXMOM = int32(299);
% MAXPHI = int32(3);
% MAXULV = int32(5);
% MAXUMU = int32(10);

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
MAXUMU = int32(256);
MAXPHI = NPHI;
%% run the disort
for ICAS = 1:1%6
	if ( ICAS == 1 ) 
        UTAU( 2 )  = 0.03125
        SSALB( 1 ) = 0.2
        FBEAM      = single(pi / UMU0);
        FISOT      = single(0.0);
    elseif ( ICAS == 2 )
        UTAU( 2 )  = 0.03125
        SSALB( 1 ) = 1.0
        FBEAM      = single(pi / UMU0);
        FISOT      = single(0.0);
    elseif ( ICAS == 3 )
        UTAU( 2 )  = 0.03125
        SSALB( 1 ) = 0.99
        FBEAM      = single(0.0);
        FISOT      = single(1.0);
    elseif( ICAS == 4 )
        UTAU( 2 )  = 32.0
        SSALB( 1 ) = 0.2
        FBEAM      = single(pi / UMU0);
        FISOT      = single(0.0);
    elseif( ICAS == 5 )
        UTAU( 2 )  = 32.0
        SSALB( 1 ) = 1.0
        FBEAM      = single(pi / UMU0);
        FISOT      = single(0.0);
    elseif( ICAS == 6 )
        UTAU( 2 )  = 32.0
        SSALB( 1 ) = 0.99
        FBEAM      = single(0.0)'
        FISOT      = single(1.0);
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
NSTR = int32(16);
NLYR = int32(1);
NMOM = NSTR;
PMOM = single(getmom( 2, 0.0, NMOM));
USRTAU    = true;
NTAU      = int32(2);
UTAU( 1 ) = 0.0;
USRANG    = true;
NUMU      = int32(6);
UMU = single(zeros(NUMU,1));
UMU( 1 )  = -0.981986;
UMU( 2 )  = -0.538263;
UMU( 3 )  = -0.018014;
UMU( 4 )  =  0.018014;
UMU( 5 )  =  0.538263;
UMU( 6 )  =  0.981986;
NPHI      = int32(1);
PHI(1) = 0;
IBCND     = int32(0);
FBEAM      = single(pi);
UMU0      = single(0.080442);
PHI0      = single(0.0);
FISOT      = single(0.0);
LAMBER    = true;
ALBEDO    = single(0.0);
PLANK     = false;
ONLYFL    = false;

% SSALB = single(zeros(NLYR,1));
% DTAUC = single(zeros(NLYR,1));

% TEMPER = single(zeros(NLYR+1,1));
% WVNMLO = single(0);
% WVNMHI = single(0);
% BTEMP = single(0);
% TTEMP = single(0);
% TEMIS = single(0);
% ACCUR = single(0.0);
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

%% Test Problem 7

UMU0 = single(0.5);
PHI0 = single(0.0);

TEMPER = single(zeros(NLYR+1,1));
TEMPER( 1 ) = 300.0;
TEMPER( 2 ) = 200.0;
WVNMLO      = single(0.0);
WVNMHI      = single(50000);
BTEMP       = single(320.0);
TTEMP       = single(100.0);
TEMIS       = single(1.0);



%% Test Problem 8:  Absorbing/Isotropic-Scattering Medium      With Two Computational Layers
clc
NSTR = int32(8);
NLYR = int32(2);
NMOM = NSTR;
PMOM = single(zeros(NMOM+1, NLYR));
PMOM(:,1) = single(getmom(1, 0.0, NMOM));
PMOM(:,2) = single(getmom(1, 0.0, NMOM));
% CALL  GETMOM( 1, 0.0, NMOM, PMOM(0,1) )
% CALL  GETMOM( 1, 0.0, NMOM, PMOM(0,2) )
USRTAU   = true;
USRANG   = true;
NUMU     = int32(4)
UMU( 1 ) = -1.0;
UMU( 2 ) = -0.2;
UMU( 3 ) =  0.2;
UMU( 4 ) =  1.0;
NPHI     = int32(1);
PHI( 1 ) = 60.0;
IBCND    = int32(0);
FBEAM    = single(0.0);
FISOT    = single(1.0 / pi);
LAMBER   = true;
ALBEDO   = single(0.0);
PLANK    = false;
ONLYFL   = false;

for ICAS = 1:3
    if ( ICAS == 1 )
        DTAUC( 1 ) = 0.25
        DTAUC( 2 ) = 0.25
        SSALB( 1 ) = 0.5
        SSALB( 2 ) = 0.3
        NTAU       = int32(3);
        UTAU( 1 )  = 0.0
        UTAU( 2 )  = 0.25
        UTAU( 3 )  = 0.5
        HEADER = 'Test Case No. 8a:  Ref. OS, Table 1, Line 4 (Two Inhomogeneous Layers)';
    elseif ( ICAS == 2 ) 
        DTAUC( 1 ) = 0.25
        DTAUC( 2 ) = 0.25
        SSALB( 1 ) = 0.8
        SSALB( 2 ) = 0.95
        NTAU       = int32(3);
        UTAU( 1 )  = 0.0
        UTAU( 2 )  = 0.25
        UTAU( 3 )  = 0.5
        HEADER = 'Test Case No. 8b:  Ref. OS, Table 1, Line 1 (Two Inhomogeneous Layers)';
    elseif ( ICAS == 3 ) 
        DTAUC( 1 ) = 1.0
        DTAUC( 2 ) = 2.0
        SSALB( 1 ) = 0.8
        SSALB( 2 ) = 0.95
        NTAU       = int32(3);
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
    DTAUC = single(DTAUC);
    SSALB = single(SSALB);
    UTAU = single(UTAU);
    UMU = single(UMU);
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






                    
%% these are the required parameters for the DISORT
NLYR = int32(10);
NMOM = int32(11);
MAXCLY = int32(3);
MAXULV = int32(12);
MAXUMU = int32(13);
MAXPHI = int32(2);
MAXMOM = int32(2);

DTAUC = single(0.21*ones(NLYR,1));
SSALB = single(0.22*ones(NLYR,1));

PMOM = single(0.23*ones(NMOM+1, NLYR));
TEMPER = single(0.24*ones(NLYR+1,1));
WVNMLO = single(12.5);
WVNMHI = single(12.6);
USRTAU = 1;%true;
NSTR = int32(16);
USRANG = 1;%true;

NPHI = int32(5);
PHI = single(0.275*ones(NPHI));
IBCND = int32(2);
FBEAM = single(0.5);
UMU0 = single(0.4);
PHI0 = single(0.3);
FISOT = single(0.2);
LAMBER = 1;%true;
ALBEDO = single(0.1);
BTEMP = single(0.6);
TTEMP = single(0.7);
TEMIS = single(0.8);
PLANK = 1;%true;
ONLYFL = 1;%true;
ACCUR = single(0.1);
PRNT = ones(5,1);%true(5,1);
HEADER = 'header sample string';

%these are the output variables
NTAU = int32(10);
UTAU = 0.25*ones(NTAU,1);

NUMU = int32(10);
UMU = (0.26*ones(NUMU,1));

%%
%this one compile DISOTESTMatlab
mex /home/osipovs/workspace/Fortran/DISORT/DISORT2.0beta/DISOTESTMatlab.F /home/osipovs/workspace/Fortran/DISORT/DISORT2.0beta/DISORT.f /home/osipovs/workspace/Fortran/DISORT/DISORT2.0beta/RDI1MACH.f /home/osipovs/workspace/Fortran/DISORT/DISORT2.0beta/ErrPackMexCompatible.f /home/osipovs/workspace/Fortran/DISORT/DISORT2.0beta/LINPAK.f
%%
try
    asd = DISOTESTMatlab()
catch error
    error.message
end
