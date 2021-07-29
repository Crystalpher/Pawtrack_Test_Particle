function [XF,YF,ZF,XX,YY,ZZ,L] = GEOPACK_TRACE_smallstep (XI,YI,ZI,DIR,RLIM,R0,IOPT,PARMOD,EXNAME,INNAME)
%GEOPACK_TRACE  TRACES A FIELD LINE FROM AN ARBITRARY POINT OF SPACE TO THE EARTH'S
%  SURFACE OR TO A MODEL LIMITING BOUNDARY.
% 
%  THE HIGHEST ORDER OF SPHERICAL HARMONICS IN THE MAIN FIELD EXPANSION USED
%  IN THE MAPPING IS CALCULATED AUTOMATICALLY. IF INNAME=IGRF_GSM, THEN AN IGRF MODEL
%  FIELD WILL BE USED, AND IF INNAME=DIP, A PURE DIPOLE FIELD WILL BE USED.
% 
%  IN ANY CASE, BEFORE CALLING TRACE, ONE SHOULD INVOKE RECALC, TO CALCULATE CORRECT
%  VALUES OF THE IGRF COEFFICIENTS AND ALL QUANTITIES NEEDED FOR TRANSFORMATIONS
%  BETWEEN COORDINATE SYSTEMS INVOLVED IN THIS CALCULATIONS.
% 
%  ALTERNATIVELY, THE SUBROUTINE RECALC CAN BE INVOKED WITH THE DESIRED VALUES OF
%  IYEAR AND IDAY (TO SPECIFY THE DIPOLE MOMENT), WHILE THE VALUES OF THE DIPOLE
%  TILT ANGLE PSI (IN RADIANS) AND ITS SINE (SPS) AND COSINE (CPS) CAN BE EXPLICITLY
%  SPECIFIED AND FORWARDED TO THE COMMON BLOCK GEOPACK1 (11th, 12th, AND 16th ELEMENTS, RESP.)
%
%
%  INPUT PARAMETER
%  [XI],[YI],[ZI]  GSM COORDS OF INITIAL POINT (IN EARTH RADII, 1 RE = 6371.2 km),
%  [DIR]           SIGN OF THE TRACING DIRECTION: IF DIR=1.0 THEN WE MOVE ANTIPARALLEL TO THE
%                  FIELD VECTOR (E.G. FROM NORTHERN TO SOUTHERN CONJUGATE POINT),
%                  AND IF DIR=-1.0 THEN THE TRACING GOES IN THE OPPOSITE DIRECTION.
%  [R0]            RADIUS OF A SPHERE (IN RE) FOR WHICH THE FIELD LINE ENDPOINT 
%  [RLIM]          UPPER LIMIT OF THE GEOCENTRIC DISTANCE, WHERE THE TRACING IS TERMINATED.
%  [IOPT]          RANGE OF THE Kp INDEXES
%  [PARMOD]        10 ELEMENT ARRAY, SOLAR WIND PRESSURE, MAGNETIC FIELD
%                  COMP., DST INDEX
%  [EXNAME]        NAME OF THE EXTERNAL FIELD MODEL
%  [INNAME]        NAME OF THE INTERNAL FIELD MODEL
%
%  OUTPUT PARAMETER
%  [XF],[YF],[ZF]  GSM COORDS OF THE LAST CALCULATED POINT OF A FIELD LINE
%  [XX],[YY],[ZZ]  ARRAYS, CONTAINING COORDS OF FIELD LINE POINTS. HERE THEIR MAXIMAL LENGTH WAS
%                  ASSUMED EQUAL TO 999.
%  [L]             ACTUAL NUMBER OF THE CALCULATED FIELD LINE POINTS. IF L EXCEEDS 999, TRACING
%                  TERMINATES, AND A WARNING IS DISPLAYED.
%
%  WRITTEN AND ADDED TO THIS PACKAGE:  APRIL 1, 2003,
%  AUTHOR:   N. A. TSYGANENKO
%
%  TRANSLATED FROM ORIGINAL FORTRAN APRIL 11, 2003
%  BY PAUL O'BRIEN
%  
%  REVISED BY PATRICK DAUM AUGUST 2005


% COMMON GEOPACK COEFFICIENTS
global GEOPACK1;
ERR=0.0001; 
L=0;
DS=0.5*DIR;
% DS=0.05*DIR;
X=XI;
Y=YI;
Z=ZI;
GEOPACK1.DS3=DIR;
AL=0.;

[R1,R2,R3] = GEOPACK_RHAND (X,Y,Z,IOPT,PARMOD,EXNAME,INNAME);

% |AD|=0.01 and its sign follows the rule:
% (1) if DIR=1 (tracing antiparallel to B vector) then the sign of AD is the same as of Br
% (2) if DIR=-1 (tracing parallel to B vector) then the sign of AD is opposite to that of Br
%  AD is defined in order to initialize the value of RR (radial distance at previous step):
AD=0.01;
if (X*R1+Y*R2+Z*R3 < 0.), AD=-0.01; end;

RR=sqrt(X^2+Y^2+Z^2)+AD;
while 1,
    breakcode = 0;
    L=L+1; 
    if L>99999,
        breakcode = 7;
        break
    end
    XX(L)=X;
    YY(L)=Y;
    ZZ(L)=Z;
    RYZ=Y^2+Z^2;
    R2=X^2+RYZ;
    R=sqrt(R2);
    
    % check if the line hit the outer tracing boundary; if yes, then terminate
    % the tracing (label 8):
    if (R > RLIM) | (RYZ > 1600.D0) | (X > 20.D0) ,
        breakcode = 8;
        break
    end
    % check whether or not the inner tracing boundary was crossed from outside,
    % if yes, then calculate the footpoint position by interpolation (go to label 6):
    if (R < R0) & (RR > R) ,
        breakcode = 6;
        break
    end
    
    % check if (i) we are moving outward, or (ii) we are still sufficiently
    % far from Earth (beyond R=5Re); if yes, proceed further:
    if ~ ( (R >= RR) | (R > 5.) ),
        
        % now we moved closer inward (between R=3 and R=5); go to 3 and begin logging
        % previous values of X,Y,Z, to be used in the interpolation (after having
        % crossed the inner tracing boundary):
        
        if ~ (R >= 3.) ,
            % we entered inside the sphere R=3: to avoid too large steps (and hence inaccurate
            % interpolated position of the footpoint), enforce the progressively smaller
            % stepsize values as we approach the inner boundary R=R0:
            FC=0.2;
            if (R-R0 < 0.05), FC=0.05; end
            AL=FC*(R-R0+0.2);
            DS=DIR*AL;
        else
            DS=DIR; 
        end
        XR=X; 
        YR=Y;
        ZR=Z;
    end
    RR=R; 
    [X,Y,Z] = GEOPACK_STEP(X,Y,Z,DS,ERR,IOPT,PARMOD,EXNAME,INNAME);
end 

% find the footpoint position by interpolating between the current and previous
% field line points:
if breakcode <=6,
    R1=(R0-R)/(RR-R); % 6
    X=X-(X-XR)*R1;
    Y=Y-(Y-YR)*R1;
    Z=Z-(Z-ZR)*R1;
else
    if breakcode <= 7,
      msgbox('trace of fieldline was unsuccessful','info','warn');
      L=99999;
    end
end
XF=X; 
YF=Y;
ZF=Z;
return;

