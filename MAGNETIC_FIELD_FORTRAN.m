function [BX, BY, BZ] = MAGNETIC_FIELD_FORTRAN(PARMOD,XGSM,YGSM,ZGSM)
    % set global GEOPACK variables
    global GEOPACK1;
    % wrapper
    t04 = @(parmod,ps,x,y,z)igrft04(1,parmod,ps,x,y,z);
    igrfgsw = @(x,y,z)igrft04(3,x,y,z);
    gswgse = @(x,y,z,d)igrft04(8,x,y,z,d);
    
    [BXGSM,BYGSM,BZGSM] = t04(PARMOD, GEOPACK1.PSI, XGSM, YGSM, ZGSM);
    
    [XGSE, YGSE, ZGSE] = GEOPACK_GSMGSE(XGSM, YGSM, ZGSM, 1);
    [XGSW, YGSW, ZGSW] = gswgse(XGSE, YGSE, ZGSE, -1);
    [HXGSW, HYGSW, HZGSW] = igrfgsw(XGSW, YGSW, ZGSW);
    [HXGSE, HYGSE, HZGSE] = gswgse(HXGSW, HYGSW, HZGSW, 1);
    [HXGSM, HYGSM, HZGSM] = GEOPACK_GSMGSE(HXGSE, HYGSE, HZGSE, -1);
    
    BX=BXGSM+HXGSM;
    BY=BYGSM+HYGSM;
    BZ=BZGSM+HZGSM;
end