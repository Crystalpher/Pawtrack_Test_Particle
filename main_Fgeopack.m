% This is the main and start function.
clear all;

%%
% solar wind parameter
sttime = datenum(2017,1,20,12,45,00);
entime = datenum(2017,1,20,12,46,00);
[time, w, by, bz, pressure, dst, vx, vy, vz] = get_TS05_parameters(sttime, entime);

GEOPACK_RECALC(2017,20,12,45,0);

% FORTRAN code wrapper & initialize
recalc = @(y,d,h,m,s,vx,vy,vz)igrft04(2,y,d,h,m,s,vx,vy,vz);
recalc(2017,20,12,45,0,vx,vy,vz);

PARMOD = zeros(10,1);
PARMOD(1) = pressure; % Solar Wind Ram Pressure, nPa
PARMOD(2) = dst; % Dst, nT
PARMOD(3) = by; % By, GSM, nT
PARMOD(4) = bz; % Bz, GSM, nT

%%
%Initialization and some parameters

particle_str = 'proton';
RE = 6371000; %Radius of the earth; unit in m
e = 1.6e-19; %Unit in C
q = e;
m = 1.6726231e-27; %mass of proton; unit in kg;
m_e = m/e;
% L = 4.0; %Unit in RE

[Info,Ti,L,MLT,MLAT,R_GSM]=read_mms_mec(sttime,entime,1,'srvy');

% initialx = L*RE;
% initialy = 0;
% initialz = 0;
initialx = R_GSM(1)*1000.0;
initialy = R_GSM(2)*1000.0;
initialz = R_GSM(3)*1000.0;

%% Find Mag-equator MLAT
kext = onera_desp_lib_kext(11); % TS04

inoptions = [0, 0, 0, 0, 0];
options  = onera_desp_lib_options(inoptions);

sysaxes = onera_desp_lib_sysaxes(2); % GSM

matlabd = datenum(2017,1,20,12,45,00);
[iyear,idoy,UT] = onera_desp_lib_matlabd2yds(matlabd);

Kp = 0;
Dst = dst;
Nsw = 0;
Vsw = 400;
Psw = pressure;
ByGSM = by;
BzGSM = bz;
AL = 0;
Bperp = sqrt(by^2+bz^2);
theta = atan(by/bz);
G1 = Vsw*(Bperp/40)^2/(1+Bperp/40)*sin(theta/2)^3;
a = 0.005;
if bz<0
   Bs = abs(bz);
else
    Bs = 0;
end
G2 = a*Vsw*Bs;
G3 = Vsw*Nsw*Bs/2000.0;
maginputs = onera_desp_lib_maginputs(Kp,Dst,Nsw,Vsw,Psw,ByGSM,BzGSM,G1,G2,G3,w(1),w(2),w(3),w(4),w(5),w(6),AL);

% % DIR_1=0.001;
% % DIR_2=-0.001;
% DIR_1=0.01;
% DIR_2=-0.01;
% RLIM=100.0;
% R0=1.0;
% IOPT=10.0; %IOPT is a dummy number for 'T96'

% [XF_1,YF_1,ZF_1,XX_1,YY_1,ZZ_1,L_1] = GEOPACK_TRACE_smallstep(initialx/RE,initialy/RE,initialz/RE,DIR_1,RLIM,R0,IOPT,PARMOD,'T01','GEOPACK_IGRF_GSM');
% [XF_2,YF_2,ZF_2,XX_2,YY_2,ZZ_2,L_2] = GEOPACK_TRACE_smallstep(initialx/RE,initialy/RE,initialz/RE,DIR_2,RLIM,R0,IOPT,PARMOD,'T01','GEOPACK_IGRF_GSM');
% 
% XX=[fliplr(XX_1),XX_2(2:end)];
% YY=[fliplr(YY_1),YY_2(2:end)];
% ZZ=[fliplr(ZZ_1),ZZ_2(2:end)];

[~,B,~,~,POSIT] = onera_desp_lib_trace_field_line(kext,options,sysaxes,matlabd,initialx/RE,initialy/RE,initialz/RE,maginputs,1);

for ii = 1:length(POSIT)
   [XX(ii),YY(ii),ZZ(ii)] = GEOPACK_GEOGSM(POSIT(ii,1),POSIT(ii,2),POSIT(ii,3),1);
end

% [Bgeo,B,gradBmag,diffB] = onera_desp_lib_get_bderivs(kext,options,sysaxes,matlabd,XX,YY,ZZ,maginputs); % nT

for kk=2:length(B)-1
   if B(kk)>B(kk-1)&&B(kk)>B(kk+1) 
       equator_index = kk;
   end
end

[~,lambda_equator,~] = cart2sph(XX(equator_index),YY(equator_index),ZZ(equator_index));

% lambda_equator = deg2rad(20);

% disp(rad2deg(lambda_equator));
fprintf('Magnetic Equator Latitude : %f\n', rad2deg(lambda_equator));

dip = -10;

initial_point = find(abs(rad2deg(atan(ZZ./(sqrt(XX.^2+YY.^2))))-rad2deg(atan(ZZ(equator_index)/(sqrt(XX(equator_index)^2+YY(equator_index)^2))))-dip)...
    ==min(abs(rad2deg(atan(ZZ./(sqrt(XX.^2+YY.^2))))-rad2deg(atan(ZZ(equator_index)/(sqrt(XX(equator_index)^2+YY(equator_index)^2))))-dip)));

initialx = XX(initial_point)*RE;
initialy = YY(initial_point)*RE;
initialz = ZZ(initial_point)*RE;

[~,lambda_initial,~] = cart2sph(initialx,initialy,initialz);

fprintf('Initial Point Latitude : %f\n', rad2deg(lambda_initial));

%%
% [Bx0,By0,Bz0] = GEOPACK_DIP(initialx/RE,initialy/RE,initialz/RE);
% [Bx0,By0,Bz0] = tsycalc(initialx/RE,initialy/RE,initialz/RE,'GEOPACK_IGRF_GSM','T01',PARMOD,2017,20,12,45,0);
[Bx0, By0, Bz0] = MAGNETIC_FIELD_FORTRAN(PARMOD,initialx/RE,initialy/RE,initialz/RE);
Bx0 = Bx0/1.0e9; By0 = By0/1.0e9; Bz0 = Bz0/1.0e9;
% [Bx0,By0,Bz0] = magnetic_field(initialx,initialy,initialz);
Btot = norm([Bx0, By0, Bz0]);
% step = -1.0*m/Btot/q/100;
step = -1.0*m/Btot/q/100;

%%
% [tmpBx,tmpBy,tmpBz] = GEOPACK_DIP(1,0,0);
% [tmpBx,tmpBy,tmpBz] = tsycalc(1,0,0,'GEOPACK_IGRF_GSM','T01',PARMOD,2017,20,12,45,0);
% [tmpBx, tmpBy, tmpBz] = MAGNETIC_FIELD_FORTRAN(PARMOD,1,0,0);
% tmpBx = tmpBx/1.0e9; tmpBy = tmpBy/1.0e9; tmpBz = tmpBz/1.0e9;
% [tmpBx,tmpBy,tmpBz] = magnetic_field(RE,0,0);
% BE = norm([tmpBx, tmpBy, tmpBz]);
% E = 10000*e;
% theta = pi/4;
% Td = -1.0*pi*q*BE*RE^2/3/L/E/(0.35+0.15*sin(theta)); %Drift periods
% Td = -2.0*pi/((6.0*L*E/(e*BE*RE^2))*(0.35+0.15*sin(theta))); %Drift periods
Td = -3.0*3600.0;
n = floor(Td/step);
n0 = ceil(n * 1.0); % total step

%%

% [~,~,Energy_Range]=read_mms_fpi_dist(sttime,entime,1,'i','fast');
% Energy_Range = Energy_Range/1000.0; % keV
% 
% Energy_tabel = exp(linspace(log(Energy_Range(28,1)),log(Energy_Range(28,2)),3)); % Energy of particles: keV
% % 1 keV-30 keV
PA = 0:3:180;%36:6:180;%0:3:180;%[0:6:84,96:6:180];%0:11.25:180; % Pitch angle: deg
Energy = [4.5207,4.8494,5.2019,5.5801];%[5.9858,6.4208,6.8875,7.3882,7.9252];%[18.3953,19.7325,21.1668,22.7054,24.3558];%[7.9252,8.5013,9.1193,9.7822,10.4933];%round(Energy_tabel,4);%7.9252;%8;
% PA = 45;

for ii=1:length(Energy)
    for jj=1:length(PA)
        %% 
        % prepare the initial condition
        % E = 50 * 1000; %Energy of particles; unit in eV
        E = Energy(ii) * 1000 * e;
        % theta = 45; %pitch angle
        theta = deg2rad(PA(jj));

%         E = E*e;
%         theta = theta / 180 * pi;

%         initialvx = sqrt(2*E/m)*sin(theta);
%         initialvy = 0;
%         initialvz = sqrt(2*E/m)*cos(theta);

        V0 = sqrt(2*E/m);

        syms vx0;
        syms vy0;
        syms vz0;
        b=[Bx0,By0,Bz0]/Btot;
        eq1 = b(1)*vx0+b(2)*vy0+b(3)*vz0-cos(theta);
        eq2 = vx0^2+vy0^2+vz0^2-1.0;
        eq3 = b(1)*vy0-b(2)*vx0;
        res = solve(eq1,eq2,eq3,vx0,vy0,vz0);
        initialvx = V0*real(double(res.vx0(1)));
        initialvy = V0*real(double(res.vy0(1)));
        initialvz = V0*real(double(res.vz0(1)));

        %%
        %++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        %改成你的文件路径
        fid = fopen(['F:\Orbit_3h_-10deg\motion_of_', num2str(E/e/1000), 'keV_', num2str(rad2deg(theta)), 'deg_', particle_str, '.txt'], 'w');

        %fprintf(fid, 'Particle: %s; Energy: %d eV\n', particle_str, E/e);
        %fprintf(fid, 'Iniatial x: %d RE; y: %d RE; z: %d RE; vx:%d km/s; vy:%d km/s; vz: %d km/s; pitch_angle: %d\n', ...
        %        initialx/RE, initialy/RE, initialz/RE, initialvx/1000, initialvy/1000, initialvz/1000, theta / pi * 180);

        %fprintf(fid, 'x(RE)  y(RE)  z(RE)  vx(km/s)  vy(km/s)  vz(km/s)  energy(eV)\n');

        Y = [initialx , initialy , initialz , initialvx , initialvy , initialvz];

        fprintf(fid, '%d  %d  %d  %d  %d  %d  %d\n', Y(1)/RE, Y(2)/RE, Y(3)/RE, Y(4)/1000, Y(5)/1000, Y(6)/1000, 0.5 * m_e * (Y(4)^2 + Y(5)^2 + Y(6)^2));
        
        fprintf('Energy : %f keV   Pitch Angle : %f deg\n', Energy(ii), PA(jj));
        
        fprintf('Pitch Angle Error : %d deg\n',rad2deg(acos(dot([Bx0,By0,Bz0],[initialvx,initialvy,initialvz]/(Btot*V0)))-theta));
        
        fprintf('Total Steps : %d\n', n0);
        for i=2:n0

%             K1 = doing_calculation(i*step, Y, q, m, PARMOD);
%             K2 = doing_calculation(i*step + step/2, Y + step/2 .* K1, q, m, PARMOD);
%             K3 = doing_calculation(i*step + step/2, Y + step/2 .* K2, q, m, PARMOD);
%             K4 = doing_calculation(i*step + step, Y + step .* K3, q, m, PARMOD);
%             Y = Y + step./6. * (K1 + 2 .* K2 + 2 .* K3 + K4);
            
%             [Bx,By,Bz] = tsycalc(Y(1)/RE,Y(2)/RE,Y(3)/RE,'GEOPACK_IGRF_GSM','T01',PARMOD,2017,20,12,45,0);
            [Bx, By, Bz] = MAGNETIC_FIELD_FORTRAN(PARMOD,Y(1)/RE,Y(2)/RE,Y(3)/RE);
            Bx = Bx/1.0e9; By = By/1.0e9; Bz = Bz/1.0e9;

            K1 = doing_calculation(i*step, Y, q, m, Bx, By, Bz);
            K2 = doing_calculation(i*step + step/2, Y + step/2 .* K1, q, m, Bx, By, Bz);
            K3 = doing_calculation(i*step + step/2, Y + step/2 .* K2, q, m, Bx, By, Bz);
            K4 = doing_calculation(i*step + step, Y + step .* K3, q, m, Bx, By, Bz);
            Y = Y + step./6 .* (K1 + 2 .* K2 + 2 .* K3 + K4);
            
            fprintf(fid, '%d  %d  %d  %d  %d  %d %d\n', Y(1)/RE, Y(2)/RE, Y(3)/RE, Y(4)/1000, Y(5)/1000, Y(6)/1000, 0.5 * m_e * (Y(4)^2 + Y(5)^2 + Y(6)^2));

            if (mod(i, 20000) == 0)
                fprintf('%f%s \n', i/n0*100.0, '%');
            end
            
            % TS04 model protection
            if Y(1)/RE < -14.99
                disp('TS04 tailwards distance out of range!');
                break;
            end

        end

        fclose(fid);
        
        if PA(jj) == 90
            figure(length(PA)*(ii-1)+jj);
            plot_trajectory(particle_str, L, E, theta);
        end

    end
end