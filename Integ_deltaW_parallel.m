clear all;

%%
particle_str = 'proton';
RE = 6371000; %Radius of the earth; unit in m
e = 1.6e-19; %Unit in C
q = e;
mp = 1.6726231e-27; %mass of proton; unit in kg;
m_e = mp/e;

%%
% solar wind parameter
sttime = datenum(2017,1,20,12,45,00);
entime = datenum(2017,1,20,12,46,00);
[time, w, by, bz, pressure, dst] = get_TS05_parameters(sttime, entime);

GEOPACK_RECALC(2017,20,12,45,0);

PARMOD = zeros(10,1);
PARMOD(1) = pressure; % Solar Wind Ram Pressure, nPa
PARMOD(2) = dst; % Dst, nT
PARMOD(3) = by; % By, GSM, nT
PARMOD(4) = bz; % Bz, GSM, nT

[Info,Ti,L,MLT,MLAT,R_GSM]=read_mms_mec(sttime,entime,1,'srvy');

initialx = R_GSM(1)*1000.0;
initialy = R_GSM(2)*1000.0;
initialz = R_GSM(3)*1000.0;

%% Find Mag-equator MLAT
kext = onera_desp_lib_kext(11);

inoptions = [0, 0, 0, 0, 0];
options  = onera_desp_lib_options(inoptions);

sysaxes = onera_desp_lib_sysaxes(2);

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

dip = -10; % deg

initial_point = find(abs(rad2deg(atan(ZZ./(sqrt(XX.^2+YY.^2))))-rad2deg(atan(ZZ(equator_index)/(sqrt(XX(equator_index)^2+YY(equator_index)^2))))-dip)...
    ==min(abs(rad2deg(atan(ZZ./(sqrt(XX.^2+YY.^2))))-rad2deg(atan(ZZ(equator_index)/(sqrt(XX(equator_index)^2+YY(equator_index)^2))))-dip)));

initialx = XX(initial_point)*RE;
initialy = YY(initial_point)*RE;
initialz = ZZ(initial_point)*RE;

[phi_initial,lambda_initial,~] = cart2sph(initialx,initialy,initialz);

fprintf('Initial Point Latitude : %f\n', rad2deg(lambda_initial));

%%
% [Bx0,By0,Bz0] = tsycalc(initialx/RE,initialy/RE,initialz/RE,'GEOPACK_IGRF_GSM','T01',PARMOD,2017,20,12,45,0);
% Bx0 = Bx0/1.0e9; By0 = By0/1.0e9; Bz0 = Bz0/1.0e9;
% Btot = norm([Bx0, By0, Bz0]);
Btot = B(initial_point)/1.0e9;
step = mp/Btot/q/100;

%%
Em = 7.5; % wave amplitude mV/m
% m = 30; % azimuthal wave number
m = 45:1:75;%[-55:5:-35,30:5:70];%-35:10:55;%5:5:50;
omega = 2.0*pi/360; % wave frequency
tau1 = 45*60.0; % growing characteristic time
tau2 = 45*60.0; % attenuation characteristic time
% E = Em*cos(lambda)*cos(m*phi-omega*t)*exp(t-t_peak/tau)(t<t_peak)(or *exp(-t+t_peak/tau)(t>t_peak))
peak = 12*60.0; % peak position of the electric field

%%
W = [13.8933,14.9033,15.9866,17.1487,18.3953]; %[4.5207,4.8494,5.2019,5.5801,5.9858]; %[10.4933,11.2560,12.0742,12.9519,13.8933]; %[5.9858,6.4208,6.8875,7.3882,7.9252]; %[18.3953,19.7325,21.1668,22.7054,24.3558];%[7.9252,8.5013,9.1193,9.7822,10.4933]; % keV
PA = 0:3:180; % deg
% PA = [78.75,90,101.25];
% PA = 90:11.25:180;

t0 = 0:10:1440; % total: 1440/60 = 24 min

% tau = tau + max(t0);

W_num = length(W);
PA_num = length(PA);
t0_num = length(t0);
m_num = length(m);

%% Parallel
core = 12;
p = parpool(core);
p.IdleTimeout = 10*60; % 10 hours

for ii = 1:W_num
    parfor jj = 1:PA_num

        orbit=[];
        x=[]; y=[]; z=[];
        vx=[]; vy=[]; vz=[];
        phi=[]; lambda=[]; vphi=[];
        t=[]; E=[];

    %         orbit = load(['D:\Research\MMS\Work\Matlab\Test_Particle_Orbit\motion_of_', num2str(W(ii)), 'keV_', num2str(PA(jj)), 'deg_', particle_str, '.txt']);
        orbit = load(['D:\Orbit_3h_-10deg\motion_of_', num2str(W(ii)), 'keV_', num2str(PA(jj)), 'deg_', particle_str, '.txt']);

        x = orbit(:,1); y = orbit(:,2); z = orbit(:,3);
        vx = orbit(:,4); vy = orbit(:,5); vz = orbit(:,6); % km/s

        [phi,lambda,~] = cart2sph(x,y,z);
        [vphi,~,~] = cart2sph(vx,vy,vz);

        for im = 1:m_num
            for kk=1:t0_num
                t = t0(kk)-step.*(0:(length(x)-1))';
    %             E = Em.*cos(lambda-lambda_equator).*cos(m(ii).*phi-omega.*t).*exp(-1.0.*(t-0.0*max(t0)).^2/tau^2);
                E = Em.*cos(lambda-lambda_equator).*cos(m(im).*phi-omega.*t).*exp((t-peak)./tau1.*(t<peak)).*exp(-1.0.*(t-peak)./tau2.*(t>peak));
                delta_W(ii,jj,kk,im)=dot(E,vphi).*step; % Energy, PA, t0, m   unit: eV
            end
            fprintf('m = %f\n', m(im));
        end
        fprintf('Pitch Angle = %f\n', PA(jj));
    end
    fprintf('Energy = %f\n', W(ii));
end
    
delete(p);

% delta_W = squeeze(delta_W);
% save('delta_W.txt','delta_W','-ascii');
% 

% %% Phase Mixing
% for ii=1:m_num
%     for kk=1:t0_num
%         delta_W_sm(ii,:,kk)=smooth(delta_W(ii,:,kk),2);
%     end
% end

for ii=1:W_num
    for jj=1:PA_num
        for im=1:m_num
            delta_W_plot(ii,jj,:,im) = delta_W(ii,jj,:,im)./max(delta_W(ii,jj,:,im));
        end
    end
end
%% Plot
%% ***********************************
% *** set plot position ***
% total panel number.
fig_num=2;
% the interval height between panels.
dh=0.015;%0.01;
% start height position.
start_height_posi=0.15;
% end height position.
end_hight_posi=0.95;
% Note: here we may modify it to plot multi column figures. 
% start left position.
left_posi=0.15;
% panel width.
width=0.65;
% calculate panel height and bottom array.
height=(end_hight_posi-start_height_posi-(fig_num-1)*dh)/fig_num;
bottom_posi=start_height_posi+(height+dh)*(0:1:(fig_num-1));
% calculate positon cell array.
position_cell=cell(fig_num,1);
for i_posi=1:fig_num
    position_cell{i_posi}=[left_posi bottom_posi(fig_num-i_posi+1) width height];
end
% colorbar postion
position_bar_cell=cell(fig_num,1);
for i_posi=1:fig_num
    position_bar_cell{i_posi}=[left_posi+width+0.01 bottom_posi(fig_num-i_posi+1) 0.01 height];
end

%%
CT_index=35;
[cmap_IDL,~] = get_CT_IDL(CT_index);
cmap_IDL(1:10,:)=1;

% for im = 1:m_num
%     figure(im);
%     figure('color',[1 1 1]);
% 
%     ip=1;
%     subplot(fig_num,1,ip);
%     plot(t0./60.0,Em*cos(lambda_initial-lambda_equator).*cos(m(im)*phi_initial-omega.*t0).*exp((t0-peak)./tau1.*(t0<peak)).*exp(-1.0.*(t0-peak)./tau2.*(t0>peak)));
%     title(['m = ',num2str(m(im))]);
%     set(gca,'xlim',[0 max(t0./60.0)],'ylim',[-8 8],'xtick',0:4:24,'xticklabel','','tickdir','in','Position',position_cell{ip});
%     ylabel('Ea [mV/m]');
% 
%     ip=2;
%     subplot(fig_num,1,ip);
%     imagesc(t0./60.0,PA,squeeze(delta_W_plot(1,:,:,im)));
%     set(gca,'xlim',[0 max(t0./60.0)],'xtick',0:4:24,'ydir','normal','tickdir','out','Position',position_cell{ip});
%     xlabel('Time [min]');
%     ylabel('PA [deg]');
%     set(gcf,'Colormap',cmap_IDL);
%     colormap(redblue);
%     cbar_axes = colorbar;
%     set(cbar_axes,'Position',position_bar_cell{ip});
%     caxis([-1,1]);
%     ylabel(cbar_axes,'PA-normalized \delta W','Rotation',90);
% end

%% Phase Mixing
delta_W_mixW = squeeze(mean(delta_W(1:1:5,:,:,:),1));

for jj=1:PA_num
    for im=1:m_num
        delta_W_mixW_plot(jj,:,im) = delta_W_mixW(jj,:,im)./max(delta_W_mixW(jj,:,im));
    end
end

% delta_W_mixW_mixPA = (delta_W_mixW(1:2:31,:,:)+delta_W_mixW([2:2:14,16,16:2:30],:,:))./2;
delta_W_mixW_mixPA = (delta_W_mixW(1:4:61,:,:)+delta_W_mixW([2:4:30,32:4:60],:,:)+...
    delta_W_mixW([3:4:31,31:4:59],:,:)+delta_W_mixW([4:4:28,30,32,34:4:58],:,:))./4.0;

PA_mix = 0:12:180;
PA_mix_num = length(PA_mix);

for jj=1:PA_mix_num
    for im=1:m_num
        delta_W_mix_plot(jj,:,im) = delta_W_mixW_mixPA(jj,:,im)./max(delta_W_mixW_mixPA(jj,:,im));
    end
end

for im = 1:m_num
    figure(im);
    figure('color',[1 1 1]);

    ip=1;
    subplot(fig_num,1,ip);
    plot(t0./60.0,Em*cos(lambda_initial-lambda_equator).*cos(m(im)*phi_initial-omega.*t0).*exp((t0-peak)./tau1.*(t0<peak)).*exp(-1.0.*(t0-peak)./tau2.*(t0>peak)));
    title(['m = ',num2str(m(im))]);
    set(gca,'xlim',[3 max(t0./60.0)*(3/4)+3],'ylim',[-8 8],'xtick',3:6:max(t0./60.0)*(3/4)+3,'xticklabel','','tickdir','in','Position',position_cell{ip});
    ylabel('Ea [mV/m]');

    ip=2;
    subplot(fig_num,1,ip);
    imagesc(t0./60.0,PA_mix,squeeze(delta_W_mix_plot(:,:,im)));
%     imagesc(t0./60.0,PA,squeeze(delta_W_mixW_plot(:,:,im)));
%     imagesc(t0./60.0,PA,squeeze(delta_W_plot(1,:,:,im)));
    set(gca,'xlim',[3 max(t0./60.0)*(3/4)+3],'xtick',3:6:max(t0./60.0)*(3/4)+3,'xticklabel',{'0','6','12','18'},'ydir','normal','tickdir','out','Position',position_cell{ip});
    xlabel('Time [min]');
    ylabel('PA [deg]');
    set(gcf,'Colormap',cmap_IDL);
    colormap(redblue);
    cbar_axes = colorbar;
    set(cbar_axes,'Position',position_bar_cell{ip});
    caxis([-1,1]);
    ylabel(cbar_axes,'PA-normalized \delta W','Rotation',90);
end