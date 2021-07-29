%{
This routine is used to read MMS/MEC-ephts04d data:
mms[]_mec_[]_l[]_ephts04d_*.cdf

Xingyu Li
2020.10.24
%}
function [Info,Ti,L_dipole,MLT,MLAT,R_GSM]=read_mms_mec(st_time,ed_time,SC,mode)
SC_str=num2str(SC);
level='2';
date_str=datestr(st_time,31);
year=date_str(1:4);
month=date_str(6:7);
day=date_str(9:10);
path=['D:\Research\MMS\Data\mms',SC_str,'\mec\',mode,'\l',level,'\ephts04d\',year,'\',month,'\'];
dirOutput=dir(fullfile(path,['mms',SC_str,'_mec_',mode,'_l',level,'_ephts04d_',year,month,day,'*.cdf']));
filename=[path,dirOutput.name];
% read the file
Info=spdfcdfinfo(filename);
time=spdfcdfread(filename,'variable',{'Epoch'});
L_dipole_all=spdfcdfread(filename,'variable',{['mms',SC_str,'_mec_l_dipole']});
MLT_all=spdfcdfread(filename,'variable',{['mms',SC_str,'_mec_mlt']});
MLAT_all=spdfcdfread(filename,'variable',{['mms',SC_str,'_mec_mlat']});
R_GSM_all=spdfcdfread(filename,'variable',{['mms',SC_str,'_mec_r_gsm']});
D=find(time < ed_time & time > st_time);
Ti=time(D);
L_dipole=L_dipole_all(D);
MLT=MLT_all(D);
MLAT=MLAT_all(D);
R_GSM=R_GSM_all(D,:);