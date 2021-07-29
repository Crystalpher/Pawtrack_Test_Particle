function [time, w, by, bz, pressure, dst, vx, vy, vz] = get_TS05_parameters(sttime, entime)

% sttime = datenum(2018, 8, 2, 0, 0, 0);
% entime = datenum(2018, 8, 3, 0, 0, 0);

strdate = datestr(sttime, 'yyyy-mm-dd');

root = ['D:\MATLABpkgs\TS05\', strdate(1:4), '_OMNI_5m_with_TS05_variables.dat'];

data = load(root);
N = size(data, 1);

year = data(:, 1);
day = data(:, 2);
hour = data(:,3);
min = data(:,4);

time = datenum(year(1), 1, 1) + (day-1) + hour / 24 + min / (60*24);

index = find(time >= sttime & time < entime);

time = time(index);
w = data(index, 18:23);
by = data(index, 6);
bz = data(index, 7);
pressure = data(index, 17);
dst = data(index, 13);
vx = data(index, 8);
vy = data(index, 9);
vz = data(index, 10);

























end