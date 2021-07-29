function plot_trajectory(particle_str, L, E, theta)

e = 1.6e-19;
RE = 6371000; %Radius of the earth; unit in m
%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
%改成你的文件路径
% fid = fopen(['D:\Research\MMS\Work\Matlab\Test_Particle_Orbit\motion of ', particle_str, '.txt'], 'r');
% % line1 = fgetl(fid);
% % line1 = fgetl(fid);
% % line1 = fgetl(fid);
% 
% filelength = ftell(fid);
% 
% xlist = zeros(1, filelength);
% ylist = zeros(1, filelength);
% zlist = zeros(1, filelength);
% %energy = zeros(1, ceil(filelength / gap));
% 
% i = 0;
% while feof(fid) == 0
%     i = i+1;
%     line1 = fgetl(fid);
%     line2 = deblank(line1);
%     line3 = regexp(line2, '\s+', 'split');
%     xlist(i) = str2num(line3{1});
%     ylist(i)  = str2num(line3{2});
%     zlist(i)  = str2num(line3{3});
% 
% end
% 
% fclose(fid);

orbit = load(['F:\Orbit_3h_-10deg\motion_of_', num2str(E/e/1000), 'keV_', num2str(rad2deg(theta)), 'deg_', particle_str, '.txt']);

xlist = orbit(:,1);
ylist = orbit(:,2);
zlist = orbit(:,3);

[x,y,z] = meshgrid([-10:0.1:10],[-10:0.1:10],[-10:0.1:10]);
[Bx,By,Bz] = magnetic_field(x * RE,y * RE,z * RE);
% [x,y,z] = meshgrid([-10:0.1:10],[-10:0.1:10],[-10:0.1:10]);
% for ii = 1:length(x)
%     for jj = length(y)
%         for kk = length(z)
%             [Bx(ii,jj,kk),By(ii,jj,kk),Bz(ii,jj,kk)] = GEOPACK_DIP(x(ii,jj,kk),y(ii,jj,kk),z(ii,jj,kk));
%         end
%     end
% end
% Bx = Bx./1.0e9; By = By./1.0e9; Bz = Bz./1.0e9;
sphere;
hold on;
grid on;
m=12;

sx = zeros(1,m);
sy = zeros(1,m);
sz = zeros(1,m);
for i=1:m
    a = 2*pi/m*i - 2*pi/m/2;
    sx(i) = sin(a)*L;
    sy(i) = cos(a)*L;
    sz(i) = 0;
end
%axis([-5 5 -5 5 -5 5]);
axis([-5 5 -5 5 -5 5] * 2);
% streamline(x,y,z,Bx,By,Bz,sx,sy,sz);
Bx = -Bx;By = -By;Bz = -Bz;
% streamline(x,y,z,Bx,By,Bz,sx,sy,sz);

plot3(xlist, ylist, zlist, 'r');

end