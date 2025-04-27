clc
clear
close all
nz=301;nx=941;nt=2001;dt=0.0005;dx=10;dz=10;nb=30;
nzb=nz+2*nb;nxb=nx+2*nb;
fd=fopen('wavefield_1.dat','r');
rcd=fread(fd,[nt nx],'float');
fclose(fd);
x=(0:nx-1)*dx;t=(0:nt-1)*dt;
figure,
imagesc(x,t,rcd);
colormap('jet');
clim([-0.5 0.5]);



fd=fopen('snap.dat','r');
seis=fread(fd,[nzb nxb],'float');
figure,
h=imagesc(seis);
colormap('jet');
clim([-0.5 0.5]);
for it=2:nt
pause(0.02);
seis=fread(fd,[nzb nxb],'float');
set(h,'Cdata',seis);
end