clc
clear
close all
nz=81;nx=201;nt=2001;dt=0.0005;dx=10;dz=10;nb=50;
nzb=nz+2*nb;nxb=nx+2*nb;

fd=fopen('vel.dat','r');
vel=fread(fd,[nz nx],'float');
fclose(fd);
figure,
imagesc(vel);
colormap('jet');





fd=fopen('snap.dat','r');
rcd=fread(fd,[nt nx],'float');
fclose(fd);
x=(0:nx-1)*dx;t=(0:nt-1)*dt;
figure,
imagesc(x,t,rcd);
colormap('gray');
clim([-1.5 1.5]);



fd=fopen('wt.dat','r');
wt=fread(fd,[nt 1],'float');
fclose(fd);
figure,
plot(wt);


fd=fopen('snap.dat','r');
seis=fread(fd,[nzb nxb],'float');
figure,
h=imagesc(seis);
colormap('jet');
clim([-0.5 0.5]);
for it=2:nt
pause(0.01);
seis=fread(fd,[nzb nxb],'float');
set(h,'Cdata',seis);
end