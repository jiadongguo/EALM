clc
clear
close all
% grid cells
par.nz=92;par.nx=334;
par.dz=10;par.dx=10;
% boundary
par.top=30;par.bot=30;par.lft=30;par.rht=30;
par.nzb=par.nz+par.top+par.bot;par.nxb=par.nx+par.lft+par.rht;
% time
par.nt=2000;par.dt=0.001;
% source position
par.sx=par.nx/2;par.sz=10;
% 地震子波
wt=ricker(par,10,1,0.1);
figure,
plot(wt);
title("雷克子波");
%速度模型
fd=fopen("vel","r");
vel=fread(fd,[par.nz par.nx],'float');
fclose(fd);
figure,
imagesc(vel);
title("速度模型");
colormap("jet");
clim([1500 5500]);
%%
vv=pad2(par,vel);
figure,
imagesc(vv);
title("填充速度模型");
colormap("jet");
clim([2500 4500]);
%%
% 吸收边界
d=absorb(par,vv);
% d(:)=0;
figure,
imagesc(d);
colormap("jet");
%%
%main loop
rcd=zeros(par.nt,par.nx);
curr=zeros(par.nzb,par.nxb);
pre=zeros(par.nzb,par.nxb);
figure,
h=imagesc(curr);
colormap("jet");
clim([-5 5]);
for it=1:par.nt
    next=fdfor2(par,pre,curr,vv,d);
    next=addsource(par,next,wt(it));
    pre=curr;curr=next;
    rcd(it,:)=curr(par.top+10,par.lft+1:par.lft+par.nx);
    set(h,'CData',curr);
    pause(0.01);
end
figure,
imagesc(rcd);
colormap("jet");
clim([-10 10]);
%%
function vv=pad2(par,v)
[nz,nx]=size(v);
nzb=par.nzb;nxb=par.nxb;
lft=par.lft;top=par.top;
rht=par.rht;bot=par.bot;
vv=zeros(nzb,nxb);
vv(top+1:top+nz,lft+1:lft+nx)=v;
vv(1:top,lft+1:lft+nx)=repmat(v(1,:),[top,1]);
vv(top+nz+1:nzb,lft+1:lft+nx)=repmat(v(nz,:),[bot,1]);
vv(1:nzb,1:lft)=repmat(vv(:,lft+1),[1,lft]);
vv(1:nzb,lft+nx+1:nxb)=repmat(vv(:,lft+nx),[1,rht]);
end
function next=fdfor2(par,pre,curr,vv,d)
nz=par.nz;nx=par.nx; nzb = par.nzb;nxb = par.nxb;top = par.top;
lft = par.lft;dz = par.dz;dx = par.dx;dt = par.dt;
ix=3:nxb-2;iz=3:nzb-2;
c0=-5/2;c1=4/3;c2=-1/12;
tmp=d(iz,ix)*dt;
% 拉普拉斯算子
lap=(c0*curr(iz,ix)+c1*(curr(iz+1,ix)+curr(iz-1,ix))+c2*(curr(iz+2,ix)+curr(iz-2,ix)))./(dz^2)...
+(c0*curr(iz,ix)+c1*(curr(iz,ix+1)+curr(iz,ix-1))+c2*(curr(iz,ix+2)+curr(iz,ix-2)))./(dx^2);
% 第一项
A=(2-tmp.^2)./(1+tmp).*curr(iz,ix);
% 第二项
B=(dt*vv(iz,ix)).^2./(1+tmp).*lap;
% 第三项
C=(tmp-1)./(tmp+1).*pre(iz,ix);

next=zeros(nzb,nxb);
next(iz,ix)=A+B+C;
end
function p=addsource(par,p,wt)
sx=par.sx+par.lft;sz=par.sz+par.top;
p(sz,sx)=p(sz,sx)+wt;
end
function d=absorb(par,vv)
r=log(1e4);
nz=par.nz;nx=par.nx;
nzb=par.nzb;nxb=par.nxb;
lft=par.lft;top=par.top;
rht=par.rht;bot=par.bot;
dx=par.dx;dz=par.dz;
d=zeros(nzb,nxb);
for ix=1:nxb
    for iz=1:nzb
        dxx=0;dzz=0;
        % left
        if ix<=lft
            l=(lft-ix+1)*dx;
            L=lft*dx;
            dxx=r*1.5*vv(iz,ix)/L*(l/L)^2;
        end
        % right
        if ix>lft+nx
            l=(ix-lft-nx)*dx;
            L=rht*dx;
            dxx=r*1.5*vv(iz,ix)/L*(l/L)^2;
        end
        % top
        if iz<=top
            l=(top-iz+1)*dz;
            L=lft*dz;
            dzz=r*1.5*vv(iz,ix)/L*(l/L)^2;
        end
        % bot
        if iz>top+nz
            l=(iz-top-nz)*dz;
            L=bot*dz;
            dzz=r*1.5*vv(iz,ix)/L*(l/L)^2;
        end
        d(iz,ix)=hypot(dxx,dzz);
    end
end
end
function wt=ricker(par,fm,amp,t0)
nt=par.nt;dt=par.dt;
t=0:dt:(nt-1)*dt;t=(t-t0)*fm*pi;
t=t.^2;wt=amp.*(1-2*t).*exp(-t);
end