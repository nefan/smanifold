%  smanifold, Copyright (C) 2009-2012, Stefan Sommer (sommer@diku.dk)
%  https://github.com/nefan/smanifold.git
% 
%  This file is part of smanifold.
% 
%  smanifold is free software: you can redistribute it and/or modify
%  it under the terms of the GNU General Public License as published by
%  the Free Software Foundation, either version 3 of the License, or
%  (at your option) any later version.
% 
%  smanifold is distributed in the hope that it will be useful,
%  but WITHOUT ANY WARRANTY; without even the implied warranty of
%  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%  GNU General Public License for more details.
% 
%  You should have received a copy of the GNU General Public License
%  along with smanifold.  If not, see <http://www.gnu.org/licenses/>.
%  

%
% example for running PGA on a 4D quadratic manifold
%

setupfile = 'examples/setupfiles/qhcatest.m';
outputDir = 'tmp/output';
tmpDir = 'tmp/';

[Vapprox,Vexact,PGAVapprox,sapprox,sexact,PGAsapprox,coordsapprox,coordsexact,PGAcoordsapprox,dataM] = runQuadraticHCA(setupfile,outputDir,tmpDir);

% plots for paper
% get setup
evalFileName = setupfile;
evalFile;

% plot
figure(1), clf
coords = Vapprox*coordsapprox;
coords = [coords(2,:); coords(3,:); coords(1,:)];
plot3(coords(1,:),coords(2,:),coords(3,:),'ro','MarkerSize',6,'MarkerFaceColor','r')
axis([-1.5 1.5 -1.5 1.5 -1.5 1.5])
grid
xx=xlim; yy=ylim; zz=zlim
arrow([xx(1) 0 0],[xx(2) 0 0],'Length',10), arrow fixlimits
arrow([0 yy(1) 0],[0 yy(2) 0],'Length',10), arrow fixlimits
arrow([0 0 zz(1)],[0 0 zz(2)],'Length',10), arrow fixlimits
text(xx(2),0.1,0,'x_2','fontsize',20,'horizontalalignment','center');
text(0,yy(2),0.15,'x_3','fontsize',20,'horizontalalignment','center');
text(0,0.1,zz(2),'x_1','fontsize',20,'horizontalalignment','center');
view(73,24)
 
% plot
figure(2), clf
coords = PGAVapprox*PGAcoordsapprox;
coords = [coords(2,:); coords(3,:); coords(1,:)];
plot3(coords(1,:),coords(2,:),coords(3,:),'ro','MarkerSize',6,'MarkerFaceColor','r')
axis([-1.5 1.5 -1.5 1.5 -1.5 1.5])
grid
xx=xlim; yy=ylim; zz=zlim
arrow([xx(1) 0 0],[xx(2) 0 0],'Length',10), arrow fixlimits
arrow([0 yy(1) 0],[0 yy(2) 0],'Length',10), arrow fixlimits
arrow([0 0 zz(1)],[0 0 zz(2)],'Length',10), arrow fixlimits
text(xx(2),0.1,0,'x_2','fontsize',20,'horizontalalignment','center');
text(0,yy(2),0.15,'x_3','fontsize',20,'horizontalalignment','center');
text(0,0.1,zz(2),'x_1','fontsize',20,'horizontalalignment','center');
view(73,24)

% plot
figure(2), hold on
v = -PGAVapprox(:,2);
v = [v(2) v(3) v(1)];
arrow([0 0 0],v,'LineWidth',5,'EdgeColor','r','FaceColor','r');
figure(1), hold on
v = -Vapprox(:,2);
v = [v(2) v(3) v(1)];
dd = d1;
dd = [dd(2) dd(3) dd(1)];
arrow(dd,dd+v,'LineWidth',5,'EdgeColor','b','FaceColor','b');
dd = d2;
dd = [dd(2) dd(3) dd(1)];
arrow(dd,dd+v,'LineWidth',5,'EdgeColor','b','FaceColor','b');

figure(3)
clf
%xrange = [-0.5+min(dataM(1,:)),0.5+max(dataM(1,:))];
%yrange = [-0.5+min(dataM(2,:)),0.5+max(dataM(2,:))];
%zrange = [-0.5+min(dataM(3,:)),0.5+max(dataM(3,:))];
xrange = [-2.5,2.5];
yrange = [-2.5,2.5];
zrange = [-2.5,2.5];

II = [2 3 4];
CC = C(II);
FF = @(x) sum(CC'.*(x.^2))-1.0;
DFF = @(x) 2*CC.*x';
D2FF = @(x) 2.0*diag(CC);
manifold = embeddedManifold(m-1,n,FF,DFF,D2FF,tol);
ImplicitPlot3D(FF,xrange,yrange,zrange,50);
hold on
dataMM = dataM(II,:);
for i=1:size(dataMM,2)
    dataMM(:,i) = manifold.toManifold(dataMM(:,i));
end
plot3(dataMM(1,:),dataMM(2,:),dataMM(3,:),'ro','MarkerSize',6,'MarkerFaceColor','r');
plot3(p(II(1)),p(II(2)),p(II(3)),'ko','MarkerSize',10,'MarkerFaceColor','k');
hold off
axis([-2.5 2.5 -2.5 2.5 -2.5 2.5])
view([51 26])
xx=xlim; yy=ylim; zz=zlim
arrow([xx(1) 0 0],[xx(2) 0 0],'Length',10), arrow fixlimits
arrow([0 yy(1) 0],[0 yy(2) 0],'Length',10), arrow fixlimits
arrow([0 0 zz(1)],[0 0 zz(2)],'Length',10), arrow fixlimits
text(xx(2),0.1,0,'x_2','fontsize',20,'horizontalalignment','center');
text(0,yy(2),0.15,'x_3','fontsize',20,'horizontalalignment','center');
text(0,0.1,zz(2),'x_4','fontsize',20,'horizontalalignment','center');
zoom(0.85)
