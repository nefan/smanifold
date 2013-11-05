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
% example for running HCA on a surface, bimodal distribution
%

setupfile = 'examples/setupfiles/surfacehca.m';
outputDir = 'tmp/output';
tmpDir = 'tmp/';

[Vapprox,Vexact,PGAVapprox,sapprox,sexact,PGAsapprox,coordsapprox,coordsexact,PGAcoordsapprox,dataM] = runQuadraticHCA(setupfile,outputDir,tmpDir);

% plots for paper
% get setup
evalFileName = setupfile;
evalFile;

V = Vexact;
s2 = sexact(1);
coordsOut = coordsexact;

% plot
figure(1), clf, hold on
coords = V*coordsOut;
coords = coords([2 3],:);
plot(coords(1,:),coords(2,:),'ro','MarkerSize',6,'MarkerFaceColor','r')
axis([-1.5 1.5 -1.5 1.5])
grid
xx=xlim; yy=ylim;
arrow([xx(1) 0],[xx(2) 0],'Length',10), arrow fixlimits
arrow([0 yy(1)],[0 yy(2)],'Length',10), arrow fixlimits

v = sqrt(s2)*V(:,2);
dd = d1;
arrow(dd,dd+v,'LineWidth',5,'EdgeColor','b','FaceColor','b');
dd = d2;
arrow(dd,dd+v,'LineWidth',5,'EdgeColor','b','FaceColor','b');

 
% plot
figure(2), clf, hold on
coords = PGAVapprox*PGAcoordsapprox;
coords = coords([2 3],:);
plot(coords(1,:),coords(2,:),'ro','MarkerSize',6,'MarkerFaceColor','r')
axis([-1.5 1.5 -1.5 1.5])
grid
xx=xlim; yy=ylim;
arrow([xx(1) 0],[xx(2) 0],'Length',10), arrow fixlimits
arrow([0 yy(1)],[0 yy(2)],'Length',10), arrow fixlimits

v = sqrt(PGAsapprox(2))*PGAVapprox(:,2);
arrow([0 0],v,'LineWidth',5,'EdgeColor','r','FaceColor','r');

figure(3), hold on
clf
%xrange = [-0.5+min(dataM(1,:)),0.5+max(dataM(1,:))];
%yrange = [-0.5+min(dataM(2,:)),0.5+max(dataM(2,:))];
xrange = [-1.1,1.1];
yrange = [-1.1,1.1];
zrange = [-1.1,1.1];

II = [1 2 3];
CC = C(II);
FF = @(x) sum(CC'.*(x.^2))-1.0;
DFF = @(x) 2*CC.*x';
D2FF = @(x) 2.0*diag(CC);
manifold = embeddedManifold(m,n,FF,DFF,D2FF,tol);
ImplicitPlot3D(FF,xrange,yrange,zrange,50);
hold on
dataMM = dataM;
plot3(dataMM(1,:),dataMM(2,:),dataMM(3,:),'ro','MarkerSize',6,'MarkerFaceColor','r');
% plot3(p(II(1)),p(II(2)),p(II(3)),'ko','MarkerSize',10,'MarkerFaceColor','k');
axis([-1.1 1.1 -1.1 1.1 -1.1 1.1])
view([80 16])
xx=xlim; yy=ylim; zz=zlim
arrow([xx(1) 0 0],[xx(2) 0 0],'Length',10), arrow fixlimits
arrow([0 yy(1) 0],[0 yy(2) 0],'Length',10), arrow fixlimits
arrow([0 0 zz(1)],[0 0 zz(2)],'Length',10), arrow fixlimits
text(xx(2),0.1,0,'x_1','fontsize',20,'horizontalalignment','center');
text(0,yy(2),0.15,'x_2','fontsize',20,'horizontalalignment','center');
text(0,0.1,zz(2),'x_3','fontsize',20,'horizontalalignment','center');
zoom(0.9)

% geodesic, principal
[xx vv sol] = manifold.Exp(p,2*pi*V(:,1));
xx = [];
for t = [0:0.01:1]
    xx = [xx manifold.getExp(sol,t)];
end
plot3(xx(1,:),xx(2,:),xx(3,:),'--k','LineWidth',2);

% components
scale = 1.0;
[x1 vv sol] = manifold.Exp(p,V(:,1));
v = scale*sqrt(s2)*V(:,2);
arrow(x1,x1+v,'LineWidth',5,'EdgeColor','b','FaceColor','b');
[x1 vv sol] = manifold.Exp(p,-V(:,1));
v = scale*sqrt(s2)*V(:,2);
arrow(x1,x1+v,'LineWidth',5,'EdgeColor','b','FaceColor','b','MarkerEdgeColor','b','MarkerFaceColor','b');
v = -scale*sqrt(PGAsapprox(2))*PGAVapprox(:,2);
arrow(p,p+v,'LineWidth',5,'EdgeColor','r','FaceColor','r');

% geodesic, illustration
w = [-1 -0.6]';
[xx vv sol] = manifold.Exp(p,B*w);
xx = [];
for t = [0:0.01:1]
    xx = [xx manifold.getExp(sol,t)];
end
plot3(xx(1,:),xx(2,:),xx(3,:),'-.g','LineWidth',2);

