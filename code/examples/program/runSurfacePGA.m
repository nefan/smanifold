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

function [Vapprox Vexact sapprox sfletcher sexact angularDiff] = runSurfacePGA(setupfile,outputDir,tmpDir,nrProcesses,exitFile,varargin)
%
% run surface pga
%

addpath(genpath('code'));
[sucess,message] = mkdir(tmpDir);

epsilon = 10e-5;

% show figures?
visible = true;
if visible
    set(0, 'DefaultFigureVisible', 'on')
else
    set(0, 'DefaultFigureVisible', 'off')
end

if exist(outputDir,'dir') == 0
    mkdir(outputDir);
end
assert(exist(outputDir,'dir') > 0);

% multicore setup
global multicoreSettings;
multicoreSettings.conf.multicoreDir = tmpDir;
multicoreSettings.conf.masterIsWorker = 1;
multicoreSettings.conf.useWaitbar = false;
multicoreSettings.conf.maxEvalTimeSingle = 10000000;
[status,result] = system('echo $OMPI_MCA_ns_nds_vpid'); 
multicoreSettings.processId = sscanf(result,'%d');
multicoreSettings.nrProcesses = sscanf(nrProcesses,'%d');
multicoreSettings.exitFile = exitFile;

if isempty(multicoreSettings.processId) % environment not setup for MPI calculations
    assert(multicoreSettings.nrProcesses == 1);
    multicoreSettings.processId = 0;
end

if startmulticore(multicoreSettings) % if true we exit when done
    return;
end

% get setup
evalFileName = setupfile;
evalFile;
% and check...
allDefined = ...
    exist('name') && ...    
    exist('m') && ...
    exist('n') && ...
    exist('C') && ...
    exist('p') && ...
    exist('data2d');
assert(allDefined,'incorrect setup: missing variables...');
global prepend
prepend = fullfile(outputDir,[name '-']);
printFig = @(name) print([prepend name '.ps'],'-dps');

% number of samples
N = size(data2d,2);
multicoreSettings.nrOfEvalsAtOnce = ceil(N/multicoreSettings.nrProcesses);

format short e;

% setup manifold: quadratic surface
dimM = m-n;
F = @(x) sum(C'.*(x.^2))-1.0;
DF = @(x) 2*C.*x';
D2F = @(x) 2.0*diag(C);
manifold = embeddedManifold(m,n,F,DF,D2F,tol);

% tangent space
if ~exist('B','var')
	B = manifold.orthFrame(p);
else
	assert(manifold.isTangent(B,p));
end

% data
figure(2)
clf
plot(data2d(1,:),data2d(2,:),'r*','MarkerSize',10);
axis equal

dataM = [];
for i = 1:N
    x2 = data2d(:,i);
    [x v] = manifold.Exp(p,B*x2);
    y3 = x;
        
    dataM(:,i) = y3;
end

if visible
    figure(1)
    clf
    %xrange = [-0.5+min(dataM(1,:)),0.5+max(dataM(1,:))];
    %yrange = [-0.5+min(dataM(2,:)),0.5+max(dataM(2,:))];
    %zrange = [-0.5+min(dataM(3,:)),0.5+max(dataM(3,:))];
    xrange = [-2.5,2.5];
    yrange = [-2.5,2.5];
    zrange = [-0.5,2.5];
    ImplicitPlot3D(F,xrange,yrange,zrange);
    hold on, plot3(dataM(1,:),dataM(2,:),dataM(3,:),'ro','MarkerSize',10,'MarkerFaceColor','r');
    hold off
    axis([-2.5 2.5 -2.5 2.5 0 2.5])
    view([16 22])
end

% example for article
if false
    % curvature
    curvature(p,[0; 1; 0],[1; 0; 0],m,Exp,dExp,tol);
    
    % injectivity radius
    [xx vv solExp] = Exp(p,[0; pi; 0]);
    [BB solDExp] = dExp(solExp,[1; 0; 0]);
    ll = []
    %tt = 0:0.06:1; % for the sphere plot
    tt = 0:0.05:1;
    for i = 1:length(tt)
        t = tt(i);
        ll(i) = norm(getDExp(solDExp,m,t));
    end
    figure(6)
    hold on, plot(tt*pi,ll);
    
    figure(1)
    hold on,    
    pp = [];
    for i = 1:length(tt)
        t = tt(i);
        p = getExp(solExp,DF,t);
        pp(:,i) = p;
        v = 2.0*getDExp(solDExp,m,t);
        myquiver(p(1),p(2),p(3),v(1),v(2),v(3));
    end
    plot3(pp(1,:),pp(2,:),pp(3,:),'k','LineWidth',3,'MarkerFaceColor','k');
    plot3([0 0],[0 0],[1 -1],'ro','MarkerSize',10,'MarkerFaceColor','r');
    view(-150,30);
end

% PGA
[Vexact Vapprox sexact sapprox sfletcher angularDiff] = runPGA(dataM,p,2,tol,manifold);

% visualization
% approximated PGA
v = Vapprox(:,1);
if visible
    figure(1)
    hold on
    quiver3(p(1),p(2),p(3),v(1,1),v(2,1),v(3,1),'g','LineWidth',3);
    hold off
end
figure(2)
hold on
v = B'*v;
if visible
    quiver(p(1),p(2),v(1,1),v(2,1),'g','LineWidth',3);
else
    quiver(p(1),p(2),v(1,1),v(2,1),'g');
end
hold off

% real PGA
v = Vexact(:,1);
if visible
    figure(1)
    hold on
    quiver3(p(1),p(2),p(3),v(1,1),v(2,1),v(3,1),'b','LineWidth',3);
    hold off
end
figure(2)
hold on
v = B'*v;
if visible
    quiver(p(1),p(2),v(1,1),v(2,1),'b','LineWidth',3);
else
    quiver(p(1),p(2),v(1,1),v(2,1),'b');
end
hold off
printFig('TM');
if ~visible
    close(2);
end

endmulticore(multicoreSettings);

end



