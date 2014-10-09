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

function [dataM,dataTM] = runSurfaceDiff(setupfile,outputDir,tmpDir,varargin)
%
% run surface diffusions
%

addpath(genpath('code'));
[sucess,message] = mkdir(tmpDir);

global epsilon;
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

% get setup
evalFileName = setupfile;
evalFile;
% and check...
allDefined = ...
    exist('name') && ...    
    exist('m') && ...
    exist('n') && ...
    exist('C') && ...
    exist('p');
assert(allDefined,'incorrect setup: missing variables...');
global prepend
prepend = fullfile(outputDir,[name '-']);
printFig = @(name) print([prepend name '.ps'],'-dps');

format short e;

% setup manifold: quadratic surface
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

% number of samples
N = 1000;
K = 20; % Brownian steps

% data
R = chol(covTM);
dataTM = reshape((randn(N*K,2)*R)',2,K,N)/K;
figure(2)
clf, hold on
for i = 1:N
%     plot([0 cumsum(dataTM(1,:,i),2)],[0 cumsum(dataTM(2,:,i),2)],'k-');
    plot(sum(dataTM(1,:,i),2),sum(dataTM(2,:,i),2),'r.','MarkerSize',10);
end
axis([-pi,pi,-pi,pi])
axis equal    
    
dataM = zeros(3,K,N);
pp = p;
parfor i = 1:N
    x = pp;
    Bx = B;
    for k = 1:K
        x2 = dataTM(:,k,i);
        [x v solExp] = manifold.Exp(x,Bx*x2);
        Bx = manifold.Pt(solExp,Bx);

        dataM(:,k,i) = x;
    end
end

% identy path ending close to x
x = manifold.Exp(p,B*xTM);
Ix = [];
for i = 1:N
    if norm(dataM(:,end,i)-x) < r;
        Ix = [Ix i];
    end
end
Ix
meanTMTx = mean(dataTM(:,:,Ix),3);
meanMTx = mean(dataM(:,:,Ix),3);

if visible
    figure(1)
    clf
    %xrange = [-0.5+min(dataM(1,:)),0.5+max(dataM(1,:))];
    %yrange = [-0.5+min(dataM(2,:)),0.5+max(dataM(2,:))];
    %zrange = [-0.5+min(dataM(3,:)),0.5+max(dataM(3,:))];
    xrange = [-2.5,2.5];
    yrange = [-2.5,2.5];
    zrange = [-2.5,2.5];
    ImplicitPlot3D(F,xrange,yrange,zrange,50);
    hold on
    for i = 1:N
%         plot3([p(1) dataM(1,:,i)],[p(2) dataM(2,:,i)],[p(3) dataM(3,:,i)],'k-');
        plot3(dataM(1,end,i),dataM(2,end,i),dataM(3,end,i),'ro','MarkerSize',4,'MarkerFaceColor','r');        
        end
    plot3(p(1),p(2),p(3),'ko','MarkerSize',10,'MarkerFaceColor','k');
    hold off
    axis([-2.5 2.5 -2.5 2.5 -2.5 2.5])
    view([16 22])
    
    figure(4)
    clf, hold on
    for ii = 1:length(Ix)
        i = Ix(ii);
        plot([0 cumsum(dataTM(1,:,i),2)],[0 cumsum(dataTM(2,:,i),2)],'k-');
        plot(sum(dataTM(1,:,i),2),sum(dataTM(2,:,i),2),'r.','MarkerSize',10);
    end
    plot(xTM(1),xTM(2),'go','MarkerSize',20,'MarkerFaceColor','g');
    plot([0 cumsum(meanTMTx(1,:),2)],[0 cumsum(meanTMTx(2,:),2)],'r-','LineWidth',10);
    plot(sum(meanTMTx(1,:),2),sum(meanTMTx(2,:),2),'ko','MarkerSize',20,'MarkerFaceColor','r');
    axis([-pi,pi,-pi,pi])
    axis equal 
    
    figure(3)
    clf
    %xrange = [-0.5+min(dataM(1,:)),0.5+max(dataM(1,:))];
    %yrange = [-0.5+min(dataM(2,:)),0.5+max(dataM(2,:))];
    %zrange = [-0.5+min(dataM(3,:)),0.5+max(dataM(3,:))];
    xrange = [-2.5,2.5];
    yrange = [-2.5,2.5];
    zrange = [-2.5,2.5];
    ImplicitPlot3D(F,xrange,yrange,zrange,50);
    hold on
    for ii = 1:length(Ix)
        i = Ix(ii);
        plot3([p(1) dataM(1,:,i)],[p(2) dataM(2,:,i)],[p(3) dataM(3,:,i)],'k-');
        plot3(dataM(1,end,i),dataM(2,end,i),dataM(3,end,i),'ro','MarkerSize',4,'MarkerFaceColor','r');        
        end
    plot3(p(1),p(2),p(3),'ko','MarkerSize',10,'MarkerFaceColor','k');
    plot3(x(1),x(2),x(3),'go','MarkerSize',10,'MarkerFaceColor','g');
    plot3([p(1) meanMTx(1,:)],[p(2) meanMTx(2,:)],[p(3) meanMTx(3,:)],'r-','LineWidth',5);
    hold off
    axis([-2.5 2.5 -2.5 2.5 -2.5 2.5])
    view([16 22])
    
    figure(5)
    clf, hold on
    [h,canvas]=cloudPlot(cumsum(dataTM(1,:,Ix),2),cumsum(dataTM(2,:,Ix),2),[-.5 2.5 -1.5 1.5],[],[20 20]); colormap gray; colormap(flipud(colormap))
    plot(sum(meanTMTx(1,:),2),sum(meanTMTx(2,:),2),'ko','MarkerSize',20,'MarkerFaceColor','r');
    plot(xTM(1),xTM(2),'go','MarkerSize',20,'MarkerFaceColor','g');
    plot([0 cumsum(meanTMTx(1,:),2)],[0 cumsum(meanTMTx(2,:),2)],'r-','LineWidth',10);
    
end


end

