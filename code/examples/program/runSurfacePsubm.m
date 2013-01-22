
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

function runSurfacePsubm(setupfile,outputDir,varargin)
%
% run surface locally principal submanifold
%

addpath(genpath('code'));

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
    exist('p') && ...
    exist('data2d');
assert(allDefined,'incorrect setup: missing variables...');
global prepend
prepend = fullfile(outputDir,[name '-']);
printFig = @(name) print([prepend name '.ps'],'-dps');

% number of samples
N = size(data2d,2);

format short e;

% setup manifold: quadratic surface
F = @(x) sum(C'.*(x.^2))-1.0;
DF = @(x) 2*C.*x';
D2F = @(x) 2.0*diag(C);
ambient = embeddedManifold(m,n,F,DF,D2F,tol);

% tangent space
if ~exist('B','var')
	B = ambient.orthFrame(p);
else
	assert(ambient.isTangent(B,p));
end

% data
figure(2)
clf
plot(data2d(1,:),data2d(2,:),'r*','MarkerSize',10);
axis equal

dataM = [];
for i = 1:N
    x2 = data2d(:,i);
    [x v] = ambient.Exp(p,B*x2);
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

% locally principal submanifold
psubm = psubManifold(dataM,1,ambient,tol);

1;

end

