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

function [Vapprox Vexact sapprox sfletcher sexact angularDiff] = runQuadraticPGA(setupfile,outputDir,tmpDir,varargin)
%
% run pga on quadratic manifold
%

addpath(genpath('code'));

epsilon = 10e-5;

% show figures?
visible = false;
if not(visible)
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
    exist('dataTM');
assert(allDefined,'incorrect setup: missing variables...');
global prepend
prepend = fullfile(outputDir,[name '-']);
printFig = @(name) print([prepend name '.ps'],'-dps');

% number of samples
N = size(dataTM,2);

format short e;

% setup manifold: quadratic surface
dimM = m-n;
F = @(x) sum(C'.*(x.^2))-1.0;
DF = @(x) 2*C.*x';
D2F = @(x) 2.0*diag(C);
[distM Exp Log LogInit dExp d2Exp] = ExpLogMaps(m,n,F,DF,D2F,tol);

% tangent space
if ~exist('B','var')
	B = null(DF(p));
else
	assert(norm(DF(p)*B) < epsilon);
end

% curvature
curvature(p,B(:,1),B(:,2),m,Exp,dExp,tol);
curvature(p,B(:,1),B(:,3),m,Exp,dExp,tol);
curvature(p,B(:,1),B(:,4),m,Exp,dExp,tol);
curvature(p,B(:,2),B(:,3),m,Exp,dExp,tol);
curvature(p,B(:,2),B(:,4),m,Exp,dExp,tol);
curvature(p,B(:,3),B(:,4),m,Exp,dExp,tol);

dataM = [];
for i = 1:N
    x2 = dataTM(:,i);
    [x v] = Exp(p,B*x2);
    y3 = x;
        
    dataM(:,i) = y3;
end

% approximated PGA
[Vapprox sapprox] = approxPGA(dataM,p,B,Log);

% real PGA
[Vexact sexact sfletcher] = exactPGA(dataM,p,dimM,B,tol,DF,Log,Exp,dExp,d2Exp,Vapprox);

Vapprox
Vexact
sapprox
sfletcher
sexact
fvalDiff = 100*(sexact-sfletcher)./sfletcher % percent
angularDiff = [];
for i = 1:size(Vexact,2)
    angularDiff(i) = acos(abs(dot(Vexact(:,i),Vapprox(:,i))))*360/(2*pi);
end
angularDiff % degress
save([prepend 'data.mat'],'dataTM','dataM','Vapprox','Vexact','sapprox','sfletcher','sexact','fvalDiff','angularDiff');

end
