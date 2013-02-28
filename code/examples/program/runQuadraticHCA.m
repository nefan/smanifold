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

function [Vapprox Vexact sapprox sexact] = runQuadraticHCA(setupfile,outputDir,tmpDir,varargin)
%
% run pga on quadratic manifold
%

addpath(genpath('code'));

global epsilon;
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

% % curvature
% curvature(p,B(:,1),B(:,2),manifold,tol);
% curvature(p,B(:,1),B(:,3),manifold,tol);
% curvature(p,B(:,1),B(:,4),manifold,tol);
% curvature(p,B(:,2),B(:,3),manifold,tol);
% curvature(p,B(:,2),B(:,4),manifold,tol);
% curvature(p,B(:,3),B(:,4),manifold,tol);

if ~exist('dataTM2','var')
    % number of samples
    N = size(dataTM,2);
    
    dataM = [];
    for i = 1:N
        x2 = dataTM(:,i);
        [x v] = manifold.Exp(p,B*x2);
        y3 = x;

        dataM(:,i) = y3;
    end
else
    % bimodal
    N1 = size(dataTM1,2);
    N2 = size(dataTM2,2);
    N = N1+N2;
    
    dataM = [];
    
    [p1,v1,solExp] = manifold.Exp(p,B*d1);
    B1 = manifold.Pt(solExp,B);
    for i = 1:N1
        x2 = dataTM1(:,i);
        [x v] = manifold.Exp(p1,B1*x2);
        y3 = x;

        dataM(:,i) = y3;
    end
    
    [p2,v2,solExp] = manifold.Exp(p,B*d2);
    B2 = manifold.Pt(solExp,B);
    for i = 1:N2
        x2 = dataTM2(:,i);
        [x v] = manifold.Exp(p2,B2*x2);
        y3 = x;

        dataM(:,i+N1) = y3;
    end

end

% approximated HCA
[Vapprox sapprox] = approxHCA(dataM,p,2,B,manifold);
 
% approximated PGA
[PGAVapprox PGAsapprox] = approxPGA(dataM,p,B,manifold);

Vapprox, sapprox
PGAVapprox, PGAsapprox

% real HCA
[Vexact sexact] = exactHCA(dataM,p,manifold.dim,B,tol,manifold,Vapprox);

Vexact, sexact

% % real PGA
% [PGAVexact PGAsexact sfletcher] = exactPGA(dataM,p,manifold.dim,B,tol,manifold,PGAVapprox);

end
