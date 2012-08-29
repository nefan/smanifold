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

function [Vexact Vapprox sexact sapprox sfletcher angularDiff] = runPGA(xi,mean,nr,tol,manifold,varargin)
% 
% runs PGA computations
%

fprintf('PGAs...\n');

m = size(mean,1);

B = manifold.orthonormalFrame(mean);
if size(varargin,2) >= 1
    B = varargin{1};
end
fprintf('   approximated PGA\n');
tic
[Vapprox sapprox] = approxPGA(xi,mean,B,manifold);
toc
sapprox % debug
c = curvature(mean,Vapprox(:,1),Vapprox(:,2),manifold,tol);
Vapprox = Vapprox(:,1:nr);
sapprox = sapprox(1:nr);
fprintf('   curvature at first approx PGA plane: %e\n',c);
fprintf('   exact PGA\n');
tic
[Vexact sexact sfletcher sfletchernoproj estimate RsDiff] = exactPGA(xi,mean,nr,B,tol,manifold,Vapprox);
toc

sapprox
sfletcher
sexact
diff = sexact-sfletcher
fvalDiff = 100*(sexact-sfletcher)./sfletcher % percent
angularDiff = [];
for i = 1:size(Vexact,2)
    angularDiff(i) = acos(abs(dot(Vexact(:,i),Vapprox(:,i))))*360/(2*pi);
end
angularDiff % degress

global prepend;
save([prepend 'data.mat'],'xi','mean','Vapprox','Vexact','sapprox','sfletcher','sexact','fvalDiff','angularDiff');

% collect results
if ~isempty(whos('global','collectfile'))
    global collectfile;
    if exist(collectfile,'file')
        load(collectfile);
    else
        collect.cs = [];
        collect.estimates = [];
        collect.fvalDiffs = [];
        collect.relfvalDiffs = [];
        collect.projDiffs = [];    
        collect.angularDiffs = [];
        collect.RsDiffs = [];
    end
    collect.cs = [collect.cs; c];
    collect.estimates = [collect.estimates; estimate];
    collect.fvalDiffs = [collect.fvalDiffs; sexact-sfletcher];
    collect.relfvalDiffs = [collect.relfvalDiffs; fvalDiff];
    collect.projDiffs = [collect.projDiffs; sfletchernoproj-sfletcher];
    collect.angularDiffs = [collect.angularDiffs; angularDiff];
    collect.RsDiffs = [collect.RsDiffs; RsDiff];
    save(collectfile,'collect');
end
