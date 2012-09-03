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

setupfile = 'examples/setupfiles/qpgatest.m';
outputDir = 'tmp/output';
tmpDir = 'tmp/';
nrProcesses = 1;
exitfile = 'tmp/exitfile';

% if running in parallel, set master = false for slave processes
assert(nrProcesses == 1);
master = true;

% loop over varying curvature contained in list Cs
if exist(exitfile,'file')
    delete exitfile;
end
if master % this branch run by master process
    [Vapprox Vexact sapprox sfletcher sexact angularDiff] = runQuadraticPGA(setupfile,outputDir,tmpDir,int2str(nrProcesses),exitfile);
else % slave processes
    runQuadraticPGA(setupfile,outputDir,tmpDir,int2str(nrProcesses),exitfile);
end
