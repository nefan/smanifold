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
% example for running PGA on a 2D surface
%

setupfile = 'examples/setupfiles/linepgatest.m';
outputDir = 'tmp/output';
tmpDir = 'tmp/';
nrProcesses = 1;
exitfile = 'tmp/exitfile';

% if running in parallel, set master = false for slave processes
assert(nrProcesses == 1);
master = true;

% loop over varying curvature contained in list Cs
Cs = [-1 -2];
l = [];
for c = Cs
    if exist(exitfile,'file')
        delete exitfile;
    end
    if master % this branch run by master process
        [Vapprox Vexact sapprox sfletcher sexact angularDiff] = runSurfacePGA(setupfile,outputDir,tmpDir,int2str(nrProcesses),exitfile,int2str(c));
        fvalDiff = 100*(sexact-sfletcher)./sfletcher; % percent
        l(:,end+1) = [fvalDiff(1); angularDiff(1)];

        set(0, 'DefaultFigureVisible', 'on')
        figure(10)
        clf
        plot(Cs(1:size(l,2)),l(1,:));
        figure(11)
        clf
        plot(Cs(1:size(l,2)),l(2,:));
    else % slave processes
        runSurfacePGA(setupfile,outputDir,tmpDir,int2str(nrProcesses),exitfile,int2str(c));
    end
end
