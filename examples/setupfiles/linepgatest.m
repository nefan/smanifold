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

% setup file for runSurfacePGA

name = 'linepgatest';

m = 3;
n = 1;

c = sscanf(varargin{1},'%e'); % get coefficient from parameters
C = [c(1) 1 1]
p = [0; 0; 1];
B = eye(m,m-n);

tol = 1e-4;

nDraws = 4;
dataRaw = [-1:1/((nDraws/2-1)/2):1; repmat(0,1,nDraws/2)]*pi*2/4;
shiftM = eye(2);
v = 1*pi/32+0*pi/2;
rotM = [cos(v) -sin(v); sin(v) cos(v)];
data2d = rotM * shiftM * dataRaw;
v = -1*pi/32+1*pi/2;
rotM = [cos(v) -sin(v); sin(v) cos(v)];
data2d = [data2d rotM * shiftM * dataRaw];

% exact PGA stuff
global mode;
mode = 'R';

global debug;
debug = false;

global collectfile;
collectfile = [name '-collect.mat'];
