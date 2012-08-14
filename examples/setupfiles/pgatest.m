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

% test setup file for runSurfacePGA

name = 'pgatest';

m = 3;
n = 1;

C = [-2 1 1];
tol = 1e-4;

p = [0; 0; 1];
%p = [0; 0; 1/sqrt(1)];

% generate data
%stream = RandStream.getDefaultStream;
%stream.reset;
nDraws = 2
randn(2,157+100*nDraws); % vary data a bit
data2d = randn(2,nDraws);
% data2d(1,1:nDraws/2) = data2d(1,1:nDraws/2)+pi/2;
% data2d(1,nDraws/2+1:end) = data2d(1,nDraws/2+1:end)-pi/2;
%data2d = trnd(2,2,nDraws);
%data2d = [7*pi/8; 0];
%data2d = [3*pi/4; 0];
%data2d = [0:0.05*pi:0.9*pi 1.1*pi:0.05*pi:1.9*pi]; data2d = [data2d; zeros(size(data2d))];
%data2d = [1.0:0.2:1]*0.5*pi; data2d = [data2d; zeros(size(data2d))];
%shiftM = eye(2);
shiftM = [1 0; 0 1/3];
v = 16*pi/32+pi/2;
rotM = [cos(v) -sin(v); sin(v) cos(v)];

data2d = rotM * shiftM * data2d;
