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

% setup file for spherePGASample

name = 'spherepgasample';

m = 3;
n = 1;

C = [1 0.5 1]
p = [0; 0; 1];
B = eye(m,m-n);

tol = 1e-5;

% generate data
nDraws = 50;
rng('default');
randn(2,-80+157*nDraws); % vary data a bit
datauniform = zeros(2,nDraws);
datauniform(1,:) = 2*pi*(rand(1,nDraws)-repmat(0.5,1,nDraws));
% data2d(1,:) = 1.2*randn(1,nDraws);
datauniform(2,:) = 0.1*randn(1,nDraws);
shiftM = eye(2);
v = 0*pi/32+0*pi/2;
rotM = [cos(v) -sin(v); sin(v) cos(v)];
datauniform = rotM * shiftM * datauniform;
B = [1 0 0; 0 1 0]';

% exact PGA stuff
global mode;
mode = 'V';

global debug;
debug = true;

global collectfile;
collectfile = [name '-collect.mat'];
