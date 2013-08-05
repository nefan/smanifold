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

% test setup file for runQuadraticHCA

name = 'qhcatest';

m = 4;
n = 1;

C = [2 -2 1 1]
tol = 1e-4;

p = [0; 0; 0; 1];

% data
d1 = [ 0 0 1]'
d2 = [ 0 0 -1]'
% generate data
rng(119); % seed
nDraws = 2*2^6
dataTM1 = randn(nDraws/2,m-n)';
dataTM2 = randn(nDraws/2,m-n)';
shiftM = diag([0.4 0.6 0.1])
%v = 16*pi/32+pi/2;
%rotM = [cos(v) -sin(v); sin(v) cos(v)];
rotM = eye(m-n)
B = eye(m,m-n)

dataTM1 = rotM * shiftM * dataTM1;
dataTM2 = rotM * shiftM * dataTM2;


% % generate data
% dataTM = [
%            1.00  0.00  1; 
%           -1.00  0.00  1;           
%            1.00  0.00 -1; 
%           -1.00  0.00 -1;
%            0.00  1.00  1; 
%            0.00 -1.00  1; 
%            0.00  1.00 -1; 
%            0.00 -1.00 -1;
%            ]';
% shiftM = diag([0.6 0.4  1])
% %v = 16*pi/32+pi/2;
% %rotM = [cos(v) -sin(v); sin(v) cos(v)];
% rotM = eye(m-n)
% B = eye(m,m-n)
% 
% dataTM = rotM * shiftM * dataTM;
