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

name = 'uniformhca';

m = 3;
n = 1;

C = [1 1 1]
tol = 1e-4;

p = [1; 0; 0];

% data
% generate data
rng(119); % seed
nDraws = 2^9
datauniform = zeros(m-n,nDraws);
datauniform(1,:) = 2*pi*(rand(nDraws,1)'-repmat(0.5,1,nDraws));
shiftM = diag([0.35])
datauniform(2:end,:) = shiftM*randn(nDraws,m-n-1)';
%v = 16*pi/32+pi/2;
%rotM = [cos(v) -sin(v); sin(v) cos(v)];
rotM = eye(m-n)
B = [0 1 0; 0 0 1]'
