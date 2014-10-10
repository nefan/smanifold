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

% test setup file for runQuadraticPGA

name = 'qdiff';

m = 3;
n = 1;

% C = [1 1 1] % sphere
C = [1 0 1] % cylinder
% C = [-2 1 1] % elliptical hyperboloid
tol = 1e-4;

p = [0; 0; 1];

% covariance
covTM = pi*diag([4^2 1.^2])
% covTM = pi*diag([4^2 1.3^2])
% covTM = pi*diag([4^2 1^2])

% xTM = [7/16*pi,3/8*pi]' % sphere and elliptical hyperboloid
xTM = [9/16*pi,2/8*pi]' % cylinder
r = .2

% number of samples
N = 50000 % hyperboloid: 50e3, sphere: 15e3
K = 20 % Brownian steps
