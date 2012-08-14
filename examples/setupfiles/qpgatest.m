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

name = 'qpgatest';

m = 5;
n = 1;

C = [1 -2 1 -2 1]
tol = 1e-4;

p = [0; 0; 0; 0; 1];

% generate data
%stream = RandStream.getDefaultStream;
%stream.reset;
nDraws = 32
randn(m-1,149+100*70); % vary data a bit
dataTM = randn(nDraws,m-n)';
shiftM = diag([2 1 2/3 1/3])
%v = 16*pi/32+pi/2;
%rotM = [cos(v) -sin(v); sin(v) cos(v)];
rotM = eye(m-n)
B = eye(m,m-n)

dataTM = rotM * shiftM * dataTM;

% S = load('/others/sommer/vertebra/vertebra1/output/quadraticPGA/mgnt03.1677644.0/qpgatest-data.mat');
% dataTM = S.dataTM;


% dataRaw = [1; 0]*pi/3;
% shiftM = eye(2);
% v = 1*pi/32+0*pi/2;
% rotM = [cos(v) -sin(v); sin(v) cos(v)];
% data2d = rotM * shiftM * dataRaw;
% v = -1*pi/32+1*pi/2;
% rotM = [cos(v) -sin(v); sin(v) cos(v)];
% data2d = [data2d rotM * shiftM * dataRaw];
% 
% dataTM = [data2d; zeros(2,2)]
