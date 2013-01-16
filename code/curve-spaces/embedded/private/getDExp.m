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

function [B P] = getDExp(sol,m,varargin)
%
% evaluate solution of intDExp
%

% time
if size(varargin,2) == 0
    assert(sol.x(end) == 1);
    t = 1;
else
    t = varargin{1};
end

% extract solution
y = deval(sol,t);
rank = numel(y)/(2*m);
B = reshape(y(1:numel(y)/2,1),m,rank);
P = reshape(y(numel(y)/2+1:end),m,rank);