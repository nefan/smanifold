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

function [B] = getD2Exp(sol,m,rankw,ranku,varargin)
%
% evaluate solution of intD2Exp
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
assert(numel(y) == 2*m*rankw*ranku);
B = reshape(y(numel(y)/2+1:end),m,rankw,ranku);