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

function [x v p] = getExp(sol,manifold,varargin)
%
% evaluate solution of intExp
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
m = manifold.ambient.m;
x = y(1:m,1);
p = y(m+1:2*m,1);
n = y(2*m+1:end);

% compute velocity
DFx = [manifold.ambient.DF(x); n'];
if numel(DFx) > 0
    GInvDFx = GInv(DFx);
    v = p-GInvDFx*DFx*p;
else
    % Euclidean
    v = p;
end