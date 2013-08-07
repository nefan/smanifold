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

function dGInvA = dGInv(A,GInvA,dA)
%
% derivative of generalized inverse, confer Decell '72
%

assert(size(A,1) == size(dA,1));
assert(size(A,2) == size(dA,2));
    
dAtrans = permute(dA,[2 1 3]);
GInvAtrans = GInvA';
% T = dAtrans*GInvAtrans*GInvA+GInvA*GInvAtrans*dAtrans;
% dGInvA = -GInvA*dA*GInvA+T-GInvA*A*T*A*GInvA;
T = mtimesx(dAtrans,GInvAtrans*GInvA)+mtimesx(GInvA*GInvAtrans,dAtrans);
dGInvA = -mtimesx(mtimesx(GInvA,dA),GInvA)+T-mtimesx(mtimesx(GInvA*A,T),A*GInvA);
