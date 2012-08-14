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

function x = isOrthonormal(B)
%
% checks if columns of B constitutes and orthonormal set
%

n = size(B,2);
x = true;

tol = 10^-12;

for i = 1:n
    vi = B(:,i);
    
    if abs(dot(vi,vi)-1) > tol
        x = false;
        break;
    end
    for j = i+1:n
        vj = B(:,j);
        
        if abs(dot(vi,vj)) > tol
            x = false
            break;
        end
    end
    
end