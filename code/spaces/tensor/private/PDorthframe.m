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

function f = PDorthframe(p)

p = PDvtoM(p);
n = sqrt(numel(p));

[u,Lambda] = schur(p);
g = u*sqrt(Lambda);

fId = zeros(n^2,n*(n+1)/2);
f = fId;
k = 1;
for i = 1:n
    for j=i:n
        M = zeros(n,n);
        M(i,j) = 1; M(j,i) = 1;
        fId(:,k) = PDMtov(M);
        k = k+1;
    end
end

fOrthId = orth(fId);

for i = 1:size(f,2)
    MId = PDvtoM(fOrthId(:,i));
    M = g*MId*g';
    f(:,i) = PDMtov(M);
end
