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

function [V D] = PCA(xi,subtractMean)
%
% Compute the PCA (Principal Component Analysis) of
% the samples xi
%
% V and D will be the eigenvalues and eigenvectors resp.
%
% mean will be computed and subtracted if subtractMean
%

N = size(xi,2); % number of points
m = size(xi,1);

u = [];
if subtractMean
    mean = 1/N*sum(xi,2);


    for i = 1:N
        u(:,i) = xi(:,i)-mean;
    end
else
    u = xi;
end
S = zeros(m,m);
for i = 1:N
    S = S + u(:,i)*u(:,i)';
end
S = 1/N*S;

[V D] = eig(S);