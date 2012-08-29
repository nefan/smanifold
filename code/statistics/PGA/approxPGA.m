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

function [V s u] = approxPGA(xi,mean,B,manifold)
%
% Compute the PGA (Principal Geodesic Analysis) of
% the samples xi in the orthonormal basis B of T_mean M
%
% V and D will be the eigenvalues and eigenvectors resp.
% of a decomposition of the tangent space T_mean R^m
%
% u contains the data projected to the tangent space of the mean
%

global multicoreSettings;

assert(isOrthonormal(B));

N = size(xi,2); % number of points
dimM = size(B,2);

parameterCell = cell(1,N);
for j = 1:N
    parameterCell{j} = {mean,xi(:,j)};                        
end        
resultCell = startmulticoremaster(manifold.Log, parameterCell, multicoreSettings.conf);
for j = 1:N
    u(:,j) = B'*resultCell{j};
end
S = zeros(dimM,dimM);
for j = 1:N
    S = S + u(:,j)*u(:,j)';
end
S = 1/N*S;

[V D] = eig(S);
V(:,end:-1:1) = V;
V = B*V; % back to RR^m
s = diag(D)';
s(1,end:-1:1) = s;
s = cumsum(s);