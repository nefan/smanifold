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

function ginvA = GInv(A)
%
% GInv Generalized inverse (confer Dedieu/Nowicki Def 2.1)
%   GInv(A) is such that if x = GInv(A)y then
%   x is the orthogonal projection of y onto Im(A)
%   applied to the map B^-1 : Im(A) -> (Ker(A))^T .
%   That is, B is the restriction of A to the orthogonal
%   complement of Ker(A) .
%
% Non-optimized version
%

[m n] = size(A);

[ U S V ] = svd(A);

tol = max(m,n) * max(max(S)) * eps(class(A)); % from orth.m
rankA = sum(sum(S > tol));

SS = diag(S); SS = SS(1:rankA); % diagonal of restriction to U'Im(A)
invSS = diag(1./SS) ;% SS invertible
ginvA = V(:,1:rankA)*invSS*U(:,1:rankA)';