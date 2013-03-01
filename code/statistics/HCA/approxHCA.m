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

function [V,s,coords] = approxHCA(data,mu,nr,B,manifold)
%
% Compute linearized HCA (Horizontal Component Analysis) of
% the samples xi in T_mean M
%
% V and D will be the eigenvalues and eigenvectors resp.
% of a decomposition of the tangent space T_mean R^m
%
%

N = size(data,2); % number of points

datak = zeros(nr,N);
V = [];
Vp = eye(manifold.dim);
s = [];
for k=1:nr
    u = zeros(manifold.dim-k+1,N);
    parfor j = 1:N
        if k > 1
            datakj = B*V*datak(1:(k-1),j);
        else
            datakj = B*zeros(manifold.dim,1);
        end
        [xj,vj,solExp] = manifold.Exp(mu,datakj);
        Bj = manifold.Pt(solExp,B*Vp);
        u(:,j) = Bj'*manifold.Log(xj,data(:,j));
    end
    S = zeros(manifold.dim-k+1,manifold.dim-k+1);
    for j = 1:N
        S = S + u(:,j)*u(:,j)';
    end
    S = 1/N*S;

    [Vk,Dk] = eig(S);
    Vk(:,end:-1:1) = Vk;
    sk = diag(Dk)';
    sk(1,end:-1:1) = sk;
    
    v = Vk(:,1);
    datak(k,:) = v'*u;
    
    V = [V Vp*v];
    Vp = null(V');
    s = [s sk(1)];
end

V = B*V; % back to coordinates or RR^m
s = cumsum(s);
coords = datak;