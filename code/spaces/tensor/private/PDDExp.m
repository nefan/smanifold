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

function DExppXW = PDDExp(p,X,W)

p = PDvtoM(p);
X = PDvtoM(X);
W = PDvtoM(W);
n = sqrt(numel(p));

[u,Lambda] = schur(p);
g = u*sqrt(Lambda);
ginv = Lambda^(-1/2)*u';
Y = ginv*X*ginv';
[v,Sigma] = schur(Y);
expSigma = exp(diag(Sigma));

% derivatives
dY = ginv*W*ginv';
dSigma = zeros(n,1);
dv = zeros(n,n);
for i=1:n
    for j=1:n
        % dSigma
        dSigma = dSigma + dY(i,j)*(v(i,:).*v(j,:))';
        
        % dv
        OmegaUij = zeros(n,n);
        for k=1:n
            for l=1:n
                if k~=l
                    OmegaUij(k,l) = v(i,k)*v(j,l)/(expSigma(l)+expSigma(k));
                end
            end
        end
        
        dv = dv + dY(i,j)*(v*OmegaUij);
    end
end



DExppXW = ... 
    (g*dv)*diag(expSigma)*(g*v)' + ...
    (g*v)*diag(dSigma.*expSigma)*(g*v)' + ...
    (g*v)*diag(expSigma)*(g*dv)';
