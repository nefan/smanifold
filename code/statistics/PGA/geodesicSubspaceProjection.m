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

function [y w R Logxy w0] = geodesicSubspaceProjection(x,mu,V,w0,tol,Exp,Log,dExp)
%
% Compute projection \pi_S(x) of x onto geodesic subspace
% S=Exp_\mu VV, where VV=\Span V
%
% V must be orthonormal
%
% Algorithm as in PGA paper
%

assert(isOrthonormal(V));

if isempty(w0)
    w0 = V*V'*Log(mu,x,tol); % initial guess (orthogonal projection)
end

function [R J] = lf(w)
    [y v solExp] = Exp(mu,V*w);
    B = dExp(solExp,V);
    assert(size(B,2) == size(V,2));
        
    if isempty(guess)
        R = Log(y,x,tol);
    else
        fprintf('Subspace projection using guess...\n');
        R = Log(y,x,tol,guess);
    end
    %guess = R;
    
    J = -2*B;
end

displayOpts = optimset('Display','off','Diagnostics','off');
%displayOpts = optimset('Display','iter','Diagnostics','on');
options = optimset(displayOpts,'Jacobian','on','TolFun',tol,'TolX',tol);
%warning off optim:fsolve:NonSquareSystem
guess = []; % might speed up things
w = lsqnonlin(@lf,V'*w0,[],[],options);
guess = [];

w = V*w;
[y v solExp] = Exp(mu,w);
Logxy = Log(x,y,tol);
R = sum(Logxy.^2);

end
