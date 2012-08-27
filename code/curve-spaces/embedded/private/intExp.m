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

function sol = intExp(x0,v0,tspan,m,n,F,DF,D2F,tol)
%
% Compute geodesic.
%
% Curve will start at x0 with initial velocity v0 and
% will be computed in the interval [0 1].
% Submanifold is defined as the preimage F^-1(0) of
% the submersion F
%
% Follows Dedieu/Nowicki.
%

epsilon = 10e-5; % shouldn't be hardcoded

dimM = m-n;

function dy = G(t,y)
    x = y(1:m); % curve
    p = y(m+1:2*m); % velocity
    
    if dimM == m % Eucledian
        dy = [p; zeros(m,1)];
        return;
    end    
    
    DFx = DF(x);
    D2Fx = D2F(x);
    GInvDFx = GInv(DFx);
    DFxtrans = DFx';
    GInvDFxtrans = GInv(DFxtrans);
    
    mu = -GInvDFxtrans*p;

    % x
    dx = p-GInvDFx*DFx*p;
    
    % p
    D2Fxdx = reshape(D2Fx*dx,n,m);
    dp = -D2Fxdx'*mu;
    
    dy = [ dx; dp];
end

if dimM ~= m
    if norm(F(x0)) > epsilon
        norm(F(x0))
        error('Initial point not on manifold.')
    end
    if norm(DF(x0)*v0) > epsilon
        DF(x0)*v0    
        error('Initial velocity not in tangent space.')
    end
end

% initial values
% x0
p0 = v0;
y0 = [ x0; p0];

% integrate
options = odeset('RelTol',min(1e-6,tol*1e-1),'AbsTol',min(1e-6,tol*1e-1));
sol = ode45(@G,tspan,y0,options);

end