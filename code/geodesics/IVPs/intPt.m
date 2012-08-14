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

function sol = intPt(B0,solExp,m,n,F,DF,D2F,tol,varargin)
%
% Compute parallel transport.
%

epsilon = 10e-5; % shouldn't be hardcoded


backwards = false; % integrate backwards
if size(varargin,2) >= 1
    backwards = varargin{1};
end

x0 = getExp(solExp,DF);
dimM = m-n;
rank = size(B0,2);

function dy = G(t,y)
    Pt = reshape(y,m,rank);
    
    [x v p] = getExp(solExp,DF,t);            
    
    DFx = DF(x);
    D2Fx = D2F(x);
    GInvDFx = GInv(DFx);    

    % x
    dx = p-GInvDFx*DFx*p;
    
    % p
    D2Fxdx = reshape(D2Fx*dx,n,m);
    %dp = -D2Fxdx'*mu;
    
    % parallel transport
    dPt = -GInv(DFx)*D2Fxdx*Pt;
    dPt = reshape(dPt,m*rank,1);
    
    dy = dPt;
end

function dy = bG(t,y)
    dy = -G(1-t,y);
end

% initial values
Pt0 = reshape(B0,m*rank,1);

y0 = Pt0;
% integrate
options = odeset('RelTol',min(1e-6,tol*1e-1),'AbsTol',min(1e-6,tol*1e-1));
if ~backwards
    sol = ode45(@G,[0 1],y0,options);
else
    sol = ode45(@bG,[0 1],y0,options);
end

end