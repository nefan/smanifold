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

function sol = intDExp(dx,dw,tspan,solExp,m,n,F,DF,D2F,quadratic,tol,varargin)
%
% Compute DExp with respect to dx (base point) and dw (tangent vector) in
% Exp_x w
%

epsilon = 10e-5; % shouldn't be hardcoded

x0 = getExp(solExp,DF,tspan(2));
dimM = m-n;
rankx = size(dx,2);
rankw = size(dw,2);
rank = rankx+rankw;

backwards = false; % integrate backwards
if size(varargin,2) >= 1
    backwards = varargin{1};
end

function dy = G(t,y)
    dex = reshape(y(1:rank*m),m,rank); % z
    dep = reshape(y(numel(dex)+1:numel(dex)+rank*m),m,rank); % y
    
    if dimM == m % Eucledian
        dy = [y(m*rank+1:end,1); zeros(m*rank,1)];
        return;
    end        
        
    [x v p] = getExp(solExp,DF,t);        
    
    DFx = DF(x);
    D2Fx = D2F(x);
    GInvDFx = GInv(DFx);
    DFxtrans = DFx';
    GInvDFxtrans = GInvDFx';
    
    mu = -GInvDFxtrans*p;

    % x
    dx = p-GInvDFx*DFx*p;
    D2Fxdx = D2Fx*dx;
    D2Fxdx = reshape(D2Fxdx,n,m);     
    
    % p
    %D2Fxdx = reshape(D2Fx*dx,n,m);
    %dp = -D2Fxdx'*mu; 
    
    % dexp            
    D2Fxdex = D2Fx*dex;
    deDFx = reshape(D2Fxdex,n,m,rank);
    deGInvDFx = dGInv(DFx,GInvDFx,deDFx);
    deGInvDFxDFx = mtimesx(deGInvDFx,DFx)+mtimesx(GInvDFx,deDFx);
    dedx = dep-reshape(mtimesx(deGInvDFxDFx,p),m,rank)-mtimesx(GInvDFx*DFx,dep);

    demu = -reshape(mtimesx(deGInvDFx,'T',p),n,rank)-mtimesx(GInvDFxtrans,dep);

    D2Fxdedx = reshape(D2Fx*dedx,n,m,rank);            
    assert(quadratic); % non-quadratic case not implemented           
    dedp = -D2Fxdx'*demu-reshape(mtimesx(D2Fxdedx,'T',mu),m,rank); % add 3rd derivative component for non-quadratic F here
            
    ddex = reshape(dedx,m*rank,1);
    ddep = reshape(dedp,m*rank,1);     
    
    dy = [ddex; ddep];
    
    if false % true for debug
        dytest = Gtest(t,y);
        err = norm(dy-dytest);
        assert(err < 1e-8);
    end
end

function dy = Gtest(t,y)
    dex = reshape(y(1:rank*m),m,rank); % z
    dep = reshape(y(numel(dex)+1:numel(dex)+rank*m),m,rank); % y
    
    if dimM == m % Eucledian
        dy = [y(m*rank+1:end,1); zeros(m*rank,1)];
        return;
    end        
        
    [x v p] = getExp(solExp,DF,t);        
    
    DFx = DF(x);
    D2Fx = D2F(x);
    GInvDFx = GInv(DFx);
    DFxtrans = DFx';
    GInvDFxtrans = GInv(DFxtrans);
    
    mu = -GInv(DFx')*p;

    % x
    dx = p-GInvDFx*DFx*p;
    D2Fxdx = D2Fx*dx;
    D2Fxdx = reshape(D2Fxdx,n,m);     
    
    % p
    %D2Fxdx = reshape(D2Fx*dx,n,m);
    %dp = -D2Fxdx'*mu;
    
    % dexp
    ddex = [];
    ddep = [];    
    for i = 1:size(dex,2)
        deix = dex(:,i); % z
        deip = dep(:,i); % y
            
        D2Fxdeix = D2Fx*deix;
        deiDFx = reshape(D2Fxdeix,n,m);
        deiGInvDFx = dGInv(DFx,GInvDFx,deiDFx);
        deiGInvDFxDFx = deiGInvDFx*DFx+GInvDFx*deiDFx;
        deidx = deip-deiGInvDFxDFx*p-GInvDFx*DFx*deip;
        
        deimu = -dGInv(DFxtrans,GInvDFxtrans,deiDFx')*p-GInvDFxtrans*deip;
                   
        D2Fxdeidx = reshape(D2Fx*deidx,n,m);            
        assert(quadratic); % non-quadratic case not implemented           
        deidp = -D2Fxdx'*deimu-D2Fxdeidx'*mu; % add 3rd derivative component for non-quadratic F here
            
        % save
        ddex(:,i) = deidx;
        ddep(:,i) = deidp;
    end
    ddex = reshape(ddex,m*rank,1);
    ddep = reshape(ddep,m*rank,1);    
    
    dy = [ddex; ddep];
end

function dy = bG(t,y)
    dy = -G(1-t,y);
end

if size(varargin,2) >= 2
    y0 = varargin{2};
else
    % dx
    dex0 = reshape(dx,m*rankx,1);    
    dep0 = reshape(zeros(m,rankx),m*rankx,1);
    
    y0 = [dex0; dep0];
    
    % dw
    dex0 = reshape(zeros(m,rankw),m*rankw,1);
    dep0 = reshape(dw,m*rankw,1);    
    
    y0 = [y0; [dex0; dep0]];
    
    assert(size(y0,1) == 2*m*rank);
end

% integrate
options = odeset('RelTol',min(1e-6,tol*1e-1),'AbsTol',min(1e-6,tol*1e-1));
if ~backwards
    sol = ode45(@G,tspan,y0,options);
else
    sol = ode45(@bG,tspan,y0,options);
end

end