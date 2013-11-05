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

function sol = intD2Exp(solExp,solw,solu,m,n,F,DF,D2F,quadratic,tol)
%
% Compute the derivative D^2Exp.
%
% solW contains solution of DExp with initial values w...
% solU contains solution of DExp with initial values u...
%

epsilon = 10e-5; % shouldn't be hardcoded

x0 = getExp(solExp,DF);
dimM = m-n;
Zw0 = getDExp(solw,m,0);
Zu0 = getDExp(solu,m,0);
rankw = size(Zw0,2);
ranku = size(Zu0,2);

function w = Xi(DM,v)
    DMv = DM*v;
    w = reshape(DMv,n,m);
end

function w = Xi2(DM,v,dv)
    assert(quadratic);
    w = Xi(DM,dv)+zeros(n,m);
end

function w = Upsilon(D2Fx,v)
    w = Xi(D2Fx,v);
end

function w = Psi(D2Fx,v1,v2)    
    w = Upsilon(D2Fx,v2)'*v1;
end

function w = Phi(x,v1,v2,v3)
    assert(quadratic);
    w = zeros(m,1);
end

function w = Rho(x,v1,v2,v3,v4)
    assert(quadratic);
    w = zeros(m,1);
end
    
function dz = dotz(x,dx,p,z,y,D2Fx,DFx,GInvDFx)
    D2Fxz = D2Fx*z;
    dDFx = reshape(D2Fxz,n,m);
    dGInvDFx = dGInv(DFx,GInvDFx,dDFx);
    dGInvDFxDFx = dGInvDFx*DFx+GInvDFx*dDFx;           
    dz = y-dGInvDFxDFx*p-GInvDFx*DFx*y;
end

function dy = G(t,y)  
    assert(numel(y) == 2*m*rankw*ranku);
    Q = reshape(y(1:m*rankw*ranku),m,rankw,ranku);
    R = reshape(y(m*rankw*ranku+1:end),m,rankw,ranku);        
    
    if dimM == m % Eucledian
        dy = zeros(size(y)); % debug, not sure if this is right
        return;
    end   
        
    [x xdot p] = getExp(solExp,DF,t);
    [Zw Yw] = getDExp(solw,m,t);
    [Zu Yu] = getDExp(solu,m,t);   
    
    DFx = DF(x);
    D2Fx = D2F(x);
    GInvDFx = GInv(DFx);
    DFxtrans = DFx';
    GInvDFxtrans = GInv(DFxtrans);
    
    f = DFx;
    g = GInvDFxtrans;    
    
    dQ = [];
    dR = [];        
    for i = 1:rankw
        zw = Zw(:,i);
        yw = Yw(:,i);
        dotzw = dotz(x,xdot,p,zw,yw,D2Fx,DFx,GInvDFx);
        
        for j = 1:ranku
            r = R(:,i,j);
            q = Q(:,i,j);
            
            zu = Zu(:,j);
            yu = Yu(:,j);
            dotzu = dotz(x,xdot,p,zu,yu,D2Fx,DFx,GInvDFx);
            
            % computations go here
            Upsilonu = Upsilon(D2Fx,zu);
            Upsilonw = Upsilon(D2Fx,zw);
            Upsilonwu = Xi2(D2Fx,zw,r);
            Upsilonutrans = Upsilonu';
            Upsilonwtrans = Upsilonw';

            dGInvDFxu = dGInv(DFx,GInvDFx,Upsilonu);
            dGInvDFxw = dGInv(DFx,GInvDFx,Upsilonw);
            dGInvDFxuTrans = dGInv(DFxtrans,GInvDFxtrans,Upsilonutrans);
            dGInvDFxwTrans = dGInv(DFxtrans,GInvDFxtrans,Upsilonwtrans);            
            
            % not optimized
            lastsum = zeros(m,1);
            D2Fxreshaped = reshape(D2Fx,n,m,m);
            for k = 1:n
                D2Fk = reshape(D2Fxreshaped(k,:,:),m,m);
                assert(quadratic); % else use line below, after modifying it to be correct :-)
                %lastsum = lastsum + g(k,1)*p(k,1)*Xi2(D2Fk,zw,zu)*xdot;
            end
            dr = -(dGInvDFxu*f+GInvDFx*Upsilonu)*yw...
                +q-GInvDFx*f*q...
                -(d2GInv(DFx,Upsilonw,Upsilonu,Upsilonwu,GInvDFx)*f+dGInvDFxw*Upsilonu...
                  +dGInvDFxu*Upsilonw+GInvDFx*Upsilonwu)*p...
                -(dGInvDFxw*f+GInvDFx*Upsilonw)*yu;                        
            assert(quadratic);
            dq = Psi(D2Fx,d2GInv(DFxtrans,Upsilonwtrans,Upsilonutrans,Upsilonwu',GInvDFxtrans)*p,xdot)...
                +Psi(D2Fx,dGInvDFxwTrans*yu,xdot)...
                +Psi(D2Fx,dGInvDFxuTrans*yw+g*q,xdot)...
                +Psi(D2Fx,dGInvDFxwTrans*p+g*yw,dotzu)...
                ...%+Phi(x,dGInvDFxwTrans*p+g*yw,xdot,zu)...
                +Psi(D2Fx,dGInvDFxuTrans*p+g*yu,dotzw)...
                +Psi(D2Fx,g*p,dr)...
                ...%+Phi(x,g*p,dotzw,zu)...
                ...%+Phi(x,dGInvDFxuTrans*p+g*yu,xdot,zw)...
                ...%+Phi(x,g*p,dotzu,zw)...
                ...%+Phi(x,g*p,xdot,r)...
                +lastsum;
            
            
            % save
            dQ(:,i,j) = dq;
            dR(:,i,j) = dr;
        end
    end
    dQ = reshape(dQ,m*rankw*ranku,1);        
    dR = reshape(dR,m*rankw*ranku,1);        
    
    dy = [dQ; dR];
end

q0 = reshape(zeros(m,rankw,ranku),m*rankw*ranku,1);
r0 = reshape(zeros(m,rankw,ranku),m*rankw*ranku,1);

y0 = [q0; r0];

% integrate
options = odeset('RelTol',min(1e-6,tol*1e-1),'AbsTol',min(1e-6,tol*1e-1));
sol = ode45(@G,[0 1],y0,options);

end