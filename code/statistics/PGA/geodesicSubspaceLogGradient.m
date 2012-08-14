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

function [J g Bx] = geodesicSubspaceLogJacobian(x,mu,y,w,Logxy,V,v,Vvp,tol,mode,DF,Exp,Log,dExp,d2Exp)
%
% Compute Jacobian and gradient of either Log_\mu(\pi_{S_v}(x)) or Log_x(\pi_{S_v}(x))
% with S_v=Exp_\mu V_v and V_v=\Span(V,v)
%
% (V,v,Vvp) must be orthonormal and constitute a full basis for T_\mu M
%
% Output in Vvp basis for gradient, RR^m,Vvp for variance Jacobian, Bx,Vvp
% for reconstruction error Jacobian
%
% mode is either 'V' for variance or 'R' for reconstruction error
%
% Algorithm as in PGA paper
%

epsilon = 10e-5; % shouldn't be hardcoded

m = size(x,1);
dimM = size(V,2)+1+size(Vvp,2);
k = size(V,2);
Vv = [V v];
VvVvp = [Vv Vvp]; % full basis for T_mu M as used in the paper
assert(all(size(VvVvp) == [m dimM]));
assert(isOrthonormal(VvVvp));

if dimM ~= m
    By = null(DF(y)); % basis for T_y M
    Bx = null(DF(x)); % basis for T_x M
else
    % Euclidean
    By = eye(m);
    Bx = eye(m);
end

% gradient of R in Vvp
[xx vv solExpxLogxy] = Exp(x,Logxy,tol);
dLogxyExpx = By'*dExp(solExpxLogxy,Bx,tol); % By, Bx
%dyLogx = inv(dLogxyExpx); % Bx, By - we use dLogxyExpx\ instead
[xx vv solExpmuw] = Exp(mu,w,tol);
[A soldwExpmu] = dExp(solExpmuw,VvVvp,tol);
dwExpmu = By'*A; % By, VvVvp
%gR = 2*(dyLogx*dwExpmu*VvVvp'*Vvp)'*Bx'*Logxy;
gR = 2*(dLogxyExpx\(dwExpmu*VvVvp'*Vvp))'*Bx'*Logxy;

% Hessian in (VvVvp Vv)
HwR = [];
D2Expmu = d2Exp(solExpmuw,VvVvp,Vv,soldwExpmu,[],tol); % RR^m, (VvVvp Vv) - recalculation of Vv-part can be avoided
%Vvx = dyLogx*dwExpmu*VvVvp'*Vv; % Bx
Vvx = dLogxyExpx\(dwExpmu*VvVvp'*Vv); % Bx
D2Expx = d2Exp(solExpxLogxy,Bx,Bx*Vvx,[],[],tol); % RR^m, (Bx Vvx)
assert(size(D2Expmu,3) == k+1);
assert(size(D2Expx,3) == k+1);
for i = 1:k+1
    dsdExpmu = By'*D2Expmu(:,:,i); % By, (VvVvp Vv)
    dsdExpx = By'*D2Expx(:,:,i); % By, (Bx Vvx)
    assert(all(size(dsdExpmu) == [dimM dimM]));
    assert(all(size(dsdExpmu) == [dimM dimM]));
    
    % in B
    %HwRi = -2*(dyLogx*dsdExpx*dyLogx*dwExpmu)'*Bx'*Logxy ...
    %    +2*(dyLogx*dsdExpmu)'*Bx'*Logxy ...
    %    +2*(dyLogx*dwExpmu)'*Vvx(:,i);
    HwRi = ... 
        -2*(dLogxyExpx\(dsdExpx*(dLogxyExpx\dwExpmu)))'*Bx'*Logxy ...
        +2*(dLogxyExpx\dsdExpmu)'*Bx'*Logxy ...
        +2*(dLogxyExpx\dwExpmu)'*Vvx(:,i);    
    
    HwR(:,i) = HwRi;
end
assert(all(size(HwR) == [dimM k+1]));

% everything here in (VvVvp Vv)
HwR = HwR'; % transpose
HwRVv = HwR(:,1:k+1);
assert(size(HwRVv,1) == size(HwRVv,2));
%invHwRVv = inv(HwRVv); - we use HwRVv\ instead
%invHwRVvend = invHwRVv(:,end);
invHwRVvend = HwRVv\[zeros(k,1); 1];
BB = HwR(:,k+2:end);

% differential of \pi_{S_v} in By, Vvp
barv = VvVvp'*Vv*invHwRVvend; % VvVvp
wk1 = VvVvp(:,k+1)'*w;
%barE = wk1*[zeros(dimM,1) [- invHwRVv*BB; eye(dimM-(k+1))]];
barE = wk1*[- HwRVv\BB; eye(dimM-(k+1))];
dvpi = dwExpmu*(-barv*gR'+barE);

if mode == 'V'
    % Jacobian in B, Vvp
    %dyLogmu = inv(dwExpmu); % VvVvp, By - we use dwExpmu\ instead
    %J = VvVvp*dyLogmu*dvpi; % RR^m, Vvp
    J = VvVvp*(dwExpmu\dvpi); % RR^m, Vvp
    g = 2*J'*w;
else if mode == 'R'
    % Jacobian in Bx, Vvp
    %J = dyLogx*dvpi; % Bx, Vvp
    J = dLogxyExpx\dvpi; % Bx, Vvp
    g = 2*J'*(Bx'*Logxy);
    else
        assert(false);
    end
end

end