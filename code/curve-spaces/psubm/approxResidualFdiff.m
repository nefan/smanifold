%  sambient, Copyright (C) 2009-2012, Stefan Sommer (sommer@diku.dk)
%  https://github.com/nefan/sambient.git
% 
%  This file is part of sambient.
% 
%  sambient is free software: you can redistribute it and/or modify
%  it under the terms of the GNU General Public License as published by
%  the Free Software Foundation, either version 3 of the License, or
%  (at your option) any later version.
% 
%  sambient is distributed in the hope that it will be useful,
%  but WITHOUT ANY WARRANTY; without even the implied warranty of
%  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%  GNU General Public License for more details.
% 
%  You should have received a copy of the GNU General Public License
%  along with sambient.  If not, see <http://www.gnu.org/licenses/>.
%  

function dn = approxResidualFdiff(xi,mu,dx,n,ambient)
%
% derivative in W directions of approxResidualF
%

% debug
if ~isempty(whos('global','debug'))
    global debug;
else
    debug = false;
end

N = size(xi,2);
m = size(xi,1);

% debug
dn = zeros(m,1);
return;

Vv = [V v];
assert(isOrthonormal(Vv));
dimV = size(V,2);
dimVv = size(Vv,2);
assert(dimVv <= ambient.dim);

Logmuys = zeros(m,N);
ys = zeros(m,N);
Rs = [];
gradmus = zeros(ambient.dim,N);
gradVs = zeros(dimV,N);
parfor j = 1:N
    x = xi(:,j);
    Logmuy = V*V'*ambient.Log(mu,x);
    Logmuys(:,j) = Logmuy;
    [y zz solExpmuLogmuy] = ambient.Exp(mu,Logmuy);
    ys(:,j) = y;
    Logxy = ambient.Log(x,y);
    Rs(j) = sum(Logxy.^2);
    
    % V- and v gradient of R in V basis
    By = ambient.orthFrame(y); % basis for T_y M
    Bx = ambient.orthFrame(x); % basis for T_x M
    
    [xx,vv,solExpxLogxy] = ambient.Exp(x,Logxy);
    dLogxyExpx = By'*ambient.DExp(solExpxLogxy,[],Bx); % By, Bx
    %dyLogx = inv(dLogxyExpx); % Bx, By - we use dLogxyExpx\ instead
    [xx,vv,solExpmuLogmuy] = ambient.Exp(mu,Logmuy);
    [Av,soldvmuLogmuyExp] = ambient.DExp(solExpmuLogmuy,[],v);
    dvmuLogmuyExp = By'*Av; % By, v
    %gradv = 2*(dyLogx*dvmuLogmuyExp)'*Bx'*Logxy;
    gradv = 2*(dLogxyExpx\dvmuLogmuyExp)'*Bx'*Logxy;
    gradVs(:,j) = -gradv*V'*Logmuy;
    
    % mu-gradient of R in B basis
    [Amu,soldmumuLogmuyExp] = ambient.DExp(solExpmuLogmuy,B,[]);
    dmumuLogmuyExp = By'*Amu; % By, B
    %gradV = 2*(dyLogx*dmumuLogmuyExp)'*Bx'*Logxy;
    gradmu = 2*(dLogxyExpx\dmumuLogmuyExp)'*Bx'*Logxy;
    gradmus(:,j) = gradmu;
end

R = sum(Rs);
gradV = sum(gradVs,2);
gradmu = sum(gradmus,2);
