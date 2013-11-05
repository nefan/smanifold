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

function [xk,Bxk,Vxk] = HCALinToMan(xV,mu,V,B,manifold)
%
% map from linearization to manifold
%
% xV - coordinates
% V - HCA basis vectors
% B - vector/matrix to be transported along path
%
% xk - mapped point
% Bxk - transported B
% Vxk - transported V

D = length(xV);
xk = mu;  
Vxk = V;
Bxk = B;
for d = 1:D % propagate trought horz. components
    [xk,vk,solExp] = manifold.Exp(xk,Vxk(:,d)*xV(d));
    M = manifold.Pt(solExp,[Vxk Bxk]);
    Vxk = M(:,1:size(Vxk,2));
    Bxk = M(:,(size(Vxk,2)+1):end);
end
