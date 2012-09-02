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

function c = curvature(p,v1,v2,manifold,tol)
%
% measure sectional curvature of manifold
%

V = orth([v1 v2]);
assert(size(V,2) == 2);
v1 = V(:,1);
v2 = V(:,2);

inttol = tol;
mint = 10e-5;

t = 10e-2;
i = 0;
c = inf;
lastc = 0;
ls = [];
ts = [];
while abs(c-lastc) > tol
    [x v solExp] = manifold.Exp(p,v1,inttol*t,t);
    [B solDExp] = manifold.DExp(solExp,v2,inttol*t,t);
    
    lastc = c;
    c = 6/t^3*(t-norm(manifold.getDExp(solDExp,t)));
    %abs(c-lastc) % debug    
    
    ts(end+1) = t;
    ls(end+1) = c;
    
    t = t/2;    
end

% ls % debug
% figure(5)
% plot(ts,ls);
