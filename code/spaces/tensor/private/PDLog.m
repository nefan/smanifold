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

function Logpx = PDLog(p,x)

p = PDvtoM(p);
x = PDvtoM(x);

[u,Lambda] = schur(p);
g = u*sqrt(Lambda);
ginv = Lambda^(-1/2)*u';
y = ginv*x*ginv';
[v,Sigma] = schur(y);
Logpx = g*v*diag(log(diag(Sigma)))*(g*v)';

Logpx = PDMtov(Logpx);
