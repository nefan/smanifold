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

function ImplitPlot3D(F,xrange,yrange,zrange,varargin)
%
% Mimic Maples implicitplot3d
%

res = 15;
if size(varargin,2) > 0
    res = varargin{1};
end

[X Y Z] = ndgrid( ...
    linspace(xrange(1),xrange(2),res), ...
    linspace(yrange(1),yrange(2),res), ...
    linspace(zrange(1),zrange(2),res));

for i = 1:res
    for j = 1:res
        for k = 1:res
            V(i,j,k) = F([X(i,j,k); Y(i,j,k); Z(i,j,k)]);
        end
    end
end

p = patch(isosurface(X,Y,Z,V,0.0));

%isonormals(X,Y,Z,V,p);
set(p,'FaceColor','b','EdgeColor','none','FaceAlpha',0.3);
daspect([1 1 1]);
axis equal;
grid on;
camlight
view(24,32);
zoom out;
zoom(2.0);
lighting gouraud;

display(p);
