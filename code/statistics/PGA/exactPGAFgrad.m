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

function [g Js gs] = exactPGAFgrad(data,v,ys,ws,Logxys,rs,Vk,B,k,mode,Fgrad,gradTol,debug)
%
% evaluate projection gradients
%

global multicoreSettings;
if debug
    global timeGrad;
end

N = size(data,2);

Vvp = null([Vk v]');
BVk = [];
Bv = B*v;
BVvp = [];
if size(Vk,2) > 0
    BVk = B*Vk;
end
if size(Vvp,2) > 0
    BVvp = B*Vvp;
end  
parameterCell = cell(1,N);
for j = 1:N
    x = data(:,j);
    y = ys(:,j);            
    w = ws(:,j);            
    Logxy = Logxys(:,j);
    parameterCell{j} = {x,y,w,Logxy,BVk,Bv,BVvp,gradTol};                        
end        
% debug
if debug
    tic
end   
resultCell = startmulticoremaster(Fgrad, parameterCell, multicoreSettings.conf);
% debug
if debug
    timeGrad = timeGrad + toc;
end   
gs = []; % mainly for visualization
dimM = size(B,2);
g = zeros(dimM,1);        
Js = zeros(dimM*N,dimM-(k+1));   
for j = 1:N
    if mode == 'V'
        Js(1+(j-1)*dimM:j*dimM,:) = B'*resultCell{j}{1}; % B, Vvp
    else if mode == 'R'
            Js(1+(j-1)*dimM:j*dimM,:) = resultCell{j}{1}; % Bx, Vvp
        else
            assert(false);
        end
    end
	gs(:,j) = Vvp*resultCell{j}{2}; % gradient in B
    g = g+gs(:,j); % accumulated gradient
            
    if mode == 'R'
        Bx = resultCell{j}{3};                
        rs(:,j) = Bx'*Logxys(:,j); % Bx
    end
end

Js = Js/N;
g = g/N;