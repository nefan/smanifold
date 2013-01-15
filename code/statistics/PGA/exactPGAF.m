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

function [fval ys ws Logxys rs linfval Rs] = exactPGAF(data,v,Vk,mode,Fproj,projTol,debug)
%
% evaluate projections
%

if debug
    global timeProj;
end

N = size(data,2);

ws = [];
ys = [];
Logxys = [];
rs = [];     
variance = 0;
R = 0;
Rs = [];
linVar = 0;
linR = 0;
% debug
if debug
    tic
end
parfor j = 1:N
    x = data(:,j);
    res = Fproj(x,v,Vk,projTol);

    ys(:,j) = res{1};
    ws(:,j) = res{2};        
    Logxys(:,j) = res{4}; 
            
    variance = variance + sum(ws(:,j).^2);
    R = R + res{3};
    Rs(j) = res{3};
    linVar = linVar + res{5};    
    linR = linR + res{6};    
end 
% debug
if debug
    timeProj = timeProj + toc;
end        

linfval = 0;
if mode == 'V'
    fval = variance/N;
    rs = ws; % B
    linfval = linVar/N;
else if mode == 'R'
        fval = R/N;
        linfval = linR/N;
    else 
        assert(false);
    end
end  
