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

global multicoreSettings;
if debug
    global timeProj;
end

N = size(data,2);

parameterCell = cell(1,N);
for j = 1:N
    x = data(:,j);
            
    parameterCell{j} = {x,v,Vk,projTol};
end        
% debug
if debug
    tic
end
resultCell = startmulticoremaster(Fproj, parameterCell, multicoreSettings.conf);        
% debug
if debug
    timeProj = timeProj + toc;
end        
ws = [];
ys = [];
Logxys = [];
rs = [];     
variance = 0;
R = 0;
Rs = [];
linVar = 0;
linR = 0;
for j = 1:N
    ys(:,j) = resultCell{j}{1};
    ws(:,j) = resultCell{j}{2};        
    Logxys(:,j) = resultCell{j}{4}; 
            
    variance = variance + sum(ws(:,j).^2);
    R = R + resultCell{j}{3};
    Rs(j) = resultCell{j}{3};
    linVar = linVar + resultCell{j}{5};    
    linR = linR + resultCell{j}{6};    
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