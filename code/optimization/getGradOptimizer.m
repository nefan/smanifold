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

function optimizer = getGradOptimizer(varargin)
 % Limitied memory BGFS optimizer

if size(varargin,2) > 0
    userOptions = varargin{1};
else
    userOptions = [];
end

options = optimset();
options.Method = 'lbfgs';
options = optimset(options,'LargeScale','off','HessUpdate','bfgs');
options = optimset(options,'GradObj','on');    
if getOption(userOptions,'derivativeCheck')
    options = optimset(options,'DerivativeCheck','on');
else
    options = optimset(options,'DerivativeCheck','off');
end

% tolerance
if isfield(userOptions,'tolFun')
    options = optimset(options,'TolFun',userOptions.tolFun);
end
if isfield(userOptions,'tolX')
    options = optimset(options,'TolX',userOptions.tolX);
end

if isfield(userOptions,'maxIter')
    options = optimset(options,'MaxIter',userOptions.maxIter);
else
    options = optimset(options,'MaxIter',1000);
end
if isfield(userOptions,'maxFunEvals')
    options = optimset(options,'MaxFunEvals',userOptions.maxFunEvals);
else
    options = optimset(options,'MaxFunEvals',1000);
end


% output options
if getOption(userOptions,'verbose')
    options = optimset(options,'Display','iter-detailed');
    options = optimset(options,'Diagnostics','on');    
    %options = optimset(options,'Display','final-detailed');    
    %options = optimset(options,'OutputFcn',outfun); 
else
    options = optimset(options,'Display','final');
end

    function res = minFuncOptimizer(initialData, F, varargin)
        if size(varargin,2) > 0
            options.OutputFcn = varargin{1}; % iteration visualizer
        end       
        res = minFunc(F,initialData,options);
    end

optimizer = @minFuncOptimizer;

end