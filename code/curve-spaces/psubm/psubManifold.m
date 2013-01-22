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

function manifold = psubManifold(xi,dim,ambient,tol)
%
% Wrapper function exporting Exp and Log maps
% for locally principal submanifold of embedding manifold 'ambient'
%

% debug
if ~isempty(whos('global','debug'))
    global debug;
else
    debug = false;
end

optimoptions = [];
if debug
    optimoptions.verbose = true;
end
optimoptions.tol = tol;
nonEucGradOptimizer = getNonEucGradOptimizer(optimoptions);

m = ambient.m;

manifold.dim = dim;
assert(manifold.dim >= 1 && manifold.dim < ambient.dim);
assert(manifold.dim == ambient.dim-1); % for now

manifold.type = 'psumb';
manifold.ambient = ambient;
assert(strcmp(ambient.type,'embedded')); % supported for now
manifold.xi = xi;

    % Exponential map
    function [x,v,sol] = Exp(x0,v0,varargin)
        ltol = tol;
        tspan = [0 1];
        if size(varargin,2) >= 1
             ltol = varargin{1};
        end
        if size(varargin,2) >= 2
             tspan = [0 varargin{2}];
        end    
        sol = intExp(x0,v0,tspan,manifold,ltol);
        [x v] = getExp(sol,manifold,tspan(2));
    end

    %
    % find point on psumb 'reasonably' close to x
    function mp = toManifold(p)
        
        % optimization curry
        function [R,grad] = f(X)
            mu = X(1:m);
            n = X((m+1):end);

            Bmu = ambient.orthFrame(mu);
            V = Bmu*null((Bmu'*n)');            
            [R,gradmu,gradn] = approxResidualF(xi,mu,V,n,Bmu,ambient);
%             fprintf('toManfiold R: %f\n',R); % debug

            grad = [Bmu*gradmu; V*gradn];
        end
        
        % project back to manifold
        function px = stepF(X)
            mu = X(1:m);            
            n = X((m+1):end);
            
            pmu = ambient.toManifold(mu);
            pn = n/norm(n);
            
            px = [pmu; pn];
        end
        
        % initial guess
        fletcherV = approxPGA(xi,p,ambient);        
        initval = [p; fletcherV(:,dim+1)];
%         initval = [p; [1 1 0]'/sqrt(2)]; % debug
        val = nonEucGradOptimizer(initval,@f,@stepF);
        mp = val(1:m);
    end

    function res = isTangent(V,p)
        res = norm(frame(p)'*V) < epsilon();
    end

    % get frame for tangent space
    % is is normal to frame
    function [fr n] = frame(p)
        
        Bp = ambient.orthFrame(p);

        % optimization curry
        function [R,grad] = f(n)
            
            V = Bp*null((Bp'*n)');            
            [R,gradmu,gradn] = approxResidualF(xi,p,V,n,Bp,ambient);

            grad = V*gradn;
        end
        
        % normalize n
        function pn = stepF(n)
            pn = n/norm(n);
        end

        % initial guess
        fletcherV = approxPGA(xi,p,ambient);        
        initval = fletcherV(:,dim+1);
%         initval = [1 1 0]'/sqrt(2); % debug
        n = nonEucGradOptimizer(initval,@f,@stepF);
        fr = Bp*null((Bp'*n)');        
    end

    function frame = orthFrame(p)              
        frame = orth(frame(p));
    end

    function res = isPoint(p)
        res = ambient.isPoint(p) ...
            && norm(toManifold(p)-p) < epsilon();
    end

manifold.isTangent = @isTangent;
manifold.isPoint = @isPoint;
manifold.frame = @frame;
manifold.orthFrame = @orthFrame;
manifold.toManifold = @toManifold;
manifold.Exp = @Exp;

end
