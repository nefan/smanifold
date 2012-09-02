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

function [dir] = shoot(p0,p1,m,n,F,DF,D2F,Exp,DExp,tol,varargin)
%
% Compute geodesic from p0 to p1 on manifold defined by
% the premimage F^-1(0) of the submersion F. The
% geodesic is produced by repeated shooting.
%

epsilon = 10e-5; % for debugging

if norm(F(p0)) > epsilon
    error('Initial point not on manifold.')
end
if norm(F(p1)) > epsilon
    error('End point not on manifold.')
end
% seems to be a good idea to do an inital projection onto manifold to
% reduce numerical errors
function [y dy] = lf(x)
    y = F(x);
    dy = DF(x);
end
projOptions = optimset('Jacobian','on','Display','off','TolFun',tol*1e-4,'TolX',tol*1e-4);
warning off optim:fsolve:NonSquareSystem
p0 = fsolve(@lf,p0,projOptions);
p1 = fsolve(@lf,p1,projOptions);
if norm(F(p0)) > epsilon
    error('Initial point not on manifold.')
end
if norm(F(p1)) > epsilon
    error('End point not on manifold.')
end

% setup
Np0 = null(DF(p0)); % basis for tangent space at p0
Np1 = null(DF(p1)); % ...
inttol = tol;

% debugging
verbose = 0;

if norm(p0-p1) < tol
    % already done
    dir = zeros(m,1);
    return;
end

% optimization currys
function [y dy] = lmf(x) % standard
    [yy v solExp] = Exp(p0,Np0*x,inttol);
    y = yy - p1;
    dy = DExp(solExp,Np0,inttol);  
end
function [y dy] = lmfb(x) % backwards
    [yy v solExp] = Exp(p1,Np1*x,inttol);
    y = yy - p0;
    dy = DExp(solExp,Np1,inttol);  
end

if verbose < 2
    displayOpts = optimset('Display','off','Diagnostics','off');
else
    displayOpts = optimset('Display','iter','Diagnostics','on');
end

function [dir residual] = doTRR(init,backwards) % Trust-region-reflective
    options = optimset(displayOpts,'Jacobian','on','TolFun',tol,'TolX',tol*1e-4);
    
    if verbose
        if ~backwards
            fprintf('Trust-region-reflective shot...\n');
        else
            fprintf('Trust-region-reflective backwards shot...\n');
        end
    end
    
    if ~backwards
        [dir,resnorm,residual,exitflag,output] = lsqnonlin(@lmf,init,[],[],options);
    else
        [dir,resnorm,residual,exitflag,output] = lsqnonlin(@lmfb,init,[],[],options);
    end
end

function [dir residual] = doLM(init,backwards) % Levenberg-Marquardt
    options = optimset(displayOpts,'Algorithm','levenberg-marquardt','Jacobian','on','TolFun',tol,'TolX',tol*1e-4);
    
    if verbose
        if ~backwards
            fprintf('Levenberg-Marquardt shot...\n');
        else
            fprintf('Levenberg-Marquardt backwards shot...\n');
        end
    end
    
    if ~backwards
        [dir,resnorm,residual,exitflag,output] = lsqnonlin(@lmf,init,[],[],options);
    else
        [dir,resnorm,residual,exitflag,output] = lsqnonlin(@lmfb,init,[],[],options);
    end
end

function [dir residual] = doGN(init,backwards) % Gauss-Newton    
    if verbose
        if ~backwards
            fprintf('Gauss-Newton shot...\n');
        else
            fprintf('Gauss-Newton backwards shot...\n');
        end
    end    
    
    if ~backwards
        [dir residual] = shootalg(p0,p1,Np0,m,n,DF,Exp,DExp,tol,init,'GN',verbose);
    else
        [dir residual] = shootalg(p1,p0,Np1,m,n,DF,Exp,DExp,tol,init,'GN',verbose);
    end
end

function [dir residual] = doGD(init,backwards) % Gradient descent        
    if verbose
        if ~backwards
            fprintf('Gradient descent shot...\n');
        else
            fprintf('Gradient descent backwards shot...\n');
        end
    end
    
    if ~backwards
        [dir residual] = shootalg(p0,p1,Np0,m,n,DF,Exp,DExp,tol,init,'GD',verbose);
    else
        [dir residual] = shootalg(p1,p0,Np1,m,n,DF,Exp,DExp,tol,init,'GD',verbose);
    end
    
end

function dir = invert(dir) % invert backwards shot
    [x v] = Exp(p1,Np1*dir);
    dir = -Np0'*v;
end

function p = shotOk(residual) % determine if shoot is ok
    p = norm(residual) < tol;
end


% initial guess
if size(varargin,2) >= 1
    init = Np0'*varargin{1};
else
    init = Np0'*(p1-p0);
end
initBackwards = Np1'*(p0-p1);

[dir residual] = doLM(init,false);
if ~shotOk(residual)
    [dir residual] = doLM(initBackwards,true);
    dir = invert(dir);
    [dir residual] = doLM(dir,false);
end

if ~shotOk(residual)
    [dir residual] = doTRR(init,false);
end
if ~shotOk(residual)
    [dir residual] = doTRR(initBackwards,true);
    dir = invert(dir);
    [dir residual] = doTRR(dir,false);
end

if ~shotOk(residual)
    [dir residual] = doGN(init,false);
end
if ~shotOk(residual)
    [dir residual] = doGN(initBackwards,true);
    dir = invert(dir);
    [dir residual] = doGN(dir,false);
end

if ~shotOk(residual)
    [dir residual] = doGD(init,false);
end
if ~shotOk(residual)
    [dir residual] = doGD(initBackwards,true);
    dir = invert(dir);
    [dir residual] = doTRR(dir,false);
end


if ~shotOk(residual) % debug
    fname = ['shooting-bailing-out-' strrep(char(java.util.UUID.randomUUID),'-','_') '.mat'];
    save(fname,'p0','p1','tol');  
    assert(false);
end

% everything ok
dir = Np0*dir;

end
