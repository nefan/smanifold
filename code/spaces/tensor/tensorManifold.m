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

function manifold = tensorManifold(n)
%
% PD(n) tensor manifold
%
% metric at eye(n) transported to entire group
%

% Exponential map
function [x,v,sol] = Exp(x0,v0,varargin)
    if size(varargin,2) >= 2
         tspan = [0 varargin{2}];
    end    
    sol = [];
    x = PDExp(x0,tspan*v0);
    v = PD
end

% DExp
% Derivatives dx and dw for Exp_x w
function [B,sol] = DExp(solExp,dx,dw,varargin)
    ltol = tol;
    tspan = [0 1];    
    if size(varargin,2) >= 1
         ltol = varargin{1};
    end
    if size(varargin,2) >= 2
         tspan = [0 varargin{2}];
    end        
    sol = intDExp(dx,dw,tspan,solExp,m,n,F,DF,D2F,true,ltol);        
    B = getDExp(sol,m,tspan(2));
end

% % DExpSolve
% function [P0] = DExpSolve(solExp,S,varargin)
%     assert(size(S,2) == 1);
%     ltol = tol;
%     if size(varargin,2) >= 1
%          ltol = varargin{1};
%     end
%     ptsol = intPt(S,solExp,m,n,F,DF,D2F,ltol,true);
%     P0 = getPt(ptsol,m);
%     i = 0;
%     while true
%         sol = intDExp(P0,solExp,m,n,F,DF,D2F,true,ltol);
%         [B P] = getDExp(sol,m);
%         res = B-S;
%         norm(res) % debug
%         if norm(res) < tol
%             break;
%         end        
%                         
%         sol = intDExp(res,solExp,m,n,F,DF,D2F,true,ltol,true,[res; P]);
%         [BB PP] = getDExp(sol,m);
%         P0 = P0 + PP;
% 
%         i = i + 1 % debug
%         assert(i <= 40);        
%     end
% end

% D2Exp in U directions
function [B sol] = D2Exp(solExp,B0,U,solw,solu,varargin)
    ltol = tol;
    if size(varargin,2) >= 1
         ltol = varargin{1};
    end
    if isempty(solw)
        solw = intDExp([],B0,[0 1],solExp,m,n,F,DF,D2F,true,ltol);
    end
    if isempty(solu)
        solu = intDExp([],U,[0 1],solExp,m,n,F,DF,D2F,true,ltol);
    end
    sol = intD2Exp(solExp,solw,solu,m,n,F,DF,D2F,true,ltol);
    B = getD2Exp(sol,m,size(B0,2),size(U,2));
end

% parallel transport
function [B sol] = Pt(solExp,B0,varargin)
    ltol = tol;
    if size(varargin,2) >= 1
         ltol = varargin{1};
    end
    sol = intPt(B0,solExp,m,n,F,DF,D2F,ltol);
    B = getPt(sol,m);
end

function d = distM(p1,p2)
    v = lLog(p1,p2);
    d = norm(v);
end

function [v t c] = LogInit(p1,p2,init,varargin)
    assert(false); % not update to new output format
    assert(norm(p1-init(:,1)) < tol);
    assert(norm(p2-init(:,end)) < tol);
    
    ltol = tol;
    if size(varargin,2) >= 1
         ltol = varargin{1};
    end
    
    [t c v] = shorten(init,F,DF,D2F,ltol*0.1); % guessing on tolerance
    
    % try shooting with v as starting vector
    [tt cc vv err] = shoot(p1,p2,F,DF,D2F,ltol,40,v);
    if isequal(err,false)
        t = tt;
        c = cc;
        v = vv;
    end
end

function [v solExp] = Log(p0,p1,varargin)
    ltol = tol;
    if size(varargin,2) >= 1
         ltol = varargin{1};
    end
    if size(varargin,2) >= 2
         lguess = varargin{2};
    end    
    
    if ~exist('lguess','var')
        [v] = shoot(p0,p1,m,n,F,DF,D2F,@Exp,@DExp,ltol);
    else
        [v] = shoot(p0,p1,m,n,F,DF,D2F,@Exp,@DExp,ltol,lguess);
    end
    
%     % debug - note that p0,p1 might be slightly off manifold causing below check to fail even if shooting worked ok
%     if norm(p1-Exp(p0,v,ltol)) > ltol
%         fname = ['Log-error-' strrep(char(java.util.UUID.randomUUID),'-','_') '.mat'];
%         save(fname,'p0','p1','v','ltol');  
%         assert(false);        
%     end

%     % preliminary multiple shooting    
%     if not(isequal(errGuess,false))
%         assert(false); % untested with new IVP structure
%         pmiddle = errGuess;
%         [tt cc vv err] = shoot(pmiddle,p2,F,DF,D2F,lExp,ldExp,ltol);        
%         if not(isequal(err,false))
%             fname = ['shooting-multiple-failure-' strrep(char(java.util.UUID.randomUUID),'-','_') '.mat'];
%             save(fname,'p1','p2','pmiddle');            
%             assert(false);
%         end
%         [v t c] = lLogInit(p1,p2,[c(:,1:end-1) cc],ltol);
%     end
end

    % find point on M reasonably close to x 
    function px = toManifold(x)
        function [y dy] = lf(x)
            y = F(x);
            dy = DF(x);
        end
        
        projOptions = optimset('Jacobian','on','Display','off','TolFun',tol*1e-4,'TolX',tol*1e-4);
        warning off optim:fsolve:NonSquareSystem
        px = fsolve(@lf,x,projOptions);
        assert(norm(F(px)) < epsilon());
    end

    function res = isTangent(V,p)
        res = norm(DF(p)*V) < epsilon();
    end

    function res = isPoint(p)
        res = norm(F(p)) < epsilon();
    end

    function frame = orthFrame(p)
        frame = null(DF(p));
    end

manifold.dim = m-n;
assert(manifold.dim >= 1);

manifold.type = 'embedded';
manifold.F = F;
manifold.DF = DF;
manifold.D2F = D2F;
manifold.m = m;
manifold.n = n;

manifold.dist = @distM;
manifold.Exp = @Exp;
manifold.DExp = @DExp;
% manifold.DExpSolve = @DExpSolve;
manifold.D2Exp = @D2Exp;
manifold.Pt = @Pt;
manifold.Log = @Log;
manifold.isTangent = @isTangent;
manifold.isPoint = @isPoint;
manifold.orthFrame = @orthFrame;
manifold.getExp = @(sol,t) getExp(sol,DF,t);
manifold.getDExp = @(sol,t) getDExp(sol,m,t);
manifold.toManifold = @toManifold;

end
