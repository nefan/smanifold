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

function [distM Exp Log LogInit DExp D2Exp Pt DExpSolve] = ExpLogMaps(m,n,F,DF,D2F,tol)
%
% Wrapper function exporting Exp and Log maps
% for submanifold defined by F
%
% Function LogInit behaves as Log except that
% it uses an initial curve joining the points
% to overcome limitations of shooting method.
%

%addpath('shooting');

epsilon = 10e-5; % this shouldn't be hardcoded

% Exponential map
function [x v sol] = lExp(x0,v0,varargin)
    ltol = tol;
    tspan = [0 1];
    if size(varargin,2) >= 1
         ltol = varargin{1};
    end
    if size(varargin,2) >= 2
         tspan = [0 varargin{2}];
    end    
    sol = intExp(x0,v0,tspan,m,n,F,DF,D2F,ltol);
    [x v] = getExp(sol,DF,tspan(2));
end

% DExp
function [B sol] = lDExp(solExp,B0,varargin)
    ltol = tol;
    tspan = [0 1];    
    if size(varargin,2) >= 1
         ltol = varargin{1};
    end
    if size(varargin,2) >= 2
         tspan = [0 varargin{2}];
    end        
    sol = intDExp(B0,tspan,solExp,m,n,F,DF,D2F,true,ltol);        
    B = getDExp(sol,m,tspan(2));
end

% % DExpSolve
% function [P0] = lDExpSolve(solExp,S,varargin)
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
function [B sol] = lD2Exp(solExp,B0,U,solw,solu,varargin)
    ltol = tol;
    if size(varargin,2) >= 1
         ltol = varargin{1};
    end
    if isempty(solw)
        solw = intDExp(B0,[0 1],solExp,m,n,F,DF,D2F,true,ltol);
    end
    if isempty(solu)
        solu = intDExp(U,[0 1],solExp,m,n,F,DF,D2F,true,ltol);
    end
    sol = intD2Exp(solExp,solw,solu,m,n,F,DF,D2F,true,ltol);
    B = getD2Exp(sol,m,size(B0,2),size(U,2));
end
% debug lines
% for orthogonal
%clf, grid on, hold on, plot(solw.x,solw.y(2,:),'r'), plot(solw.x,solw.y(2,:)./solw.x,'k'), plot(0.5*(solw.x(2:end)+solw.x(1:end-1)),diff(solw.y(2,:)./solw.x)./(4*pi*diff(solw.x)),'m'), plot(sol.x,sol.y(5,:)./sol.x.^2,'g'), plot(sol.x,sol.y(5,:),'b')
% for parallel
%global ii, global ll, global llw, global llt, global dotll, ll = []; llw = []; llt = []; dotll = []; for ii = 1:length(sol.x); ll(ii) = norm(sol.y(4:6,ii)); llt(ii) = ll(ii)/sol.x(ii)^2; dotll(ii) = dot(sol.y(4:6,ii),getDExp(solw,m,sol.x(ii))); end, for ii = 1:length(solw.x); llw(ii) = norm(solw.y(1:3,ii)); end, clf, grid on, hold on, plot(solw.x,llw,'r'), plot(solw.x,llw./solw.x,'k'), plot(0.5*(solw.x(2:end)+solw.x(1:end-1)),diff(llw./solw.x)./(4*pi*diff(solw.x)),'m'), plot(sol.x,llt.^(-1).*dotll,'g'), plot(sol.x,ll,'b')
% untested
%global ii, global ll, global llw, global llt, global dotll, global llwt, ll = []; llw = []; llt = []; dotll = []; llwt = []; for ii = 1:length(sol.x); ll(ii) = sum(sol.y(4:6,ii).^2); llt(ii) = ll(ii)/sol.x(ii)^4; dotll(ii) = dot(sol.y(4:6,ii),getDExp(solw,m,sol.x(ii))); end, for ii = 1:length(solw.x); llw(ii) = sum(solw.y(1:3,ii).^2); llwt(ii) = llw(ii)/solw.x(ii)^2; end, clf, grid on, hold on, plot(solw.x,llw,'r'), plot(solw.x,llwt,'k'), plot(0.5*(solw.x(2:end)+solw.x(1:end-1)),diff(llwt)./(4*pi*diff(solw.x)),'m'), plot(sol.x,2*dotll,'g'), plot(sol.x,ll,'b')

% parallel transport
function [B sol] = lPt(solExp,B0,varargin)
    ltol = tol;
    if size(varargin,2) >= 1
         ltol = varargin{1};
    end
    sol = intPt(B0,solExp,m,n,F,DF,D2F,ltol);
    B = getPt(sol,m);
end

function d = ldistM(p1,p2)
    v = lLog(p1,p2);
    d = norm(v);
end

function [v t c] = lLogInit(p1,p2,init,varargin)
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

function [v solExp] = lLog(p0,p1,varargin)
    ltol = tol;
    if size(varargin,2) >= 1
         ltol = varargin{1};
    end
    if size(varargin,2) >= 2
         lguess = varargin{2};
    end    
    
    if ~exist('lguess','var')
        [v] = shoot(p0,p1,m,n,F,DF,D2F,@lExp,@lDExp,ltol);
    else
        [v] = shoot(p0,p1,m,n,F,DF,D2F,@lExp,@lDExp,ltol,lguess);
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

distM = @ldistM;
Exp = @lExp;
DExp = @lDExp;
%DExpSolve = @lDExpSolve;
D2Exp = @lD2Exp;
Pt = @lPt;
Log = @lLog;
LogInit = @lLogInit;

end
