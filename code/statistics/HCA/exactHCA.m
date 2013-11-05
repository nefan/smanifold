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

function [V,s,coords] = exactHCA(data,mu,nr,B,tol,manifold,Vapprox)
%
% Compute exact HCA
%
% data contains points on manifold
% mu is basepoint (usually mean)
% nr indicates nr of principal directions to compute
% B orthonormal basis for T_\mu M
%

epsilon = 10e-5;

assert(isOrthonormal(B));

% debug
if ~isempty(whos('global','debug'))
    global debug;
else
    debug = false;
end
global prepend;
printFig = @(name) print([prepend name '.ps'],'-dpsc');

% mode -reconstruction error
mode = 'R'; % reconstruction error R or variance V
sign = 1;
descent = 'GN'; % Gauss-Newton GN or gradient descent GD

N = size(data,2);
assert(nr <= manifold.dim);

% tolerances
projTol = 1e-3*tol/N;
gradTol = 1e-3*tol/N;
minProjTol = 1e-7*tol/N;
minGradTol = 1e-7*tol/N;

% projection
function r = Fproj(x,v,Vk,datakBV,ltol)
    if ~isempty(Vk)
        BVk = B*Vk;
    else 
        BVk = [];
    end    
    [xk,Bxkv,BxkVk] = HCALinToMan(datakBV,mu,BVk,B*v,manifold);
    [y,w,R,Logxky,w0] = geodesicSubspaceProjection(x,xk,Bxkv,[],ltol,manifold);
    w = Bxkv'*w; % to v basis
    linR = sum(manifold.Log(manifold.Exp(xk,w0),x).^2);
    r = {y,w,R,Logxky,linR};
end

% gradient
function r = Fgrad(x,datakBV,y,w,Logxy,Vk,v,Vvp,ltol)   
    if ~isempty(Vk)
        BVk = B*Vk;
    else 
        BVk = [];
    end    
    [xk,M,BxkVk] = HCALinToMan(datakBV,mu,BVk,[B*v B*Vvp],manifold);    
    Bxkv = M(:,1);
    BxkVvp = M(:,2:end);
    [J,g,Bx] = geodesicSubspaceLogGradient(x,xk,y,Bxkv*w,Logxy,BxkVk,Bxkv,BxkVvp,ltol,mode,manifold);
    r = {J,g,Bx};
end

Logx = [];
parfor j = 1:N
    Logx(:,j) = B'*manifold.Log(mu,data(:,j),tol);
end

dataBV = zeros(nr,N);
V = []; % orthonormal basis for V_k
s = []; % variances
% debug
if debug
    global timeProj;
    global timeGrad;
    timeProj = 0;
    timeGrad = 0;
end
for k = 0:nr-1
    Vk = V(:,1:k); % actually equal to V, but we use Vk in the loop    
    Vkp = []; % basis for V_k^\perp
    if k == 0
        Vkp = eye(manifold.dim);
    else
        Vkp = null(Vk');
    end
    
    LogxVkp = Vkp'*Logx; % project to Vkp
    [VVkp skp] = PCA(LogxVkp,false); % regular PCA
    v = Vkp*VVkp(:,end); % initial guess    
%     v = [0 1]'; % for illustration
%     v = 1/sqrt(2)*[1 1]'; % for illustration
%     if k == 1
%         v = 1/sqrt(2)*[1 1 0]'; % debug
%     end
    
    p = 0.5; % descent direction decrease factor
    c = 0.25; % for Armijo condition
    alphahat = 0.25; % 1; often works better
%     alphahat = 4.0; % for illustration
    alpha = alphahat;
    limitFactor = 1.0;
    firstRun = true;
    reRun = false;
    stepErr = inf;
    valErr = inf;
    prevfval = inf;
    prevv = [];
    prevg = [];
    gn = inf;
    prevgn = 0;
    i = 0;
    maxIter = 400;
    minIter = 1;
%     maxIter = 1; % for illustration
%     minIter = 1;   
    
    while (i < minIter || (~firstRun && prevgn < gn) || min(stepErr,valErr) > tol) && i < maxIter % loop until convergence    
        
        % projections
        if k > 0
            datak = B*V*dataBV(1:k,:);
        else
            datak = B*zeros(manifold.dim,N);
        end
        [fval ys ws Logxys rs linfval xxRs] = exactHCAF(data,dataBV(1:k,:),v,Vk,@Fproj,projTol,debug);

        if not(firstRun || reRun) % Armijo condition            
            assert(c*alpha*limitFactor*prevg'*descentDir < epsilon); % descent condition
            if sign*fval > sign*prevfval+c*alpha*limitFactor*prevg'*descentDir               
                % backtrack
                alpha = p*alpha % debug output
                NN = 8;
                if alpha <= p^NN*alphahat;
                    projTol = projTol*0.1;
                    gradTol = gradTol*0.1;
                    if projTol < minProjTol || gradTol < minGradTol
						figure(4)
                        testres(6,prevv,Vk,p,alpha,sign,descentDir);
                        printFig('TMassert');
						close(4)
                        assert(false);
                    end
                    alpha = alphahat;
                    reRun = true;
                    prevfval = inf;
                    fprintf('lowering tolerances: %e, %e\n',projTol,gradTol);
                    v = prevv;
                    continue;
                end
                v = prevv + alpha*limitFactor*sign*descentDir;
                v = v/norm(v); % project to unit sphere        
                continue;
            else
                alpha = alphahat;
            end
        end
        
        if k == manifold.dim-1 % nothing to do
            stepErr = 0;
            valErr = 0;
            prevv = v;
            prevws = ws;
            prevfval = fval;
            minIter = 0;
            continue;
        end

        % gradients
        [g Js gs Vvp rs] = exactHCAFgrad(data,dataBV(1:k,:),v,ys,ws,Logxys,rs,Vk,B,k,mode,@Fgrad,gradTol,debug);  

        g = 1/N*g; % seems to work nicely :-)
        descentDir = zeros(manifold.dim,1);
        if strcmp(descent,'GD')
            descentDir = -g;
        else if strcmp(descent,'GN')
                [UU SS VV] = svd(Js);
                S = [];
                for r = 1:min(size(SS)) % extract non-zero singular values
                    ss = SS(r,r);
                    if abs(ss) < 10e-5
                        break;
                    end
                    S(r) = ss;
                end
                UU1 = UU(:,1:length(S));
                descentDir = -1/N*Vvp*VV*diag(1./S)*UU1'*reshape(rs,manifold.dim*N,1); % debug on sign
            else
                assert(false);
            end
        end
        
        % debug
        limitCorrection = true;
        limitLength = 1.0;
        limitFactor = 1.0;
        if limitCorrection && norm(descentDir) > limitLength
            limitFactor = limitLength/norm(descentDir);
        end
        
        % error update
        valErr = abs(fval-prevfval);
        
        % save        
        prevfval = fval;
        prevv = v;
        prevg = g;
        prevgn = gn;
        prevws = ws;
        
        % update        
        v = v + alpha*sign*limitFactor*descentDir;
        v = v/norm(v); % project to unit sphere        
        gn = norm(descentDir);
        
        % error update
        stepErr = norm(v-prevv);        
        
        % control
        firstRun = false;
        reRun = false;
        i = i+1;
        % debug
        if true % debug
            fprintf('it %d,%d: %e, %e, %e, %e\n',...
                k,i,fval,gn,valErr,stepErr);
        end
    end
    
    dataBV(k+1,:) = prevws;
    V(:,k+1) = prevv;
%     V(:,k+1) = v; % for illustration
    s(k+1) = prevfval;

    assert(isOrthonormal(V));
    global outputDir;
    save([outputDir '/hca-red-it-' int2str(k) '.mat'],'V','s','dataBV')
end

% shift back to R^m
V = B*V;

% coordinates
coords = dataBV;

% debug
if debug
    fprintf('timeProj %f\n',timeProj);
    fprintf('timeGrad %f\n',timeGrad);
end

end
