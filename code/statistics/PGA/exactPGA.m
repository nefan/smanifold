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

function [V,s,sapprox,sfletcher,estimate,RsDiff] = exactPGA(data,mu,nr,B,tol,manifold,Vapprox)
%
% Compute exact PGA
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

% mode - variance or reconstruction error
if ~isempty(whos('global','mode'))
    global mode;
else
    mode = 'V'; % reconstruction error R or variance V
end
fprintf('PGA mode: %s\n',mode);
sign = 1;
if mode == 'V'
    sign = -1;
end
descent = 'GD'; % Gauss-Newton GN or gradient descent GD

N = size(data,2);
assert(nr <= manifold.dim);

% tolerances
projTol = 1e-3*tol/N;
gradTol = 1e-3*tol/N;
minProjTol = 1e-7*tol/N;
minGradTol = 1e-7*tol/N;

% projection
function r = Fproj(x,v,Vk,ltol)    
    [y,w,R,Logxy,w0] = geodesicSubspaceProjection(x,mu,B*[Vk v],[],ltol,manifold);
    w = B'*w; % to B basis
    linVar = sum(w0.^2);
    linR = sum(manifold.Log(manifold.Exp(mu,w0),x).^2);
    r = {y,w,R,Logxy,linVar,linR};   
end

% gradient
function r = Fgrad(x,y,w,Logxy,BVk,Bv,BVvp,ltol)    
    [J,g,Bx] = geodesicSubspaceLogGradient(x,mu,y,B*w,Logxy,BVk,Bv,BVvp,ltol,mode,manifold);
    r = {J,g,Bx};
end

% expectation
function e = expectedPGADiff(mu,x,lLogx,VV)
    assert(isOrthonormal(VV));
    BVV = B'*VV;
    w = lLogx;
    pw = BVV*BVV'*w; % projection to subspace
    rw = w-pw; % residual
    
    [y,v,solExp] = manifold.Exp(mu,B*pw);

    % eR
    LeJ = norm(rw); % length of Euclidean Jacobi field
            
    LJ1 = norm(manifold.DExp(solExp,B*rw));
    [xx,v,solExp2] = manifold.Exp(mu,B*w);
    LJ2 = norm(manifold.DExp(solExp,B*rw));
    %lexpR = 0.5*(LJ1+LJ2);
    lexpR = norm(manifold.Log(y,x));

    % eV
    F = manifold.DExp(solExp,VV);
    g = -BVV*2*F'*manifold.Log(y,x);
    
    lexpV = norm(pw+g);
    
    e = {norm(rw), lexpR, norm(pw), lexpV, norm(g)};
end

Logx = [];
parfor j = 1:N
    Logx(:,j) = B'*manifold.Log(mu,data(:,j),tol);
end

% measure expected difference
fprintf('Expectation on difference...\n');
Vexpect = Vapprox(:,1);
expR = [];
expV = [];
gV = [];
R = [];
V = [];
fun = @expectedPGADiff;
parfor j = 1:N
    res = fun(mu,data(:,j),Logx(:,j),Vexpect);
    R(j) = res{1};
    expR(j) = res{2};
    V(j) = res{3};
    expV(j) = res{4};
    gV(j) = res{5};
end
R % debug
expR % debug
%diffR = (R.^2-expR.^2)./sum(R.^2);
diffR = 1-expR./R;
%diffV = (V.^2-expV.^2)./sum(V.^2);
%diffV = gV.^2/sum(V.^2);
diffV = gV./V;
estimate = [mean(diffR),std(diffR),mean(diffV),std(diffV),mean(R-expR),std(R-expR),mean(gV),mean(R),mean(V)];
% fprintf('Expectation (rel diff R/rel std diff R/rel diff V/diff R/std diff R/diff V/R/V):\n\t%e/%e/%e/%e/%e/%e/%e/%e\n', ...
%     mean(diffR),std(diffR), ...
%     mean(diffV), ...    
%     mean(R-expR),std(R-expR), ...
%     mean(gV),mean(R),mean(V));
fprintf('Expectation (rho/sigma):\n\t%e/%e\n',mean(gV),std(R-expR));
clear('expR','expV','gV','R','V');

% measure approximated s
fprintf('Evaluating fletcher PGA...\n');
sapprox = []; % true projection Fletcher method variances
sfletcher = []; % plain Fletcher method variances
for k = 0:nr-1    
    Vk = B'*Vapprox(:,1:k);
    v = B'*Vapprox(:,k+1);
    [fval ys ws Logxys rs linfval xxRs] = exactPGAF(data,v,Vk,mode,@Fproj,projTol,debug);
    if k == 0
        linRs = xxRs;
    end
    sapprox(k+1) = fval;
    sfletcher(k+1) = linfval;
end
sapprox % debug
sfletcher % debug

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
    v = 1/sqrt(2)*[1 1]'; % for illustration
    
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
    
    % initialization from previous runs
    if k == 0 && ~isempty(whos('global','exactPGAInit'))
        global exactPGAInit;
        i = exactPGAInit.i;
        v = exactPGAInit.v;
        prevv = exactPGAInit.prevv;
        g = exactPGAInit.g;
        prevg = exactPGAInit.prevg;
        fval = exactPGAInit.fval;
        prevfval = exactPGAInit.prevfval;
        alpha = exactPGAInit.alpha;
        limitFactor = exactPGAInit.limitFactor;
        descentDir = exactPGAInit.descentDir;
        projTol = exactPGAInit.projTol;
        gradTol = exactPGAInit.gradTol;
		approxPGAfval = exactPGAInit.approxPGAfval;
		approxPGAv = exactPGAInit.approxPGAv;
        assert(all(all(exactPGAInit.data == data)));
        assert(exactPGAInit.tol == tol);
        assert(exactPGAInit.sign == sign);
        assert(exactPGAInit.mode == mode);
		firstRun = false;
		reRun = false;
		fprintf('Restarting old run iteration %d, old fval: %e\n',i,exactPGAInit.prevfval);
    end    
    
    while (i < minIter || (~firstRun && prevgn < gn) || min(stepErr,valErr) > tol) && i < maxIter % loop until convergence    
        
        % projections
        [fval ys ws Logxys rs linfval xxRs] = exactPGAF(data,v,Vk,mode,@Fproj,projTol,debug);
        if k == 0
            exactRs = xxRs;
        end

        if (firstRun) % debug
            approxPGAfval = fval;
            approxPGAv = v;
        end
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
            minIter = 0;
            continue;
        end

        % gradients
        [g Js gs Vvp] = exactPGAFgrad(data,v,ys,ws,Logxys,rs,Vk,B,k,mode,@Fgrad,gradTol,debug);
        
        % debug        
        if debug && manifold.dim == 2
            %if firstRun % for illustration
                exactPGA2DimVis(mu,B,v,ws,gs,Logx,i,N,manifold)
            %end
        end     

        %Js = sqrt(2/N)*Js; % seems to work better...
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
                descentDir = -Vvp*VV*diag(1./S)*UU1'*reshape(rs,manifold.dim*N,1); % debug on sign
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
            fprintf('it %d,%d: %e, %e, %e, %e, %e, %e, %e\n',...
                k,i,fval,gn,valErr,stepErr,abs(fval-approxPGAfval),abs(fval-approxPGAfval)*100/approxPGAfval,acos(abs(dot(prevv,approxPGAv)))*360/(2*pi));
%             fname = [prepend 'exact-PGA-k-' int2str(k) '-i-' int2str(i) '.mat'];
%             save(fname,'data','mu','i','Vk','v','B','ys','ws','prevv','g','prevg','fval','prevfval','descentDir','sign','alpha','limitFactor','tol','approxPGAfval','approxPGAv','projTol','gradTol','mode');                         
        end
    end
        
    % debug
    if debug && manifold.dim == 2 && k < nr-1
		figure(4)
        if manifold.dim == 2
            dir = [v(2); -v(1)];
            dir = dir/norm(dir);
        else
            dir = descentDir/norm(descentDir);
        end
%         testres(5,prevv,Vk,p,alpha,sign,dir);
%         printFig(['TMdone' num2str(k)]);
		close(4)
    end 
    
    V(:,k+1) = prevv;
%     V(:,k+1) = v; % for illustration
    s(k+1) = fval;

    assert(isOrthonormal(V));
end

% shift back to R^m
V = B*V;

% statistics
%RsDiff = (sqrt(linRs)-sqrt(exactRs))./sqrt(linRs);
RsDiff = (linRs-exactRs)./linRs;

% debug
if debug
    fprintf('timeProj %f\n',timeProj);
    fprintf('timeGrad %f\n',timeGrad);
end

end
