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

function expectedPGADiff(data,mu,Logx,nr,B,tol,manifold,Vapprox)
%
% Compute expectation on difference (ECCV2010 paper)
%

assert(isOrthonormal(B));

% debug
if ~isempty(whos('global','debug'))
    global debug;
else
    debug = false;
end

% mode - variance or reconstruction error
if ~isempty(whos('global','mode'))
    global mode;
else
    mode = 'V'; % reconstruction error R or variance V
end

N = size(data,2);

% expectation
function e = ediff(mu,x,lLogx,VV)
    assert(isOrthonormal(VV));
    BVV = B'*VV;
    w = lLogx;
    pw = BVV*BVV'*w; % projection to subspace
    rw = w-pw; % residual
    
    [y v solExp] = manifold.Exp(mu,B*pw);

    % eR
    LeJ = norm(rw); % length of Euclidean Jacobi field
            
    LJ1 = norm(manifold.DExp(solExp,[],B*rw));
    [xx v solExp2] = manifold.Exp(mu,B*w);
    LJ2 = norm(manifold.DExp(solExp,[],B*rw));
    %lexpR = 0.5*(LJ1+LJ2);
    lexpR = norm(manifold.Log(y,x));

    % eV
    F = manifold.DExp(solExp,[],VV);
    g = -BVV*2*F'*manifold.Log(y,x);
    
    lexpV = norm(pw+g);
    
    e = {norm(rw), lexpR, norm(pw), lexpV, norm(g)};
end


% measure expected difference
fprintf('Expectation on difference...\n');
Vexpect = Vapprox(:,1);
expR = [];
expV = [];
gV = [];
R = [];
V = [];
fun = @ediff;
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

end