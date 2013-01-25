function m = intrinsicMean(xi,manifold,tol,varargin)
%
% Compute the intrinsic (Frechet) mean of xi
%
% At this point we assure nothing for non-localized data.
% Returns empty matrix on error.
%

step = 1.0; % we may need to change it during iteration
maxIter = 400;

N = size(xi,2); % number of points
if size(varargin,2) >= 1
    m = varargin{1}; % initial guess
else
    m = xi(:,1); % initial guess
end
i = 1;

while i <= maxIter
    pm = m;
    
    % iterate
    ss = zeros(manifold.m,N);
    parfor j = 1:N
        ss(:,j) = manifold.Log(m,xi(:,j));
    end
    s = sum(ss,2);
    m = manifold.Exp(m,step/N*s);
    
    err = norm(pm - m);
    fprintf('mean computation err (iteration %d): %f\n',i,err);
    if err < tol
        break;
    end
    i = i + 1;
    
    if i == 30 % primitive, change to real line search
        step = 0.5;
    end
end

assert(i <= maxIter);
