function [dir residual] = shootalg(p0,p1,Np0,m,n,DF,Exp,DExp,tol,init,alg,verbose)
%
% shooting algorithm
% (basically just Gauss-Newton or gradient descent)
%

maxIter = 200;

i = 0;

inttol = tol;
reduceStep = false;
step = 1.0;
minStep = step/2^8;
dimM = m-n;


% detect divergent oscillations - we reduce step size if decrease in length
% is less than 1-stepDecreateFactor percent and angle varies a lot
stepDecreaseFactor = 1-0.02;
angleVariationFactor = cos(3*pi/4);

[guess v solExp] = Exp(p0,Np0*init,inttol);
frame = DExp(solExp,Np0,inttol);
residual = p1-guess; % residual in embedding space

lenRes = norm(residual);
lenStep = 0;
angleVariation = 0;
correction = zeros(dimM,1);
lastCorrection = correction;
dir = init;

% for debugging
if verbose
    formatstr          = ' %5.0f       %5.0f     %13.6g  %12.3g %12.6g       %12.6g\n';
    fprintf( ...
        ['\n ', alg, ' optimization:' ...
        '\n                                                                           Norm of \n', ...
           '   Iteration  Func-count    Residual       angleVar.       Step            step\n']);
    
    fc = 1;
end

while lenRes > tol && i <= maxIter
    if verbose
        fprintf(formatstr,i,fc,lenRes,angleVariation,step,lenStep);
    end
   
    if reduceStep         
        if lenRes > lastLenRes || (step == 1.0 && lenRes > lastLenRes*stepDecreaseFactor && angleVariation < angleVariationFactor)
            % reduce step size          
            step = step*0.5;
            if step < minStep % no luck
                if verbose
                    fprintf('warning: exiting %s optimization since lower step limit reached\n',alg);
                end
                return;
            end
            guess = lastGuess;
            frame = lastFrame;
            dir = lastDir;
            correction = lastCorrection;
            lenRes = lastLenRes;
            %inttol = 1e-3*tol;
            reduceStep = false;
            i = i + 1;            
            continue;
        else
            step = 1.0; % reset step            
            %inttol = tol;
        end
    else
        reduceStep = true;
    end
    
    if strcmp(alg,'GD')
        pRes = frame'*residual; % gradient descent
    else % GN
        % coordinates in tangent frame
        [U S V] = svd(frame);
        U1 = U(:,1:size(S,2));
        pRes = V * diag(1./diag(S)) * U1'*residual;        
    end
    lenStep = norm(pRes);
    
    if lenStep < tol
        % choose random direction?
        if verbose % residual in embedding space
            fprintf('warning: exiting %s optimization since lower gradient limit reached (zero projection)\n',alg);
        end        
        return;
    end
    
    lastCorrection = correction;
    correction = pRes; % residue transported to Np0    
    lastDir = dir;
    dir = dir + step*correction;

    lastGuess = guess;    
    lastFrame = frame;
    [guess v solExp] = Exp(p0,Np0*dir,inttol);
    frame = DExp(solExp,Np0,inttol);
    if verbose
        fc = fc + 1;
    end
    
    residual = p1 - guess;
    lastLenRes = lenRes;
    lenRes = norm(residual);    
    angleVariation = dot(correction,lastCorrection)/(norm(correction)*norm(lastCorrection));    
    
    i = i + 1;
end

if verbose
    fprintf(formatstr,i,fc,lenRes,angleVariation,step,lenStep);    
    
    if lenRes < tol
        fprintf('%s optimization: convergence after %d iterations\n',alg,i);
    end
    if i > maxIter
        fprintf('warning: exiting %s optimization since iteration limit reached\n',alg);
    end
end
