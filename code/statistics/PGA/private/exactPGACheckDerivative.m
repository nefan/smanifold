function exactPGACheckDerivative(v,V,Vk,B,mode,data,Fproj,projTol,fval,g,tol)
%
% finite difference derivative check
%

epsilon = 10e-5;

fprintf('Forward difference check. Gradient norm: %e\n',norm(g));

function y = gamma(t,w,v)
    assert(abs(norm(v) - 1) < epsilon);
    assert(abs(norm(w) - 1) < epsilon);  
    y = cos(t)*v+sin(t)*w;
end
    
dt = 1e-4; % step    
Vvp = null([Vk v]');
gfdVvp = zeros(size(Vvp,2),1);
gVvp = Vvp'*g; % debug
for i = 1:size(Vvp,2)
	w = Vvp(:,i); % direction        
    assert(norm(gamma(0,w,v) - v) < epsilon);
        
    % forward difference
    gfdVvp(i,1) = (exactPGAF(data,gamma(dt,w,v),V,Vk,B,mode,Fproj,projTol,true) - fval)/dt;        
    fprintf('Forward difference error %d (step %e, absolute/relative): %e/%e\n',i,dt,abs(gfdVvp(i)-gVvp(i)),abs((gfdVvp(i)-gVvp(i))/gVvp(i)));
end    

gfd = Vvp*gfdVvp;
    
g % debug
gfd % debug
    
fprintf('Forward difference error (step %e, absolute/relative/angular): %e/%e/%e\n',dt,max(abs(g-gfd)),max(abs(g-gfd)./gfd),acos(dot(g,gfd)/(norm(g)*norm(gfd)))*360/(2*pi));

end
