function exactPGA2DimVis(mu,B,v,ws,gs,Logx,it,N,manifold)

global prepend;
printFig = @(name) print([prepend name '.ps'],'-dpsc');

figure(3), clf, hold on
plot(Logx(1,:),Logx(2,:),'r*');
quiver(0,0,v(1),v(2),'k','LineWidth',1);
            
for j = 1:N
    w = ws(:,j);
    %yyy = ys(:,j);
    %xxx = data(:,j);
    g = gs(:,j);
    %RRR = Rs(:,j);
             
    %vv = Log(y,x);
    %assert(abs(norm(vv)-Rs(:,j)) < tol);
    %[xx vv solExp] = Exp(y,vv);
    %for t = 0:1/10:1
    %    p = B'*Log(mu,getExp(solExp,DF,t));
    %    plot(p(1),p(2),'k.');
    %end                
                
    plot(w(1),w(2),'kx');
    quiver(w(1),w(2),g(1),g(2));           
end

axis equal
hold off
% printFig(['TMit' num2str(i)]);
%close(3);

figure(1), hold on
axis([-2.5 2.5 -2.5 2.5 -2.5 2.5])
% for PGAalgIllustration
if it == 25        
    for i = 1:N
        xi = manifold.Exp(mu,B*Logx(:,i));
        pxi = manifold.Exp(mu,B*ws(:,i));
        plot3(pxi(1,:),pxi(2,:),pxi(3,:),'bo','MarkerSize',10,'MarkerFaceColor','b');
        [xx vv sol] = manifold.Exp(xi,manifold.Log(xi,pxi));
        xx = [];
        for t = [0:0.01:1]
            xx = [xx manifold.getExp(sol,t)];
        end
        plot3(xx(1,:),xx(2,:),xx(3,:),'--k','LineWidth',1);
        
        [xx vv sol] = manifold.Exp(mu,B*2*pi*v);
        xx = [];
        for t = [0:0.05:1]
            xx = [xx manifold.getExp(sol,t)];
        end
        plot3(xx(1,:),xx(2,:),xx(3,:),'k','LineWidth',2.5);        
    end
    quiver3(0,0,1,1,0,0,'r','LineWidth',4.0)
    quiver3(0,0,1,0,1,0,'k','LineWidth',3.0)
end
if it == 0 % sphere illustration
    [xx vv sol] = manifold.Exp(mu,B*2*pi*[0 1]');
    xx = [];
    for t = [0:0.01:1]
        xx = [xx manifold.getExp(sol,t)];
    end
    plot3(xx(1,:),xx(2,:),xx(3,:),'--k','LineWidth',2.5);    
end

end