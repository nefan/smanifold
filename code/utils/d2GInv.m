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

function d2GInvA = d2GInv(A,dtA,dsA,dsdtA,GInvA)
%
% second derivative of generalized inverse, confer Decell '72
%

assert(all(size(A) == size(dtA)));
assert(all(size(A) == size(dsA)));
assert(all(size(A) == size(dsdtA)));

function R = Lambda(A,B,GInvA)
    R = dGInv(A,GInvA,B);
end

function R = tLambda(A,B,C,D,GInvA,LambdaAC)
    XXAB = XX(A,B,GInvA);
    YYABCD = YY(A,B,C,D,GInvA,LambdaAC);
    AGInvA = A*GInvA;
    GInvAA = GInvA*A;
    R=-LambdaAC*B*GInvA-GInvA*D*GInvA-GInvA*B*LambdaAC+YYABCD...
        -(LambdaAC*A+GInvA*C)*XXAB*AGInvA...
        -GInvAA*YYABCD*AGInvA...
        -GInvAA*XXAB*(C*GInvA+A*LambdaAC);    
end

function R = XX(A,B,GInvA)
    R = B'*GInvA'*GInvA+GInvA*GInvA'*B';
end

function R = YY(A,B,C,D,GInvA,LambdaAC)
    R = D'*GInvA'*GInvA+B'*(LambdaAC'*GInvA+GInvA'*LambdaAC)...
        +(LambdaAC*GInvA'+GInvA*LambdaAC')*B'+GInvA*GInvA'*D';
end
    
d2GInvA = tLambda(A,dtA,dsA,dsdtA,GInvA,Lambda(A,dsA,GInvA));
% debug - check commutativity
%assert(all(all(d2GInvA == tLambda(A,dsA,dtA,dsdtA,GInvA,Lambda(A,dtA,GInvA)))));

end