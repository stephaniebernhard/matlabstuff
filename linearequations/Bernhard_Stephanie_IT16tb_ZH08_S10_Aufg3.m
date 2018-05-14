function [ xn, n, n2 ] = Bernhard_Stephanie_IT16tb_ZH08_S10_Aufg3( A, b, x0, tol, opt )
% Given a matrix A, a vector b, a starting vector x0, a maximum error tol
% and an opt (Jacobi, Gauss-Seidel), this function returns a vector xn
% which is the approximate solution after n steps and return n2, the number
% of iterations according to a-priori estimate
% Funktionsaufruf:
% [xn,n,n2]=Bernhard_Stephanie_IT16tb_ZH08_S10_Aufg3([8,5,2;5,9,1;4,2,7],[19;5;34],[1;-1;3],10^-4,"Jacobi")
% [xn,n,n2]=Bernhard_Stephanie_IT16tb_ZH08_S10_Aufg3([8,5,2;5,9,1;4,2,7],[19;5;34],[1;-1;3],10^-4,"Gauss-Seidel")
% transpose vector if needed
[m1,m2] = size(b);
dim=length(A);
if m1 == 1
    b=b';
end
[m3,m4]= size(x0);
if m3 == 1
    x0=x0';
end
% initialize decomposition
D = zeros(dim);
L = zeros(dim);
R = zeros(dim);
% decompose A = L + D + R
for i=1:dim
    D(i,i) = A(i,i);
    for j=1:i-1
        L(i,j) = A(i,j);
    end
    for k=i+1:dim
        R(i,k) = A(i,k);
    end
end
xold = x0;
n=1;
if opt=="Jacobi"
    % calculate x1
    x1 = -D^-1*(L+R)*x0+ D^-1*b;
    xnew = x1;
    while norm(xold-xnew,inf) > tol
        xold = xnew;
        xnew = -D^-1*(L+R)*xnew + D^-1*b;
        n=n+1;
    end
    B = -D^-1*(L+R);
elseif opt=="Gauss-Seidel"
    x1 = -(D+L)^-1*R*x0+(D+L)^-1*b;
    xnew = x1;
    while norm(xold-xnew,inf) > tol
        xold = xnew;
        xnew = -(D+L)^-1*R*xnew+(D+L)^-1*b;
        n=n+1;
    end
    B = -(D+L)^-1*R;
else
    error('option given is not valid');
end
xn=xnew;
% estimate number of iterations
n2=log((tol-norm(B,inf)*tol)/norm(x1-x0,inf))/log(norm(B,inf));
end

