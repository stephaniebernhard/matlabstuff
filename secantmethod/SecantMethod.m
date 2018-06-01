function [ y ] = SecantMethod( f, x0, x1, tol)
% [y] = SecantMethod(@(x) exp(x.^2)+x^-3-10, 2, 1.8, 10^-3)
    sol = fzero(f,x0);
    y=x0;
    while (abs(y-sol)>tol)
        y = x1-((x1-x0)/(f(x1)-f(x0)))*(f(x1));
        x0=x1;
        x1=y;
    end
end

