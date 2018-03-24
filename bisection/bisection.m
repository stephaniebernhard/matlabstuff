function [ root,xit,n] = Bernhard_Stephanie_IT16tb_ZH08_S4_Aufg2( func, a, b, tol)
% Funktionsaufruf:
    if func(a)*func(b)>=0
        return
    end
    xit = [];
    n = 0;
    while abs(a-b)> tol
        n=n+1;
        xi = (a+b)/2;
        xit = [xit,xi];
        if func(a)*func(xi)<0
            b = xi;
        else
            a = xi;
        end
        
    end
    root = (a+b)/2;
end

