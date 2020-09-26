function [y,iterations] = regulafalsi(f,a,b,e)

%this function attempts to find the root using the regula falsi algorithm
%write the function in f = 0 form
%provide initial interval in the form [a,b]
%specify absolute error tolerance e

if f(a)*f(b) > 0 %check if the interval contains a root of f
    disp('there is no root in the specified interval');
else 
    y = (f(a)*b - f(b)*a)/(f(a)-f(b)); %implement regula falsi estimate of root
    error = abs(f(y)); %calculate absolute error
    count = 0; %initialize count of iterations
    while error > e; %check if error is below set tolerance
        count = count + 1;
        if f(a)*f(y)<0
            b=y;
        else 
            a=y;
        end
        y = (f(a)*b - f(b)*a)/(f(a)-f(b)); %output estimate for root
        error = abs(f(y));
        iterations = count;
    end
end
end
