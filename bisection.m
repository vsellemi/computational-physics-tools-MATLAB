function [y,iterations] = bisection(f,a,b,e)

%this function attempts to find the root using the bisection algorithm
%write the function in f = 0 form
%provide initial interval in the form [a,b]
%specify absolute error tolerance e

if f(a)*f(b) > 0 %check if the interval contains a root of f
    disp('there is no root in the specified interval');
else
    y = (a+b)/2; %calculate the midpoint
    error = abs(f(y)); %calculate the absolute error
    count = 0; %initialize iterations count
    while error > e; %check if the error is less than specified tolerance
        count = count + 1; %count number of required iterations
        if f(a)*f(y)<0 %build next subinterval based on sign of f at midpoint
            b=y;
        else
            a=y;
        end
        y = (a+b)/2; %set final estimate once error tolerance is achieved
        error = abs(f(y));
        iterations = count; %output the number of iterations required to achieve tolerance
    end
end
end
        
