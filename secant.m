function [y,iterations] = secant(f,a,b,e)

%this function attempts to find the root using the secant method
%write the function in f = 0 form
%provide initial interval in the form [a,b]
%specify absolute error tolerance e

y = (f(a)*b - f(b)*a)/(f(a)-f(b)); %implement secant method iterations
error = abs(f(y)); %calculate absolute error
count = 0; %initialize count of iterations

while error > e
    %check if error is below set tolerance
    count = count + 1; 
    a = b;
    b = y;
    y = (f(a)*b - f(b)*a)/(f(a)-f(b)); %output estimate for root
    error = abs(f(y));
    iterations = count;
end
end
 
