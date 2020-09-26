%Victor Sellemi

%The function returns the Lagrange interpolating polynomial for the
%function f using the nodes specified by the values in input x1. The output
%polynomial is a Matlab function. Define input function f as an anonymous
%function for intended results

function Q = lagrangeinterp(f,x1)
syms x; %x is a symbolic variable for function building
y1 = f(x1); %evaluate the function at the nodes
n = length(x1); %n denotes the number of nodes
Q = 0; %initialize the interpolating polynomial

for k = 1:n; %loop over nodes to build the interpolating polynomial
    l = 1;
    for j = [1:k-1 k+1:n];
        l = l*(x - x1(j))/(x1(k) - x1(j));
    end  
    Q = Q+l*y1(k);
end

Q = matlabFunction(Q); %convert symbolic function to Matlab function

end
