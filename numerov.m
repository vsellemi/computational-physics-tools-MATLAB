function y = numerov(E)

h = 1; m = 1; dx = .001; h = dx; N = 200; %set parameters

W = zeros(1,N); psi = zeros(1,N);
V = 0; f = (2*m/h)*(V-E); 
psi(1) = 0; psi(2) = 0.1;
W(1) = (1-(h^2/12)*f)*psi(1); 
 
for n = 2:N-1; 
W(n) = (1-(h^2/12)*f)*psi(n); 
W(n+1) = h^2*f*psi(n) + 2*W(n) - W(n-1);
psi(n+1) = W(n+1)/(1-(h^2/12)*f);
end

y = psi(N)

end
