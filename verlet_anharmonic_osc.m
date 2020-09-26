%% Verlet algorithm (anharmonic driven oscillator)
clf; 

k = 1; m=1; gamma=.95; %set conditions on k and m
dt = .001; h = dt; %set size of time step
x = [];  a = []; 
x(1) = 1; x(2) = 1; a(1) = (-k/m)*x(1)- (gamma/m)*x(1)^3; %initial values
N = 1:10^5/2; %number of time steps

for i = 2:length(N)-1; %iterate through verlet algorithm to build x(t)
    a(i) = (-k/m)*x(i) + (gamma/m)*x(i)^3; %added forcing term
    x(i+1) = 2*x(i) - x(i-1) + (h^2)*a(i);
end
E = 0.5*k*x.^2 + (gamma/4)*x.^4; %total energy 
subplot(2,1,1); plot(N*dt,x); xlabel('time'); ylabel('x(t)');
title('harmonic oscillator with driving force gamma*x^3');
subplot(2,1,2); plot(N*dt,E); xlabel('time'); ylabel('total energy');
title('total energy of driven harmonic oscillator');