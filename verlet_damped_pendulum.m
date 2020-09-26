%% Verlet algorithm (damped and driven pendulum)
clf; %clear all previous figures

omega = 2; Q = 1; gamma = 0.1; K = 1; %set constant values
dt = 0.001; h=dt; %set time step
x1 = []; x=[]; v = []; a = []; %build arrays for x,v, and a
%in this case x = theta, v = angular velocity, and a = angular acceleration
x(1) = 1; x(2) = 1; x1(1) = 1; x1(2) = 1; %set initial values
v(1)= 0 ;a(1)=-K*sin(x(1))-gamma*v(1)+Q*sin(omega*1*dt); 
N = 1:10^5; %total number of time steps
for i = 2:length(N)-1; %use verlet twice to simulate damped driven pendulum
    a(i) = -K*sin(x(i)) - (gamma)*v(i-1) + Q*sin(omega*i*dt); 
    x1(i+1) = 2*x1(i) - x1(i-1) + (h^2)*a(i); 
    v(i+1) = (x1(i+1)-x(i-1))/(2*h); 
    x(i+1) = x(i) + h*v(i); 
end
subplot(2,1,1); plot(N*dt,x); xlabel('time'); ylabel('theta(t)');
title('angle of damped driven pendulum');
subplot(2,1,2); plot(x,v); xlabel('theta(t)'); ylabel('w(t)')
title('phase-space of damped driven pendulum');