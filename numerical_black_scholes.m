%% Victor Sellemi, 05/15/2018 

%% Crank-Nicholson Algorithm to solve the Black-Scholes equation %%
%We attempt to solve for the time evolution of the Black-Scholes equation
%using the Crank-Nicholson iteration. Here we solve for the time evolution
%of the value of a call option with strike price K
clear all; close all; 

dx = .01; %step size in x, x corresponds to stock price as S = Kexpx
dt = 0.001; %step size in time, corresponding to years
sigma = 0.1; %volatility, taken as an average over 1 year
r = 0.12; %spot interest rate, taken as a constant for 1 year
k = 2*r/sigma^2; %constant used in modified diffusion equation
Nx = 100; %number of x-steps
Nt = 10; %number of time steps, Nt*dt < 1, given the time scale

%Build finite differences operators for the Black-Scholes equation
D2x = 1/dx^2*(spdiags(-2*ones(Nx,1),0,sparse(Nx,Nx))+ ... 
    spdiags(ones(Nx,1),1,sparse(Nx,Nx))+ ...
    spdiags(ones(Nx,1),-1,sparse(Nx,Nx)));
D2x(1,1) = -1/dx^2; D2x(end,end) = -1/dx^2;  %include insulating boundary conditions

%First derivative must be specified in backwards difference form for k-1>0
D1x = -(1/dx)*spdiags(-ones(Nx,1),0,sparse(Nx,Nx))+ ... 
    spdiags(ones(Nx,1),-1,sparse(Nx,Nx));
D1x(end,end) = 0; 

HBS = D2x + (k-1)*D1x - k*speye(Nx); %full Black Scholes Hamiltonian operator
x = linspace(log(0.6),log(1.6),Nx); %vector of stock prices S = exp(x)
v = []; %initialize solution array

%specify the payoff function as the initial value for iterations, as the
%Black-Scholes diffusion equation transformation proceeds in tau time, with
%tau defined as tau = T-t
v(:,1) = max(exp(x) - 1,0); 


for n = 1:Nt-1 %loop over time steps
    %Crank-Nicholson explicit iteration scheme
    v(:,n+1) = ((inv(eye(Nx) - (0.5*dt*HBS)))*(eye(Nx) + (0.5*dt*HBS)))*v(:,n);
end

%plot results
figure(1); plot(100*exp(x(1:2:end)),Nx.*v(1:2:end,1:3:end),'o'); xlim([70 130]); hold on;
%set(gca,'fontsize',11); legend('t = T','t = 2','t = 1','t = 0', 'Location', 'northwest');  

%solve for the analytic solution%
K = 100; sigma = 0.1; r = 0.12; 
T = linspace(0,.4,10); 
S = linspace(70,130,100);
V = []; 
for n = 1:length(T); 
    t = T(n); 
    V(:,n) = S.*normcdf((log(S./K)+(r+sigma^2/2)*t)/sigma/sqrt(t)) - ...
        K*exp(-r*t).*normcdf((log(S./K)+(r-sigma^2/2)*t)/sigma/sqrt(t)); 
end
%add analytic results to the plot
plot(S,V(:,1:3:end),'-'); xlabel('Stock price S'); ylabel('Value of a call option C');



%% Hamiltonian time evolution of Black-Scholes equation
%here, we transform the Hamiltonian operator into "momentum space" to solve
%for the time evolution using a split exponentiated operator method,
%simulated again for a call option with strike price K
clear all; close all; 

%specify the same constants as above (keep in mind time is in years so time
%scale should be small
dx = .01; dt = .25; sigma = 0.1; r = 0.12; 
k = 2*r/sigma^2; Nx = 256; 

D2x = 1/dx^2*(spdiags(-2*ones(Nx,1),0,sparse(Nx,Nx))+ ... 
    spdiags(ones(Nx,1),1,sparse(Nx,Nx))+ ...
    spdiags(ones(Nx,1),-1,sparse(Nx,Nx)));
D2x(1,1) = -1/dx^2;
D2x(end,end) = -1/dx^2;

D1x = -(1/dx)*spdiags(-ones(Nx,1),0,sparse(Nx,Nx))+ ... 
    spdiags(ones(Nx,1),-1,sparse(Nx,Nx));
D1x(end,end) = 0;

HBS = -D2x + (k-1)*D1x - k*speye(Nx); %Black-Scholes Hamiltonian

N = 256; L = 3; 
G = (2*pi/L)*((-N/2+1):N/2); %Hamiltonian operator in matlab ordered basis

%build the Hamiltonian operator in momentum space
T = fftshift(G.^2*sigma^2/2 + (1i*0.5*sigma^2-r).*G); T=[0, T(1:end-1)]; 
V = r*ones(1,N); 
x = linspace(log(0.2),log(1.8),N); 
C0 = max(exp(x) - 1,0); %value of call option at maturity
g = @(x1) max(exp(x1) - 1,0); %payoff function 


v = []; 
time = 0;ii=1; %initialize for loop
while time < 1.75; %propagate in time
    v(:,ii) = C0; %update the output every iteration
    %plot results
    figure(1); plot(100*exp(x(1:2:end)),100.*C0(1:2:end),'o'); hold on; xlim([70 130]);
    %time evolution is specified by
    C = exp((dt/2).*V).*ifft(exp((-dt).*T).*fft(exp((dt/2).*V).*C0)); C(1:10) = 0;
    C0 = abs(C); %solution must be real
    time = time+dt; ii = ii + 1; 

end 
xlabel('Stock price S'); ylabel('Value of a call option C');

%Plot analytic results
K = 100; sigma = 0.1; r = 0.12; 
T = linspace(0,.4,10); 
S = linspace(70,130,100); 

V = []; 
for n = 1:length(T); 
    t = T(n); 
    V(:,n) = S.*normcdf((log(S./K)+(r+sigma^2/2)*t)/sigma/sqrt(t)) - ...
        K*exp(-r*t).*normcdf((log(S./K)+(r-sigma^2/2)*t)/sigma/sqrt(t)); 
end

plot(S,V(:,1:2:end),'-'); xlabel('Stock price S'); ylabel('Value of a call option C');

%plot in three dimensions
figure(2); h = surf(linspace(0,1,ii-1),100*exp(x),Nx.*v); %set(h,'LineStyle','none');
xlabel('time'); ylim([70 130]); zlim([0,100]); ylabel('X ; S = exp(X)'); zlabel('Price of a call option C'); 



%% Barrier options and potentials
clear all; close all; 

%specify constants
dt = .25; sigma = 0.1; r = 0.12; 
N = 256*3; L = 3; 
G = (2*pi/L)*((-N/2+1):N/2); %KE operator in matlab ordered basis
T = fftshift(G.^2*sigma^2/2 + (1i*0.5*sigma^2-r).*G); T=[0, T(1:end-1)]; 

%specify different potentials
VV = []; 
VV(:,1) = r*ones(1,N); 
VV(:,2) = [1*ones(1,N/2), 0*ones(1,N/2)];
VV(:,3) = [.5*ones(1,N/2), 0*ones(1,N/2)];
VV(:,4) = [0*ones(1,N/2), -.2*ones(1,N/2)];

ii = 1; 
while ii <= 4; 
    V = VV(:,ii)';

x = linspace(log(0.2),log(1.8),N); 
K = 1;
C0 = max(exp(x) - K,0); 
g = @(x1) max(exp(x1) - 1,0); 

time = 0; %initialize for loop
while time < 1.5 ; %propagate in time
    if ii == 1; figure(1); h1 = plot(100*exp(x(1:2:end)),100.*C0(1:2:end),'-b'); hold on; %xlim([85 110]);
    elseif ii == 2 && time > .5; figure(1); h2 = plot(100*exp(x(1:2:end)),100.*C(1:2:end),'-*r'); hold on; %xlim([80 120]);
    elseif ii == 3; figure(1); h3 = plot(100*exp(x(1:2:end)),100.*C(1:2:end),'-g'); hold on; xlim([80 120]);
    elseif ii == 4 && time > .5; figure(1); h4 = plot(100*exp(x(1:2:end)),100.*C(1:2:end),'--','Color',[.9 .6 0]); 
        hold on; xlim([90 110]);
    end
C = exp((dt/2).*V).*ifft(exp((-dt).*T).*fft(exp((dt/2).*V).*C0));
C(1:10) = 0; C0 = abs(C); C = abs(C);
time = time+dt; %ii = ii + 1; 
end 
xlabel('Stock price S'); ylabel('Value of a call option C');
ii = ii + 1; 
end

legend([h1,h2,h3,h4],'Black-Scholes','V = 1','V = 0.5','V = -0.2','Location','Northwest');



%% Analytic solution for the stochastic discount factor
clear all; close all; 

%set constants
sigma = .1; 
dx = .01; 
dt = .01; 
r = 0.12; 
Nx = 100;
Nt = 100;

T = (1:Nt).*dt; 
X = (-Nx/2+1:1:Nx/2).*dx; 
X1 = [-.3,0,.35]; 
for k = 1:length(X1); 
    x1 = X1(k); 
    pBS = zeros(Nx,Nt);
for n = 1:length(T); 
    for i = 1:length(X);
        t = T(n); x = X(i); 
        pBS(i,n) = (exp(-r*t)/sqrt(2*pi*t*sigma^2)) * ...
            exp(-(1/2/t/sigma^2)*(x-x1+t*(r-0.5*sigma^2))^2);
    end
end        
figure(1); plot(X,2.5.*pBS(:,1:20:Nx)./Nx); hold on; 
xlabel('X ; S = exp(X)'); ylabel('Pricing kernel');
end 

%% Monte-Carlo Solution for the stochastic discount factor
clear all; close all; 

p = []; %intialize solution array for the pricing kernel
XX = linspace(-.5,.5,5); %specify initial x array
for i = 1:length(XX) %iterate over different intial values of x
    
X0 = XX(i); 
mu = 0; %no drift term in Ito-Weiner evolution
var = .1^2; %variance
T = (0:.02:3); %time scale (years)
N = length(T); 

for n = 1:N %time evolution according to Ito-Weiner solution
    t = T(n); sigma = sqrt(var); 
    X(n) = X0 + normrnd(t*mu,sqrt(t)*sigma); 
end

%plot time evolution according to ito-weiner process
figure(1); plot(T,X); hold on; xlabel('Time (years)'); ylabel('X(t)');

edges = linspace(-1,1,75); 
times = [1,round(N/4),round(N/2),N]; %specify time values of interest 
X1 = X(time(1)); X2 = X(time(2)); X3 = X(time(3)); X4 = X(time(4)); 
%Quantum fluctuations (10^6 iterations)
X1_fluctuations = X1 + normrnd(0*mu,0*sigma,[1,1e6]);
X2_fluctuations = X2 + normrnd(T(25)*mu,sqrt(T(25))*sigma,[1,1e6]); 
X3_fluctuations = X3 + normrnd(T(75)*mu,sqrt(T(75))*sigma,[1,1e6]);
X4_fluctuations = X4 + normrnd(T(end)*mu,sqrt(T(end))*sigma,[1,1e6]);
%Calculate the probablity distribution
[N1,edges1] = histcounts(X1_fluctuations,edges); 
[N2,edges2] = histcounts(X2_fluctuations,edges); 
[N3,edges3] = histcounts(X3_fluctuations,edges); 
[N4,edges4] = histcounts(X4_fluctuations,edges); 
X = edges(2:end);

%plot normalized pricing kernel
figure(2); plot(X,N1./2e6,X,N2./1e6,X,N3./1e6,X,N4./1e6); hold on; 
 
end


%% Geometric Brownian Motion Example
clear all; close all; 

mu = 0; %drift 
var = [0,.1,.2,.3]; %variance vector
X0 = 1; 
T = (0:.03:1); 
for i = 1:length(var);
for n = 1:length(T); 
    t = T(n); sigma = sqrt(var(i)); 
    X(n) = X0 + normrnd(t*mu,sqrt(t)*sigma); 
end
plot(T,exp(X)); hold on; 
end
xlabel('time'); ylabel('S(t)'); 
title('time evolution of stock price subject to geometric brownian motion'); 
legend('sigma = 0', 'sigma = .1', 'sigma = .2', 'sigma = .3'); 