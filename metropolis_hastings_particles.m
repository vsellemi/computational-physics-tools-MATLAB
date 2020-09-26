%% 03/07, Victor Sellemi

%% position distribution of particles due to random diffusion

N = 10^4; %number of particles
Y = []; %initialize array of random positions
for k = 1:N; %loop for each particle
    %generate vector of N random positions after 1e2,1e3,1e4 timesteps
    Y(k,1) = sum(-1 + 2*round(rand(10^2,1))); 
    Y(k,2) = sum(-1 + 2*round(rand(10^3,1))); 
    Y(k,3) = sum(-1 + 2*round(rand(10^4,1)));
end

%plot results
X = -300:10:300;
Y1 = hist(Y(:,1),X); Y2 = hist(Y(:,2),X); Y3 = hist(Y(:,3),X); 
plot(X,Y1,X,Y2,X,Y3);
legend('100','1000','10000'); xlabel('position'); ylabel('density');
title('Density of 10^4 particles at various timesteps after random diffusion');

%% Metropolis-Hastings sampling histogram from exponential distribution
close all;

lambda = 1; Y = []; N = 10^4; %lamba = parameter, Y = array; N = number of points
f = @(x) lambda*exp(-lambda*x); %define desired distribution pdf
k = 1; %initialize index for result array

while length(Y) < N %Metropolis Hastings algorithm 
    y = 10*rand(1,1); %random y s.t. 0<y<10
    f1 = rand(1,1); %random f1 s.t. 0<f1<1
    if f(y) >= f1 %discard y if f(y) < f1, keep otherwise
        Y(k) = y; %store points that we want
        k = k+1; %increase index for final array of points
    end
end

%plot results
X = 0.1:0.1:10; X1 = 0:0.1:10;
Y = hist(Y,X); plot(X,Y,X1,f(X1).*(N*3/20));
legend('Metropolis-Hastings','Analytic'); xlabel('x'); ylabel('frequency');
title('Metropolis-Hastings algorithm applied to exponential distribution');

%%  several 2D ?random-walk? trajectories with decaying exponential jump-length distribution 
close all;

f = @(x) lambda*exp(-lambda*x); %exponential decay function 
X = 0.0001:0.0001:1; %timesteps
jumplength = f(X); %exponential decay in jump length over time

for i = 1:5; %generate 5 different possible paths
x = jumplength'.*(-1+2*rand(length(X),1)); %random x motion each time step 
y = jumplength'.*(-1+2*rand(length(X),1)); %random y motion each time step

%plot the displacement from origin for each path
plot(cumsum(x),cumsum(y)); hold on; %plot the displacement from origin
end

%%  the distribution of particle positions in 2D after 100 and 1000 timesteps
close all;

lambda = 1;
f = @(x) lambda*exp(-lambda*x); %exponential decay function 
X = 0.01:0.01:1; X1 = 0.001:0.001:1; %tf = 100 and tf1 = 1000
jumplength = abs(f(X)); jumplength1 = abs(f(X1)); %exponential decay in jump length
x = []; y = []; x1 = []; y1 = []; %initialize arrays
N = 10^3; %number of particles

%for N particles, we calculate the final random 2D displacement by taking
%the sum of 100 or 1000 random movements of magnitude [-6,6] in both the x
%and y directions. We choose the magnitude to be [-6,6] so that our figures
%most closely resemble those from lecture

for i = 1:N; %loops over each particle
x(i) = sum(jumplength'.*(-6+12*rand(length(X),1))); 
y(i) = sum(jumplength'.*(-6+12*rand(length(X),1)));
x1(i) = sum(jumplength1'.*(-6+12*rand(length(X1),1)));
y1(i) = sum(jumplength1'.*(-6+12*rand(length(X1),1)));
end

figure(1); plot(x,y,'o'); axis([-300 300 -300 300]);
figure(2); plot(x1,y1,'o'); axis([-300 300 -300 300]);

%% convergence of position distribution to gaussian, regardless of single-jump distribution in 1D
close all;

N = 10^4; %number of particles
f = @(x) exp(-0.2*(x-5).^2); 


Z = [];
w = [-1,1]; 
axislim = [20,30,40,50,50];

for i = 1:10
    Y = [];
    k=1;
    while length(Y) < N; %Metropolis Hastings algorithm 
        y = 20*rand(1,1); 
        f1 = rand(1,1);
        if f(y(1:length(y))) >= f1(1:length(y));
            index = randi([1,2],[1,length(y)]);
            ww = [];
            for j = 1:length(index);
                ww(j) = w(index(j));
            end
            Y(k) = sum(y.* ww);
            k = k+1; 
        end
    end
    Z(:,i) = Y;
end

figure(1); plot(-20:1:20,hist(Z(:,1),-20:1:20)./N);
W = []; 
Nt = [1,2,3,5,10];
for i = 2 : length(Nt);
    T = Nt(i);
    X = -axislim(i):1:axislim(i);
    W = sum(Z(:,1:T)');
    figure(i); plot(X,hist(W,X)./N);
end


%%  calculate the density at equilibrium inside an ?insulated? 1-dimensional domain
close all;

%initialize constants and arrays
D = 1; dx = .1; N = 100; Y = []; Z=[]; 
V = [0.1,0.5,1]; X = .1:.1:10; 

for i = 1:length(V); %loop over possible values of velocity
v = V(i);

%d2/dx2 operator
sLdx2=1/dx^2*(spdiags(-2*ones(N,1),0,sparse(N,N))+spdiags(ones(N,1),1,sparse(N,N))+spdiags(ones(N,1),-1,sparse(N,N)));
sLdx2(1,1) = -1/dx^2; sLdx2(end,end) = -1/dx^2; 
sLdx2=[ones(1,N);sLdx2];

%d/dx operator
sLdx1 = spdiags(-1*ones(N,1),0,sparse(N,N))+spdiags(ones(N,1),-1,sparse(N,N));
sLdx1(end,end) = 0;
sLdx1 = [zeros(1,N);sLdx1];

sL = D*sLdx2 + (v/dx)*sLdx1; %full operator

M = zeros(N,1); M = [1;M];

Y(:,i) = (sL\M)./dx; %left division

Z(:,i) = (v/(D*(exp(v*N*dx*(1/D))-1))).*exp((v/D).*X); %analytic sol. 

end

%plot results
plot(X,Y(:,1),'ob',X,Z(:,1),'-r',X,Y(:,2),'ob',X,Y(:,3),'ob',X,Z(:,2),'-r',X,Z(:,3),'-r'); 
legend('calc','analytic'); title('drift and diffusion at equilibrium');
xlabel('position'); ylabel('density');

%%  coupled drift, diffusion, and spin-flip with opposite drift velocity and particle conservation
close all;

D = 1; dx = 1e-6; N = 100; P = []; v = 10^3;
T = [1e-9,1e-10,1e-11];

for i = 1:length(T);
    
tau = T(i); 

%d2/dx2 operator
sLdx2=1/dx^2*(spdiags(-2*ones(N,1),0,sparse(N,N))+spdiags(ones(N,1),1,sparse(N,N))+spdiags(ones(N,1),-1,sparse(N,N)));
sLdx2(1,1) = -1/dx^2; sLdx2(end,end) = -1/dx^2; 

%d/dx operator for spin up
sLdx1_up = spdiags(-1*ones(N,1),0,sparse(N,N))+spdiags(ones(N,1),-1,sparse(N,N));
sLdx1_up(end,end) = 0;

%d/dx operator for spin down
sLdx1_down = spdiags(-1*ones(N,1),0,sparse(N,N))+spdiags(ones(N,1),1,sparse(N,N));
sLdx1_down(end,end) = 0;

%1/tau decay matrix
sLdecay = (1/tau)*speye(N);

%drift,diffusion,and -1/t on diagonals operator for spin up
sLD_up = D*sLdx2 + (v/dx)*sLdx1_up - sLdecay; 
sLD_down = D*sLdx2 - (v/dx)*sLdx1_down - sLdecay;

sL = [sLD_up sLdecay; sLdecay sLD_down]; %full matrix operator 
sL = [ones(1,2*N); sL]; %add conservation of particle number
M = zeros(2*N,1); M = [(2*N);M]; %RHS of matrix equation
Y = (sL\M); %solution by left division
n_up = Y(1:N); n_down = Y(N+1:end); %spin up and spin down
P(:,i) = (n_up - n_down)./(n_up+n_down); %polarization vector

end

%plot results
X = -(N*dx)/2+dx:dx:(N*dx)/2;
plot(X,P(:,1),'o',X,P(:,2),'o',X,P(:,3),'o'); 
legend('tau = 1e-9','tau = 1e-10','tau = 1e-11');
xlabel('x'); ylabel('polarization'); title('Polarization at equlibrium');


%%  coupled drift, diffusion, and relaxation with a source term for ?n
close all;

D = 1; dx = .1; N = 100; Y = []; %set parameters
v=2; T = [0.1,0.25,.5,1]; 

for i = 1:length(T);
tau = T(i);

%d2/dx2 operator
sLdx2=1/dx^2*(spdiags(-2*ones(N,1),0,sparse(N,N))+spdiags(ones(N,1),1,sparse(N,N))+spdiags(ones(N,1),-1,sparse(N,N)));
sLdx2(1,1) = -1/dx^2; sLdx2(end,end) = -1/dx^2; 

%d/dx operator
sLdx1 = spdiags(-1*ones(N,1),0,sparse(N,N))+spdiags(ones(N,1),-1,sparse(N,N));
sLdx1(end,end) = 0;

%1/tau decay matrix
sLdecay = (1/tau)*speye(N); 

%drift diffusion and decay operator
sL = D*sLdx2 + (v/dx)*sLdx1 - sLdecay; 

M = zeros(N,1); M(round(N/2)) = -1;

Y(:,i) = (sL\M); %backwards division to get numerical solution

end

X = .1:.1:10;
plot(X,Y(:,1),X,Y(:,2),X,Y(:,3),X,Y(:,4));
legend('.1','.25','.5','1');
title('coupled drift, diffusion,and relaxation with a source term');

%%  drift, diffusion, and relaxation with a source term and precession
close all;

D = 1; dx = .05; N = 100; Y = []; v=0; tau = 5; %set parameters
wx = 1; wy = 0; wz = 0;

%d2/dx2 operator
sLdx2=1/dx^2*(spdiags(-2*ones(N,1),0,sparse(N,N))+spdiags(ones(N,1),1,sparse(N,N))+spdiags(ones(N,1),-1,sparse(N,N)));
sLdx2(1,1) = -1/dx^2; sLdx2(end,end) = -1/dx^2; 

%d/dx operator
sLdx1 = spdiags(-1*ones(N,1),0,sparse(N,N))+spdiags(ones(N,1),-1,sparse(N,N));
sLdx1(end,end) = 0;

%1/tau decay matrix
sLdecay = (1/tau)*speye(N);

%omega times identity matrices
sL_wx = wx*speye(N); sL_wy = wy*speye(N); sL_wz = wz*speye(N);

%drift diffusion and decay matrix
H = D*sLdx2 + (v/dx)*sLdx1 - sLdecay;
Hx = D*sLdx2 + (1/dx)*sLdx1 - sLdecay;

%full matrix
sL = [Hx sL_wz -sL_wy; -sL_wz H sL_wx; sL_wy -sL_wx H];
sL = [ones(1,3*N);sL];

M = zeros(3*N,1); M=[3*N;M];

Y = (sL\M)./(N);

Sx = Y(1:N); 
Sy = Y(N+1:2*N);
Sz = Y(2*N+1:end);


%plot results
X = .05:.05:10;
plot(X,[Sx' fliplr(Sx')],X,[Sy' fliplr(Sy')],X,[Sz' fliplr(Sz')]);
xlabel('position'); ylabel('density');
