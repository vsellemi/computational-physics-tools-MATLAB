%% 04/06, Victor Sellemi

%% the bandstructure of an infinite stack of dielectric slabs
%perpendicular propagation to layers

a = 100; %lattice constant
dz = 1; N = 100; c = 1; 

%forward difference and backward difference d/dz operators
D1f = (1/dz)*(spdiags(-ones(N,1),0,sparse(N,N))+spdiags(ones(N,1),-1,sparse(N,N)));
D1b = (1/dz)*(spdiags(ones(N,1),0,sparse(N,N))+spdiags(-ones(N,1),1,sparse(N,N)));

%Wave vectors in the z and y directions
Ky = linspace(0,1/(2*pi),N/5);
Kz = linspace(0,1/(2*pi),N/5);

e = [ones(1,.4*a) 13*ones(1,.2*a) ones(1,.4*a)]; %epsilon values
Ez = sparse(diag(e)); %E(z) operator


%M = -(D1b + 1i*k*speye(N))*(Ez\(D1f+1i*k*speye(N)));

Es = []; E = zeros(3,length(Kz));

for w = 1:length(Kz);
    for j = 1:length(Ky);
        ky1 = Ky(j); kz1 = Kz(w); 
        M1 = (D1b + 1i*kz1*speye(N))*inv(Ez)*(D1f+1i*kz1*speye(N));
        M1(1,1) = M(1,1) - 1/dz^2; 
        sL1 = inv(Ez)*((ky1^2).*speye(N)) - M1;
        Es(:,j) = eigs(sL1,3,'SM');  
    end
   E(1:3,w) = Es(:,1);
end

%plot results

X = Kz.*(a/(2*pi));
Y = sqrt(Es).*(a/(2*pi*c));

plot(X,Y); 
xlabel('k_za/2pi'); ylabel('omega a/2pi c'); 
title('Band structure of infinite stack of dielectric slabs'); 

%% the temperature dependence of magnetization for a 25×25 2-dimensional array
close all;

%METROPOLIS ALGORITHM%

n = 25; %set size of nxn array
N = n^4; %number of iterations of metropolis algorithm
J = 1; %exchange constant

S = ones(n,n); %initialize ordered spin ups
M = []; %initialize magnetization array

T = 0:0.2:4; %temperature vector 

for t = 1:length(T); %loop over different temperature values    
for k = 1:N; %iterations of Metropolis algorithms
    i = randi(n); j = randi(n); %generate random indices
    %enforce periodic boundary conditions
    if i == 1; im1 = n; else im1 = i-1; end;
    if i == n; ip1 = 1; else ip1 = i+1; end;
    if j == 1; jm1 = n; else jm1 = j-1; end;
    if j == n; jp1 = 1; else jp1 = j+1; end;
    %calculate total energy from alignment of nearest neighbors
    E = -J*(S(ip1,j)+S(im1,j)+S(i,jp1)+S(i,jm1))*S(i,j); 
    %if E is positive flip spin, if E is negative then flip spin with
    %Boltzmann-Arhenius probability
    if E > 0 || E<=0 && rand(1,1) <= exp(2*E/T(t)); 
        S(i,j) = -1*S(i,j);
    end
   
end
M(t) = sum(sum(S))./(n*n); %calculate magnetization
end

%ANALYTIC ONSAGER RESULT%
f  = @(x) (1 - sinh(2*J./x).^(-4)).^(1/8);
M_analytic = f(T(2:end)); 

for i = 1:length(M_analytic); %set M=0 above Tc
    if isreal(M_analytic(i)) == 0; 
        M_analytic(i) = 0; 
    end
end

%PLOT RESULTS%
plot(T,M,T(2:end),M_analytic); axis([0 4 min(M) 1.1]); 
legend('Metropolis 25x25', 'Onsager 1944');
xlabel('Temperature [J/k]'); ylabel('Magnetization [mu N]');
title('Magnetization of a 25x25 array'); 

%% the magnetization at nonzero temperature for a 32×32 array 
close all; 

n = 32; %set size of nxn array
N = n^4; %number of iterations of metropolis algorithm
J = 1; %exchange constant
kT = 0.3;  %Boltzmann constant times Temperature

S = (4*pi).*rand(n,n) - 2*pi; %initialize random spin array
M = []; %initialize magnetization array

  
for k = 1:N; %iterations of Metropolis algorithms
    i = randi(n); j = randi(n); %generate random indices
    theta_r = .4*rand(1,1) - 0.2; % generate random rotation in (-.2,.2) 
    Sf = S(i,j) + theta_r; %new orientation of spin at (i,j) location
    %periodic boundary conditions
    if i == 1; im1 = n; else im1 = i-1; end;
    if i == n; ip1 = 1; else ip1 = i+1; end;
    if j == 1; jm1 = n; else jm1 = j-1; end;
    if j == n; jp1 = 1; else jp1 = j+1; end;
    %calculate change in energy from nearest neighbor interactions
    Ei = -J*(cos(S(i,j) - S(im1,j)) + cos(S(i,j) - S(ip1,j)) ...
        + cos(S(i,j) - S(i,jm1)) + cos(S(i,j) - S(i,jp1))); 
    Ef = -J*(cos(Sf - S(im1,j)) + cos(Sf - S(ip1,j)) ...
        + cos(Sf - S(i,jm1)) + cos(Sf - S(i,jp1))); 
    dE = Ef - Ei; %change in energy dE
    %if dE is negative apply the rotation and if dE is positive, apply the
    %rotation with probability P=exp(-dE/kT)
    if dE <= 0 || dE > 0 && rand(1,1) <= exp(-dE/kT); 
        S(i,j) = S(i,j) + theta_r;
    end
end

%plot results
[x,y] = meshgrid(1:n,1:n); X = cos(S); Y = sin(S); 
quiver(x,y,X,Y); axis([0 32 0 32]); title('orientations of a 32x32 spin lattice');
