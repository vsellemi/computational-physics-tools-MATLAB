%% 05/04, Victor Sellemi

%% pz-orbital bandstructure of graphene
clear all; close all;

d = 9.4486299429; %conversion of .5nm to bohrs
e = 0.0367493; %conversion from eV to hartree
Vpp = 3*e; %Vpp = 3eV
Ep = 0; %Ep = 0; 
a = 2*d; %a = 1 bohr

%define k along IBZ perimeter in reciprocal space
kx1 = zeros(1,100); ky1 = linspace(0,4*pi/3/a,100); 
kx2 = linspace(0,pi/sqrt(3)/a,100); ky2 = linspace(4*pi/3/a,pi/a,100);
kx3 = linspace(pi/sqrt(3)/a,0,100); ky3 = linspace(pi/a,0,100);
kx = [kx1, kx2, kx3]; ky = [ky1, ky2, ky3]; 

%define n values
n1 = [a/sqrt(3),0];
n2 = [-a/sqrt(3)/2,a/2];
n3 = [-a/sqrt(3)/2,-a/2]; 

%Find the energies at different values of k
E1 = [];E2 = [];
for i = 1:length(kx)
    k = [kx(i),ky(i)];
    fk = exp(1i*dot(k,n1)) + exp(1i*dot(k,n2)) + exp(1i*dot(k,n3)); 
    fk1 = conj(fk);
    H = [Ep, -Vpp*fk; -Vpp*fk1, Ep]; %2x2 matrix
    [~,D] = eig(H); 
    E1(i) = D(end,end); E2(i) = D(1,1);
end

%plot results
mx=[200,200]; my=[-10,10]; Kx=[100,100]; Ky=[-10,10]; X = 1:length(kx); 
figure(1); plot(X,E1./e,X,E2./e,mx,my,'--k',Kx,Ky,'--k'); ylabel('Energy [eV]');
set(gca,'xtick',[0,100,200,300]); set(gca,'xticklabel',{'G';'K';'M';'G'});
title('Pz-orbital bandstructure of graphene'); 

%% 20fs time evolution of an electron wavefunction in a quadratic potential 
% ground state, hw = 1eV
clear all; close all; 

d = 9.4486299429; %conversion of .5nm to bohrs

%set constants
h = 1; m = 1; L = 2*d;
dx = .1; 
N = 256; dt = .01; %dt  = 0.01fs
omega = (1/h);

x = (d/2).*linspace(-2,2,N); %position vector
V = 0.5*m*omega^2.*x.^2; %potential vector
G = (2*pi/L)*((-N/2+1):N/2); %KE operator in matlab ordered basis
T=fftshift(G.^2*h^2/2/m); T=[0, T(1:end-1)]; 

%find the eigenvectors of the Hamiltonian matrix
D2x = 1/dx^2*(spdiags(-2*ones(N,1),0,sparse(N,N))+ ...
    spdiags(ones(N,1),1,sparse(N,N))+spdiags(ones(N,1),-1,sparse(N,N)));
H = -D2x + sparse(diag(V)); [F,~] = eig(full(H)); 
psi0 = -F(:,1)'; %ground state

YY = []; %solution array for psi^*psi
time = 0;ii=1; %initialize for loop
while time < 20 %time evolution in fs
    psi=exp(-1i.*V*dt/2/h).*ifft(exp(-1i.*T*dt/h).*fft(exp(-1i.*V*dt/2/h).*psi0)); 
    psi0 = psi; YY(ii,:) = abs(psi).^2; %store results each timestep
    time = time + dt; ii = ii+1;
end

%plot results
figure(1); imagesc(x./(d/2),linspace(0,20,20/dt),YY); colormap('jet'); cc = colorbar; 
xlabel('position [nm]'); ylabel('time [fs]'); ylabel(cc, 'psi*psi'); 
title('20fs time evolution of an electron ground state in quadratic V'); 


%%  20fs time evolution of an electron wavefunction in a quadratic potential 
%Equal and symmetric superposition of ground state and first excited state
clear all; close all; 

d = 9.4486299429; %conversion of .5nm to bohrs

%set constants
h = 1; m = 1; dx = .1; L = 2*d; 
N = 256; dt = 0.01; %dt  = 0.01fs
omega = (1/h);

x = (3*d/4)*linspace(-2,2,N); %position vector

V = 0.5*m*omega^2.*x.^2; %potential vector
G = (2*pi/L)*((-N/2+1):N/2); %KE operator in matlab ordered basis
T=fftshift(G.^2*h^2/2/m); T=[0, T(1:end-1)]; 

%find the eigenvectors of the Hamiltonian matrix
D2x = 1/dx^2*(spdiags(-2*ones(N,1),0,sparse(N,N))+ ...
    spdiags(ones(N,1),1,sparse(N,N))+spdiags(ones(N,1),-1,sparse(N,N)));
H = -D2x + sparse(diag(V)); [F,D] = eig(full(H)); 
psi00 = F(:,1)'; %ground state
psi11 = F(:,2)'; %first excited state

psi0 = (psi00 + psi11)./sqrt(2); %superposition of states

YY = []; %solution array for psi^*psi
time = 0;ii=1; %initialize for loop
while time < 20; %propagate in time
    psi = exp((-1i*dt/2/h)*V).*ifft(exp((-1i*dt/h)*T).*fft(exp((-1i*dt/2/h)*V).*psi0));
    psi0 = psi; 
    YY(ii,:) = abs(psi).^2;
    time = time + dt; ii = ii+1;
end

%plot results
figure(1); imagesc(x./(3*d/4),linspace(0,20,20/dt),YY); colormap('jet'); cc = colorbar; 
xlabel('position [nm]'); ylabel('time [fs]'); ylabel(cc, 'psi*psi'); 
title('Superposition of electron ground and first excited state, 20fs time evolution in quadratic V'); 



%% 20fs time evolution of an electron wavefunction in a quadratic potential 
%ground state, including self-consistent repulsion
clear all; close all; 

d = 9.4486299429; %conversion of .5nm to bohrs

%set constants
h = 1; m = 1; L = 2*d;
dx = .1; 
N = 256; dt = .01; %dt  = 0.01fs
omega = (1/h);

x = (d/2).*linspace(-2,2,N); %position vector
V0 = 0.5*m*omega^2.*x.^2; %potential vector
G = (2*pi/L)*((-N/2+1):N/2); %KE operator in matlab ordered basis
T = fftshift(G.^2*h^2/2/m); T=[0, T(1:end-1)]; 

%find the eigenvectors of the Hamiltonian matrix
D2x = 1/dx^2*(spdiags(-2*ones(N,1),0,sparse(N,N))+ ...
    spdiags(ones(N,1),1,sparse(N,N))+spdiags(ones(N,1),-1,sparse(N,N)));
H = -D2x + sparse(diag(V0)); [F,~] = eig(full(H)); 
psi0 = -F(:,1)'; %ground state

YY = []; %solution array for psi^*psi
time = 0;ii=1; %initialize for loop
while time < 20; %time evolution in fs
    V = (N/d)*abs(psi0).^2 + V0; 
    psi=exp(-1i.*V*dt/2/h).*ifft(exp(-1i.*T*dt/h).*fft(exp(-1i.*V*dt/2/h).*psi0)); 
    psi0 = psi; YY(ii,:) = abs(psi).^2; %store results each timestep
    time = time + dt; ii = ii+1;
end

%plot results
figure(1); imagesc(x./(d/2),linspace(0,20,20/dt),YY); colormap('jet'); cc = colorbar; 
xlabel('position [nm]'); ylabel('time [fs]'); ylabel(cc, 'psi*psi'); 
title('20fs time evolution of electron ground state, Gross-Pitaevsky eqn'); 

%% 20fs time evolution of an electron wavefunction in a quadratic potential 
%ground state subject to sudden additional inclusion of a constant Coulomb force
clear all; close all; 

d = 9.4486299429; %conversion of .5nm to bohrs

%set constants
h = 1; m = 1; L = 2*d;
dx = .1; 
N = 256; dt = .01; %dt  = 0.01fs
omega = (1/h);

x = (d/2).*linspace(-2,2,N); %position vector
V0 = 0.5*m*omega^2.*x.^2; %potential vector
G = (2*pi/L)*((-N/2+1):N/2); %KE operator in matlab ordered basis
T=fftshift(G.^2*h^2/2/m); T=[0, T(1:end-1)]; 

%find the eigenvectors of the Hamiltonian matrix
D2x = 1/dx^2*(spdiags(-2*ones(N,1),0,sparse(N,N))+ ...
    spdiags(ones(N,1),1,sparse(N,N))+spdiags(ones(N,1),-1,sparse(N,N)));
H = -D2x + sparse(diag(V0)); [F,~] = eig(full(H)); 
psi0 = -F(:,1)'; %ground state

YY = []; %solution array for psi^*psi
time = 0;ii=1; %initialize for loop
while time < 20; %time evolution in fs
    if time < 3; V = V0; elseif time >= 3; V = V0 + (3).*x; end; 
    psi=exp(-1i.*V*dt/2/h).*ifft(exp(-1i.*T*dt/h).*fft(exp(-1i.*V*dt/2/h).*psi0)); 
    psi0 = psi; YY(ii,:) = abs(psi).^2; %store results each timestep
    time = time + dt; ii = ii+1;
end

%plot results
figure(1); imagesc(x./(d/2),linspace(0,20,20/dt),YY); colormap('jet'); cc = colorbar; 
xlabel('position [nm]'); ylabel('time [fs]'); ylabel(cc, 'psi*psi'); 
title('20fs time evolution of an electron ground state in quadratic V with Coulomb force after 3fs'); 

%%  Dynamic scattering of a wavepacket
clear all; close all; 

%set constants
h = 1; m = 1; L = 300;
N = 512; dt = .1; %dt  = 0.01fs
x0 = 20; %x0 = 20nm
s0 = 1; %s0 = 1nm
k = 6;%5; %k = 5nm^-1
V0 = 1; %V0 = 1eV 
sV = 1; %sV = 1nm

x = linspace(-40,40,N); %position vector
V = V0*exp((-(x.^2))/2/sV^2); 
G = (2*pi/L)*((-N/2+1):N/2); %KE operator in matlab ordered basis
T = fftshift(G.^2*h^2/2/m); T=[0, T(1:end-1)]; 
psi0 = exp(-(0.5/2/s0^2).*(x+x0).^2).*exp(1i*k*x); 

YY = []; %solution array for psi^*psi
time = 0;ii=1; %initialize for loop
while time < 100; %time evolution in fs
    psi=exp(-1i.*V*dt/2/h).*ifft(exp(-1i.*T*dt/h).*fft(exp(-1i.*V*dt/2/h).*psi0)); 
    psi0 = psi;
    YY(ii,:) = abs(psi).^2;
    time = time + dt; ii = ii+1;
end


%plot results
figure(1); imagesc(x,linspace(0,100,100/dt),YY); colormap('jet'); cc = colorbar; 
xlabel('position [nm]'); ylabel('time [fs]'); ylabel(cc, 'psi*psi'); 
title('Dynamic scattering with single barrier potential');


%%  Dynamic scattering
clear all; close all; 

%set constants
h = 1; m = 1; L = 300;
N = 256; dt = .1; %dt  = 0.01fs
x0 = 20; %x0 = 20nm
s0 = 1; %s0 = 1nm
k = 6; %k = 5nm^-1
V0 = 4; %V0 = 4eV 
sV = .1; %sV = .1nm
xV = 5; %xV = 5nm

x = linspace(-40,40,N); %position vector
V = V0*(exp((-(x.^2))/(2*sV^2)) + exp(-(x-xV).^2/2/sV^2)); 
G = (2*pi/L)*((-N/2+1):N/2); %KE operator in matlab ordered basis
T = fftshift(G.^2*h^2/2/m); T=[0, T(1:end-1)]; 
psi0 = exp(-(0.5/2/s0^2).*(x+x0).^2).*exp(1i*k*x); 

YY = []; %solution array for psi^*psi
time = 0;ii=1; %initialize for loop
while time < 100; %time evolution in fs
    psi=exp(-1i.*V*dt/2/h).*ifft(exp(-1i.*T*dt/h).*fft(exp(-1i.*V*dt/2/h).*psi0)); 
    psi0 = psi;
    YY(ii,:) = abs(psi).^2;
    time = time + dt; ii = ii+1;
end


%plot results
figure(1); imagesc(x,linspace(0,50,50/dt),YY); colormap('jet'); cc = colorbar; 
xlabel('position [nm]'); ylabel('time [fs]'); ylabel(cc, 'psi*psi'); 
title('Dynamic scattering with double barrier potential');  


