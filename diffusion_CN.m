%% 03/16, Victor Sellemi

%% the time evolution of a gaussian distribution inside an ?insulated? 1-dimensional domain, subject to drift and diffusion using the Crank-Nicolson method

%set system parameters and initialize solution array
D = -1; v=1; dt = 1e-2; dx = 1e-1; Nx = 100; Nt = 1000; Y = []; 

%d2/dx2 operator
sLdx2=1/dx^2*(spdiags(-2*ones(Nx,1),0,sparse(Nx,Nx))+spdiags(ones(Nx,1),1,sparse(Nx,Nx))+spdiags(ones(Nx,1),-1,sparse(Nx,Nx)));
sLdx2(1,1) = -1/dx^2; sLdx2(end,end) = -1/dx^2; 

%d/dx operator
sLdx1 = spdiags(-1*ones(Nx,1),0,sparse(Nx,Nx))+spdiags(ones(Nx,1),-1,sparse(Nx,Nx));
sLdx1(end,end) = 0;

%full drift and diffusion operator
sL = D*sLdx2 + (v/dx)*sLdx1; 

Y(:,1) = normpdf(dx:dx:Nx*dx,2,1); %gaussian initial conditions
X = dx:dx:Nx*dx; %vector of positions

for n = 1:Nt-1; %loop over time steps
    %Crank Nicholson iteration scheme 
    Y(:,n+1) = ((inv(eye(Nx) - (0.5*D*dt*sL)))*(eye(Nx) + (0.5*D*dt*sL)))*Y(:,n);
    %plot results at every 50 timesteps and last 200 timesteps
    if rem(n,50) == 1 || rem(n,1000) > 800; 
    plot(X,Y(:,n)); xlabel('distance'); ylabel('density'); hold on;  
    title('time evolution of gaussian in 1d subject to drift and diffusion');
    end
end

%% The reflection spectrum from a 100 nm-thick soap membrane (n=1.33), using the transfer-matrix method applied to front and back interfaces with air 
close all;
n_soap = 1.33; n_air = 1; d = 100e-9; %set parameter values
z = [0,d]; %position of each layer

lambda = (380:1:750)*1e-9; %wavelength in nm
k_air1 = (2*pi*n_air)./lambda; %wave number values for air
k_soap1 = (2*pi*n_soap)./lambda; %wave number values for soap

r=[]; %initialize reflectance array
for i = 1:length(lambda); %loop over different values for wavelength    
    %set wavenumbers for corresponding wavelength
    k_air = k_air1(i); k_soap = k_soap1(i);     
    %build transfer matrices
    T_air_soap = [(0.5 + n_soap/(2*n_air))*exp(1i*(k_soap - k_air)*z(1)),...
        (0.5 - n_soap/(2*n_air))*exp(-1i*(k_soap + k_air)*z(1)); ...
        (0.5 - n_soap/(2*n_air))*exp(1i*(k_soap + k_air)*z(1)), ...
        (0.5 + n_soap/(2*n_air))*exp(-1i*(k_soap - k_air)*z(1))];
    T_soap_air = [(0.5 + n_air/(2*n_soap))*exp(1i*(k_air - k_soap)*z(2)),...
        (0.5 - n_air/(2*n_soap))*exp(-1i*(k_air + k_soap)*z(2)); ...
        (0.5 - n_air/(2*n_soap))*exp(1i*(k_air + k_soap)*z(2)),...
        (0.5 + n_air/(2*n_soap))*exp(-1i*(k_air - k_soap)*z(2))];
    M = T_air_soap * T_soap_air; %product of transfer matrices
    r(i) = (abs(M(2,1))/abs(M(1,1)))^2; %calculate reflectance
end

%plot results
plot(lambda.*10^9,r ,'-'); axis([300 800 .05 .08]); xlabel('wavelength(nm)');
ylabel('reflection'); title('reflection spectrum of a 100nm soap membrane');

%% The transmission spectrum from air into germanium, through a multilayer dielectric stack
close all;
n_Zns = 2.2; n_Ge = 4.2; n_air = 1; %set index of refraction values

d = [0; 97.33; 48.60; 761.47; 412.85; 720.06; 382.28; 705.03; 370.42;...
709.26; 358.23; 718.52; 353.08; 724.86; 360.01; 710.47; 398.52; 564.95;...
40.79; 224.72; 125.31; 133.58; 98.28; 268.21; 138.25; 238.01; 125.48; ...
232.65; 68.54; 168.55; 150.14; 254.28; 125.25; 307.19; 165.16; ...
256.22; 133.04; 289.60; 147.63; 266.04; 134.34; 265.60; 156.86; ...
294.15; 123.17; 250.12; 178.96; 528.64]*1e-9; %thickness of each layer

z=cumsum(fliplr(d')); %position of each layer

lambda = (2000:2:7000)*1e-9; %wavelength in nm
k_Zns1 = (2*pi*n_Zns)./lambda; %wavenumber for Zns
k_Ge1 = (2*pi*n_Ge)./lambda; %wavenumber for Ge
k_air1 = (2*pi*n_air)./lambda; %wavenumber for air

T = []; t = []; 

for i = 1:length(lambda);  
    k_Zns = k_Zns1(i); k_Ge = k_Ge1(i); k_air = k_air1(i); %set k
    for k = 1:length(z);
    %build transfer matrices
    T_air_Zns{k} = [(0.5 + n_Zns/(2*n_air))*exp(1i*(k_Zns - k_air)*z(k)),...
    (0.5 - n_Zns/(2*n_air))*exp(-1i*(k_Zns + k_air)*z(k)); ...
    (0.5 - n_Zns/(2*n_air))*exp(1i*(k_Zns + k_air)*z(k)), ...
    (0.5 + n_Zns/(2*n_air))*exp(-1i*(k_Zns - k_air)*z(k))];

    T_Zns_Ge{k} = [(0.5 + n_Ge/(2*n_Zns))*exp(1i*(k_Ge - k_Zns)*z(k)),...
    (0.5 - n_Ge/(2*n_Zns))*exp(-1i*(k_Ge + k_Zns)*z(k)); ...
    (0.5 - n_Ge/(2*n_Zns))*exp(1i*(k_Ge + k_Zns)*z(k)), ...
    (0.5 + n_Ge/(2*n_Zns))*exp(-1i*(k_Ge - k_Zns)*z(k))];

    T_Ge_Zns{k} = [(0.5 + n_Zns/(2*n_Ge))*exp(1i*(k_Zns - k_Ge)*z(k)),...
    (0.5 - n_Zns/(2*n_Ge))*exp(-1i*(k_Zns + k_Ge)*z(k)); ...
    (0.5 - n_Zns/(2*n_Ge))*exp(1i*(k_Zns + k_Ge)*z(k)), ...
    (0.5 + n_Zns/(2*n_Ge))*exp(-1i*(k_Zns - k_Ge)*z(k))];    
    end
    
    %product of transfer matrices in order
    M = T_air_Zns{1}*T_Zns_Ge{2}*T_Ge_Zns{3}*T_Zns_Ge{4}*T_Ge_Zns{5}*T_Zns_Ge{6} ...
        *T_Zns_Ge{7}*T_Ge_Zns{8}*T_Zns_Ge{9}*T_Ge_Zns{10}*T_Zns_Ge{11} ...
        *T_Ge_Zns{12}*T_Zns_Ge{13}*T_Ge_Zns{14}*T_Zns_Ge{15} ...
        *T_Ge_Zns{16}*T_Zns_Ge{17}*T_Ge_Zns{18}*T_Zns_Ge{19} ...
        *T_Ge_Zns{20}*T_Zns_Ge{21}*T_Ge_Zns{22}*T_Zns_Ge{23} ...
        *T_Ge_Zns{24}*T_Zns_Ge{25}*T_Ge_Zns{26}*T_Zns_Ge{27} ...
        *T_Ge_Zns{28}*T_Zns_Ge{29}*T_Ge_Zns{30}*T_Zns_Ge{31} ...
        *T_Ge_Zns{32}*T_Zns_Ge{33}*T_Ge_Zns{34}*T_Zns_Ge{35} ...
        *T_Ge_Zns{36}*T_Zns_Ge{37}*T_Ge_Zns{38}*T_Zns_Ge{39} ...
        *T_Ge_Zns{40}*T_Zns_Ge{41}*T_Ge_Zns{42}*T_Zns_Ge{43} ...
        *T_Ge_Zns{44}*T_Zns_Ge{45}*T_Ge_Zns{46}*T_Zns_Ge{47};
    t(i) = (n_Ge/n_air)*(abs(1/M(1,1)))^2; %calculate transmission
end

%plot results
plot(lambda.*10^9,t,'-'); xlabel('wavelength[nm]'); ylabel('Transmission');



%% The mode dispersion for a dielectric waveguide with ? = 11.4 and thickness a
close all;

a = 30; e_slab = 11.4; N = 200; dx=1; w = .1; c = 1; 
   
e = [ones(1,N/2-a/2), 11.4*ones(1,a), ones(1, N/2-a/2)]; %E(z) values
Ez = sparse(diag(e)); %E(z) operator
ky = sparse(diag((pi/(a))*ones(1,N))); %ky operator

%forward difference and backward difference d/dz operators
sLdz_fd = (1/dx)*(spdiags(-1*ones(N,1),0,sparse(N,N))+spdiags(ones(N,1),-1,sparse(N,N)));
sLdz_bd = (1/dx)*(spdiags(1*ones(N,1),0,sparse(N,N))+spdiags(-ones(N,1),1,sparse(N,N)));

%full master equation operator
sL = (ky.^2)*inv(Ez) - .5*sLdz_fd*inv(Ez)*sLdz_bd - ((w/c)^2)*speye(N);

[v,d] = eig(full(sL));

%plot results
figure(1); x = dx:dx:N;
plot(x,abs(v(:,1)).^2,x,abs(v(:,2)).^2,x,abs(v(:,5)).^2); 
legend('mode 1', 'mode 2', 'mode 5');

%plot dispersion relation
k_vector = (pi/(10*a):pi/(10*a):pi/a); 
for i = 1:length(k_vector);
    sL = (k_vector(i).^2).*inv(Ez) - sLdz_fd*inv(Ez)*sLdz_bd - ((w/c)^2)*speye(N);
    [a,b] = eig(full(sL));
    figure(2); plot(diag(b)); hold on; xlabel('ka/2pi'); 
    ylabel('omega a/2pi c');
end
