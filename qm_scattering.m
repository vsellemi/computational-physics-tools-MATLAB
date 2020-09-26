%% 04/27, Victor Sellemi

%% QM scattering transmission probability through 1eV high .5nm wide barrier

h = 1; m = 1; %constants
d = 9.4486299429; %conversion of .5nm to bohrs
e = 0.0367493; %conversion from eV to hartree
z = [0,d]; %position of each layer in bohrs
V0 = 1*e; %potential in barrier in hartree units
E = (0:0.001:5).*e; %energy values of particles in hartree units
kL1 = sqrt(2*m*E/h^2); %free space wavenumber
kR1 = sqrt(2*m*(E-V0)/h^2); %wavenumber in the barrier
nL1 = kL1; nR1 = kR1; %analogous n-values

t=zeros(length(E),1); %initialize transmission array
for i = 1:length(E); %loop over different values for wavelength    
    %set wavenumbers for corresponding wavelength
    kL = kL1(i); kR = kR1(i); 
    nL = nL1(i); nR = nR1(i); 
    M = eye(2); 
    for j = 1:length(z);
        if mod(j,2) == 1; k1 = kL; k2 = kR; n1 = nL; n2 = nR;
        elseif mod(j,2) == 0; k1 = kR; k2 = kL; n1 = nR; n2 = nL; end 
        T = [(0.5 + n2/(2*n1))*exp(1i*(k2 - k1)*z(j)),...
        (0.5 - n2/(2*n1))*exp(-1i*(k2 + k1)*z(j)); ...
        (0.5 - n2/(2*n1))*exp(1i*(k2 + k1)*z(j)), ...
        (0.5 + n2/(2*n1))*exp(-1i*(k2 - k1)*z(j))];
    M = M * T;  
    end
    t(i) = (abs(1/M(1,1))).^2; %calculate transmission probability
end

%plot results
plot(E./e,t,'-'); xlabel('Energy [eV]'); ylabel('Transmission'); axis([0 5 0 1]);
title('QM scattering through a single 1eV high .5nm wide barrier');

%%  QM scattering transmission probability through two 1eV high .5nm wide barriers
clear all; close all; 

h = 1; m = 1; %constants
d = 9.4486299429; %conversion of .5nm to bohrs
z = [0,d,3*d,4*d]; %position of barriers
e = 0.0367493; %conversion from eV to hartree
V0 = 1*e; %potential in barrier
E = (0:0.001:5).*e; %energy values of particles
kL1 = sqrt(2*m.*E/h^2); %free space wavenumber
kR1 = sqrt(2*m.*(E-V0)/h^2); %wavenumber in the barrier
nL1 = kL1; nR1 = kR1; 

t=zeros(length(E),1); %initialize transmission array
for i = 1:length(E); %loop over different values for wavelength    
    %set wavenumbers for corresponding wavelength
    kL = kL1(i); kR = kR1(i); 
    nL = nL1(i); nR = nR1(i); 
    %build transfer matrices
    M = eye(2); 
    for j = 1:length(z);
        if mod(j,2) == 1; k1 = kL; k2 = kR; n1 = nL; n2 = nR;
        elseif mod(j,2) == 0; k1 = kR; k2 = kL; n1 = nR; n2 = nL; end 
        T = [(0.5 + n2/(2*n1))*exp(1i*(k2 - k1)*z(j)),...
        (0.5 - n2/(2*n1))*exp(-1i*(k2 + k1)*z(j)); ...
        (0.5 - n2/(2*n1))*exp(1i*(k2 + k1)*z(j)), ...
        (0.5 + n2/(2*n1))*exp(-1i*(k2 - k1)*z(j))];
    M = M * T;  
    end
    t(i) = (abs(1/M(1,1))).^2; %calculate transmission probability
end

%plot results
plot(E./e,t,'-'); xlabel('Energy [eV]'); ylabel('Transmission'); axis([0 5 0 1]);
title('QM scattering through two 1eV high .5nm wide barriers with 1nm separation');

%% QM scattering through three, five, ten, and twenty barriers
clear all; close all; 

Ns = [3,5,10,20]; %number of barriers 
for b = 1:length(Ns);
N = Ns(b);
h = 1; m = 1; %constants
d = 9.4486299429; %conversion of .5nm to bohrs
pos = d.*cumsum(repmat([1,2],1,N)); 
z = [0, pos(1:2*N - 1)]; 
e = 0.0367493; %conversion from eV to hartree
V0 = 1*e; %potential in barrier
E = (0:0.001:5).*e; %energy values of particles
kL1 = sqrt(2*m.*E/h^2); %free space wavenumber
kR1 = sqrt(2*m.*(E-V0)/h^2); %wavenumber in the barrier
nL1 = kL1; nR1 = kR1; 

t=zeros(length(E),1); %initialize transmission array
for i = 1:length(E); %loop over different values for wavelength    
    %set wavenumbers for corresponding wavelength
    kL = kL1(i); kR = kR1(i); 
    nL = nL1(i); nR = nR1(i); 
    %build transfer matrices
    M = eye(2); 
    for j = 1:length(z);
        if mod(j,2) == 1; k1 = kL; k2 = kR; n1 = nL; n2 = nR;
        elseif mod(j,2) == 0; k1 = kR; k2 = kL; n1 = nR; n2 = nL; end 
        T = [(0.5 + n2/(2*n1))*exp(1i*(k2 - k1)*z(j)),...
        (0.5 - n2/(2*n1))*exp(-1i*(k2 + k1)*z(j)); ...
        (0.5 - n2/(2*n1))*exp(1i*(k2 + k1)*z(j)), ...
        (0.5 + n2/(2*n1))*exp(-1i*(k2 - k1)*z(j))];
    M = M * T;  
    end
    t(i) = (abs(1/M(1,1))).^2; %calculate transmission probability
end

%plot results
figure(b); plot(E./e,t,'-'); xlabel('Energy [eV]'); ylabel('Transmission'); axis([0 5 0 1]);
title(['QM scattering through ', num2str(N),' 1eV high .5nm wide barriers with 1nm separation']);

end

%% transmission probability through the double barrier potential for several values of the barrier width
clear all; close all; 

width = .2.*(1:5);
for b = 1:length(width);
w = width(b);
h = 1; m = 1; %constants
d = 9.4486299429; %conversion of .5nm to bohrs
pos = d.*cumsum(repmat([w,2],1,2)); 
z = [0, pos(1:3)]; %position of barriers
e = 0.0367493; %conversion from eV to hartree
V0 = 1*e; %potential in barrier
E = (0:0.001:5).*e; %energy values of particles
kL1 = sqrt(2*m.*E/h^2); %free space wavenumber
kR1 = sqrt(2*m.*(E-V0)/h^2); %wavenumber in the barrier
nL1 = kL1; nR1 = kR1; 

t=zeros(length(E),1); %initialize reflectance array
for i = 1:length(E); %loop over different values for wavelength    
    %set wavenumbers for corresponding wavelength
    kL = kL1(i); kR = kR1(i); 
    nL = nL1(i); nR = nR1(i); 
    %build transfer matrices
    M = eye(2); 
    for j = 1:length(z);
        if mod(j,2) == 1; k1 = kL; k2 = kR; n1 = nL; n2 = nR;
        elseif mod(j,2) == 0; k1 = kR; k2 = kL; n1 = nR; n2 = nL; end 
        T = [(0.5 + n2/(2*n1))*exp(1i*(k2 - k1)*z(j)),...
        (0.5 - n2/(2*n1))*exp(-1i*(k2 + k1)*z(j)); ...
        (0.5 - n2/(2*n1))*exp(1i*(k2 + k1)*z(j)), ...
        (0.5 + n2/(2*n1))*exp(-1i*(k2 - k1)*z(j))];
    M = M * T;  
    end
    t(i) = (abs(1/M(1,1))).^2; %calculate transmission probability
end

%plot results
plot(E./e,t,'-'); hold on; 
xlabel('Energy [eV]'); ylabel('Transmission'); axis([0 5 0 1]);
title('Transmission through double barrier for several values of barrier width');
if b == length(width); legend('.1nm','.2nm', '.3nm', '.4nm', '.5nm'); end

end

%%  transmission through 2-barrier potential using finite differences
clear all; close all; 

N = 50; dz = 1; h = 1; m = 1; A = (h^2/2/m/dz^2);

d = 9.4486299429; %conversion of .5nm to bohrs
z = [0,.5*d,2.5*d,3*d]; %position of barriers
e = 0.0367493; %conversion from eV to hartree
V0 = 1*e; %potential in barrier
E = (.001:0.001:5).*e; %energy values of particles
kL1 = sqrt(2*m.*E/h^2); %free space wavenumber
kR1 = sqrt(2*m.*(E-V0)/h^2); %wavenumber in the barrier
nL1 = kL1; nR1 = kR1; 

%second derivative finite differences operator in 1D
D2z=1/dz^2*(spdiags(-2*ones(N,1),0,sparse(N,N))+spdiags(ones(N,1),1,sparse(N,N))+spdiags(ones(N,1),-1,sparse(N,N)));
%potential energy operator
Vz = [V0*ones(1,.2*N), zeros(1,.6*N), V0*ones(1,.2*N)];
V = sparse(diag(Vz));

%%FINITE DIFFERENCES METHOD%%
T_fd = zeros(length(E),1); 
for i = 1:length(E); 
    kL = kL1(i); kR = kR1(i); 
    H = E(i)*speye(N) + A*D2z - V;
    H(1,1) = H(1,1) + (A)*exp(1i*kL*dz);
    H(end,end) = H(end,end) + (A)*exp(1i*kR*dz);
    b = zeros(N,1); b(1,1) = (1i*2*A)*sin(kL*dz); 
    psi = H\b;
    T_fd(i) = abs(psi(end,1))^2 * (sin(kR*dz)/sin(kL*dz));   
end


%%TRANSFER MATRICES METHOD%%
t=zeros(length(E),1); %initialize reflectance array
for i = 1:length(E); %loop over different values for wavelength    
    %set wavenumbers for corresponding wavelength
    kL = kL1(i); kR = kR1(i); 
    nL = nL1(i); nR = nR1(i); 
    %build transfer matrices
    M = eye(2); 
    for j = 1:length(z);
        if mod(j,2) == 1; k1 = kL; k2 = kR; n1 = nL; n2 = nR;
        elseif mod(j,2) == 0; k1 = kR; k2 = kL; n1 = nR; n2 = nL; end 
        T = [(0.5 + n2/(2*n1))*exp(1i*(k2 - k1)*z(j)),...
        (0.5 - n2/(2*n1))*exp(-1i*(k2 + k1)*z(j)); ...
        (0.5 - n2/(2*n1))*exp(1i*(k2 + k1)*z(j)), ...
        (0.5 + n2/(2*n1))*exp(-1i*(k2 - k1)*z(j))];
    M = M * T;  
    end
    t(i) = (abs(1/M(1,1))).^2; %calculate transmission probability
end


%plot results
plot(E(1:50:end)./e,t(1:50:end),'o',E./e,t); axis([0 5 0 1]);
xlabel('Energy[eV]'); ylabel('Transmission'); 
title('transmission probability through two barrier potential using FD'); 
legend('finite differences','transfer matrices','Location','southeast'); 


%% 1-D bandstructure, V(x)=Vcos(2pix/a), for V=0.5 eV. 
clear all; close all;

h = 1; m = 1; N = 100; A = h^2/2/m;

d = 9.4486299429; %conversion of .5nm to bohrs
e = 0.0367493; %conversion from eV to hartree
a = 2*d; 
x = linspace(0,a,N); 
V0 = 0.5*e; 
Vx = V0*cos((2*pi/a).*x);
k = pi/a;
G = [-2*pi/a, 0, 2*pi/a];
H = [A*(k+G(1))^2+V0, V0/2, V0/2; ...
    V0/2, A*(k+G(2))^2+V0,  V0/2; ...
    V0/2, V0/2, A*(k+G(3))^2 + V0];
[F,D] = eig(H);
%E = A*(k+G).^2;
Cs1 = F(:,1)'; Cs2 = F(:,2)'; 
fs1 = zeros(length(G),length(x)); fs2 = zeros(length(G),length(x)); 
for i = 1:length(G); 
    fs1(i,:) = exp(1i*G(i).*x);  
    fs2(i,:) = exp(1i*G(i).*x);
end
psi1 = Cs1*fs1;
psi2 = Cs2*fs2; 

%plot results
figure(1);
[AX, Y1, Y2] = plotyy(x./a,[abs(psi1).^2; abs(psi2).^2],x./a,Vx./e);
ylim(AX(2),[-1,1]); ylabel(AX(1),'|psi|^2'); ylabel(AX(2),'V(x) [eV]'); xlabel('position[a]'); 
title('Real-space eigenfunctions of the two lowest bands at k = ?/a'); 


%bandstructure
N = 7; 
a = 2*d; %a = 1nm
V0 = 0*e; %V0 = .5eV 
x = (-(N-1)/2:1:(N-1)/2)*a/pi;
Vx = V0*cos((2*pi/a).*x);
kk = linspace(-pi/a,pi/a,100); %k values from -pi/a to pi/a
GG = (-N+1:2:N-1).*(pi/a); %-6*pi/a:2*pi/a:6*pi/a; %G values 
EE = []; %initialize energy array

for i = 1:length(kk);
    
k=kk(i);
H = diag(A*(k+GG).^2+Vx(4));
for ii = 1:length(Vx); 
    VV = Vx([1,2,3,5,6,7]);
    H(ii,ii+1:end) = VV(1:end-ii+1);
    if ii > 1; H(ii,1:ii-1) = VV(end-ii+2:end); end;

end
[F,D] = eig(H); 
EE(i,:) = diag(D); 

end

figure(2); plot(kk*a/pi,EE./e,'k'); xlabel('k [pi/a]'); ylabel('Energy [eV]');
title('bandstructure at V0=0'); 

%%  scaling of the lowest three bandgaps with 0 < V ? < 1eV for the potential above
clear all; close all; 


h = 1; m = 1; N = 7; A = h^2/2/m;

e = 0.0367493;
V0s = (0:0.01:1).*e;
gap1 = []; gap2 = []; gap3 = []; 

for ee = 1:length(V0s);
    
V0 = V0s(ee);
d = 9.4486299429;
a = 2*d; %a = 1nm
x = (-(N-1)/2:1:(N-1)/2)*a/pi;
Vx = V0*cos((2*pi/a).*x);
kk = linspace(-pi/a,pi/a,100); %k values from -pi/a to pi/a
GG = (-N+1:2:N-1).*(pi/a); %-6*pi/a:2*pi/a:6*pi/a; %G values 
EE = []; %initialize energy array
for i = 1:length(kk);
    
k=kk(i);
H = diag(A*(k+GG).^2+Vx(4));
for ii = 1:length(Vx); 
    VV = Vx([1,2,3,5,6,7]);
    H(ii,ii+1:end) = VV(1:end-ii+1);
    if ii > 1; H(ii,1:ii-1) = VV(end-ii+2:end); end;

end
[F,D] = eig(H); 
EE(i,:) = diag(D); 

end

gap3(ee) = (abs(EE(end,1) - EE(end,2))'); 
gap1(ee) = (abs(EE(length(kk)/2,4) - EE(length(kk)/2,5))');
gap2(ee) = (abs(EE(end,3) - EE(end,4))'); 

end; 

plot(V0s./e,-fliplr(gap1./e),V0s./e,-fliplr(gap2./e),V0s./e,-fliplr(gap3./e)); 
legend('gap1','gap2','gap3');
ylabel('gap[eV]'); xlabel('V0[eV]');
title('Gap scaling of lowest three bandgaps with 0<V0<1'); 

%%  comparison of 1-D bandstructure of an infinite lattice and
%transmission coefficient for twenty 0.5eV-high barriers separated by 1nm
clear all; close all; 

N = 20;
h = 1; m = 1; %constants
d = 9.4486299429; %conversion of .5nm to bohrs
pos = d.*cumsum(repmat([1,2],1,N)); 
z = [0, pos(1:2*N - 1)]; 
e = 0.0367493; %conversion from eV to hartree
V0 = 1*e; %potential in barrier
E = (0:0.001:5).*e; %energy values of particles
kL1 = sqrt(2*m.*E/h^2); %free space wavenumber
kR1 = sqrt(2*m.*(E-V0)/h^2); %wavenumber in the barrier
nL1 = kL1; nR1 = kR1; 

t=zeros(length(E),1); %initialize transmission array
for i = 1:length(E); %loop over different values for wavelength    
    %set wavenumbers for corresponding wavelength
    kL = kL1(i); kR = kR1(i); 
    nL = nL1(i); nR = nR1(i); 
    %build transfer matrices
    M = eye(2); 
    for j = 1:length(z);
        if mod(j,2) == 1; k1 = kL; k2 = kR; n1 = nL; n2 = nR;
        elseif mod(j,2) == 0; k1 = kR; k2 = kL; n1 = nR; n2 = nL; end 
        T = [(0.5 + n2/(2*n1))*exp(1i*(k2 - k1)*z(j)),...
        (0.5 - n2/(2*n1))*exp(-1i*(k2 + k1)*z(j)); ...
        (0.5 - n2/(2*n1))*exp(1i*(k2 + k1)*z(j)), ...
        (0.5 + n2/(2*n1))*exp(-1i*(k2 - k1)*z(j))];
    M = M * T;  
    end
    t(i) = (abs(1/M(1,1))).^2; %calculate transmission probability
end


%add comparison with scattering
h = 1; m = 1; %constants
d = 9.4486299429;
e = 0.0367493;
A = h^2/2/m;
N = 7; 
a = 2*d; %a = 1nm
V0 = 1*e; %V0 = .5eV 
x = (-(N-1)/2:1:(N-1)/2)*a/pi;
Vx = V0*cos((2*pi/a).*x);
kk = linspace(-pi/a,pi/a,100); %k values from -pi/a to pi/a
GG = (-N+1:2:N-1).*(pi/a); %-6*pi/a:2*pi/a:6*pi/a; %G values 
EE = []; %initialize energy array

for w = 1:length(kk);
    
k=kk(w);
H = diag(A*(k+GG).^2+Vx(4));
for ii = 1:length(Vx); 
    VV = Vx([1,2,3,5,6,7]);
    H(ii,ii+1:end) = VV(1:end-ii+1);
    if ii > 1; H(ii,1:ii-1) = VV(end-ii+2:end); end;
end
[F,D] = eig(H); 
EE(w,:) = diag(D); 

end

figure(1);[AX, Y1, Y2] = plotyy(EE./5/e,kk*a/pi,E./e,t);
ylabel(AX(1),'k[pi/a]'); ylabel(AX(2),'Transmission'); xlabel('Energy [eV]'); 
title('Comparison with scattering of 20 barrier potential'); 

