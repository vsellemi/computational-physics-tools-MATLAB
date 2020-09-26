%% 04/20, Victor Sellemi

%%  the ground-state wavefunction and eigenenergy for the hydrogen atom
%variational method with gaussian basis exp(?bnr?2), bn = {13, 2, 0.44, 0.12}

close all; clear all; 

b = [13,2,0.44,0.12]; %vector of bn values
N = length(b); %dimension of basis
r = 0:0.01:3; %r vector

phi = []; 
for i = 1:N; %construct the Gaussian basis with bn values
    bn = b(i); 
    phi(i,:) = exp(-bn*r.^2); 
end

S = zeros(N,N); T = zeros(N,N); V = zeros(N,N); 
for m = 1:N;
    for n = 1:N;
        S(m,n) = (pi/(b(m) + b(n)))^1.5; %overlap matrix
        T(m,n) = (3*b(m)*b(n)*pi^1.5)/(b(m)+b(n))^2.5; %KE operator
        V(m,n) = (-2*pi)/(b(m)+b(n)); %PE operator
    end
end

[F,D] = eig(T+V,S); %eigenvectors are coefficients of basis
Y = (F'*phi)'; Y = (1/27.2)*abs(Y(:,end)).^2; %build ground state wavefunction
E0 = D(end,end); %ground state energy

%plot results
YY = (8/3^3)*exp(-r./.5); %exact value of |psi|^2 with E = -0.5 Hartree
plot(r,Y,r,YY); xlabel('radial position'); ylabel('|psi|^2'); 
title('Ground state wavefunction of hydrogen atom'); 
legend(['variational groundstate, E = ', num2str(E0),' Hartree'], ...
    'exact, E = -0.5 Hartree');

%% the self-consistent single-electron wavefunction and binding energy of neutral helium

close all; clear all; 

%set constants and parameters
h = 1; m = 1; Z = 2; A  = 1.44; %A = e^2/4/pi/e0
N = 200; dr = 0.01; 
r = linspace(0.01,0.5,N);

%finite differences second derivative operator with respect to r
D2r=1/dr^2*(spdiags(-2*ones(N,1),0,sparse(N,N))+spdiags(ones(N,1),1,sparse(N,N))+spdiags(ones(N,1),-1,sparse(N,N)));
H1 = -(h^2/2/m)*D2r - sparse(diag((Z)*A*(1./r))); 

diffU = 1; count = 1; %initialize values
U_scf = zeros(1,N); %initial guess for U_scf

while diffU > 1e-8; %repeat until convergence
    U1 = U_scf;
    H = H1 + sparse(diag(U1)); %full Hamiltonian operator
    [F,D] = eig(full(H)); %caluclate eigenvectors and eigenvalues
    f = (F(:,1)'); %ground state eigenvector
    sigma = abs(f).^2; %self consistent field
    
    %plot results
    figure(1); plot(r,N.*f.*f); hold on; xlabel('distance[nm]');ylabel('Psi*Psi [nm^-1]'); 
    title(['Single electron wavefunction of neutral helium ',num2str(count), 'iterations' ]);
    figure(2); plot(r,-1./r,r,U_scf,'--'); xlabel('distance[nm]');ylabel('Energy (eV)'); 
    title('Single electron binding energy of neutral helium'); axis([0 .5 -100 20]);
    
    for i = 1:length(r); %calculate U_scf in the iteration
        U_scf(i) = (A/r(i))*sum(sigma(1:i)*dr) + ...
            A*sum(sigma(i:end)./r(i:end));       
    end
    U2 = U_scf; 
    diffU = max(abs(U2 - U1));
    count = count +1; 
end


%%  the self-consistent electron density for a 9nm QW

clear all; close all; 

%set constants
N = 500; h = 1; m = 1; a = .5*N; V0 = -3; B = 30; 
mu = 0; dz = .01; e = 1; e0 = 1; p0 = 0; %positive background charge
z = -6+12/N:12/N:6;
epsilon = 4;
alpha = 0.1; 

%second derivative finite differences operator in 1D
D2z=1/dz^2*(spdiags(-2*ones(N,1),0,sparse(N,N))+spdiags(ones(N,1),1,sparse(N,N))+spdiags(ones(N,1),-1,sparse(N,N)));
Vz = [zeros(1,.25*N),V0*ones(1,a),zeros(1,.25*N)]; 
V = sparse(diag(Vz));
%forward difference and backward difference d/dz operators
D1f = (1/dz)*(spdiags(-ones(N,1),0,sparse(N,N))+spdiags(ones(N,1),-1,sparse(N,N)));
D1b = (1/dz)*(spdiags(ones(N,1),0,sparse(N,N))+spdiags(-ones(N,1),1,sparse(N,N)));

Uz = zeros(1,N); %initial guess
diffp = 1; count = 1;  
while diffp > 1e-4; %1meV tolerance
U1 = Uz; 
H = -(h^2/2/m)*D2z + V + sparse(diag(U1)); %Hamiltonian FD operator
[F,D] = eig(full(H)); %calculate eigenvalues and eigenvectors

En = diag(D); %energy eigenvalues of solutions
fn = 1./(1+exp(B*(En-mu))); %Fermi Dirac function
pr = 2*fn'*abs(F').^2; %electron density

M = -D1f*epsilon*D1b - (.01)*(sparse(diag(pr - p0))); %e^2/e0
b = zeros(N,1); b(.25*N-1:.75*N+1,1) = .25;  
Uz = (M\b)'; %calculate U(z) by solving poisson's equation
Uz = U1 + alpha*(Uz-U1); %numerical damping
Uz = [fliplr(Uz(N/2+1:end)),Uz(N/2+1:end)];
U2 = Uz;

diffp = max(abs(U2-U1)); count = count + 1;

end

%PLOT RESULTS%
figure(1); plot(z,pr*166,z,-V0+Vz+Uz); xlabel('position[nm]'); ylabel('density'); 
title(['Self consistent electron density for a 9nm QW - ', num2str(count), ' iterations']);

%% delta doping case

clear all; close all; 

%set constants
N = 500; h = 1; m = 1; a = .5*N; V0 = -3; B = 30; 
mu = 0; dz = .01; e = 1; e0 = 1; p0 = 1; %positive background charge
z = -10+12/N:12/N:2;
epsilon = 4;
alpha = 0.1; 

%second derivative finite differences operator in 1D
D2z=1/dz^2*(spdiags(-2*ones(N,1),0,sparse(N,N))+spdiags(ones(N,1),1,sparse(N,N))+spdiags(ones(N,1),-1,sparse(N,N)));
Vz = [0*ones(1,.75*N),V0*ones(1,.25*N)]; 
V = sparse(diag(Vz));
%forward difference and backward difference d/dz operators
D1f = (1/dz)*(spdiags(-ones(N,1),0,sparse(N,N))+spdiags(ones(N,1),-1,sparse(N,N)));
D1b = (1/dz)*(spdiags(ones(N,1),0,sparse(N,N))+spdiags(-ones(N,1),1,sparse(N,N)));

Uz = zeros(1,N); %initial guess
diffp = 1; count = 1;  
while diffp > 1e-4; %1meV tolerance
U1 = Uz; 
H = -(h^2/2/m)*D2z + V + sparse(diag(U1)); %Hamiltonian FD operator
[F,D] = eig(full(H)); %calculate eigenvalues and eigenvectors

En = diag(D); %energy eigenvalues of solutions
fn = 1./(1+exp(B*(En-mu))); %Fermi Dirac function
pr = 2*fn'*abs(F').^2; %electron density

M = -D1f*epsilon*D1b - (.01)*(sparse(diag(pr - p0))); %e^2/e0
b = zeros(N,1); b(1:.75*N,1) = .25; 
Uz = (M\b)'; %calculate U(z) by solving poisson's equation
Uz = U1 + alpha*(Uz-U1); %numerical damping
U2 = Uz; 

diffp = max(abs(U2-U1)); count = count + 1;

end

%PLOT RESULTS%
UU = 3.5-(-V0+Vz+Uz); UU = 1+[-UU(1:.75*N),5.5-UU(.75*N+1:end)];
figure(1); plot(z-2,pr*166,z,UU); xlabel('position[nm]'); ylabel('density');
axis([-10 2 0 4]); title(['Self consistent electron density for a 9nm QW with delta doping - ', num2str(count), ' iterations']);
