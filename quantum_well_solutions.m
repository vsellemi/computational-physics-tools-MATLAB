%% 04/13, Victor Sellemi

%% lowest eigenvalues of a 1D potential for an electron
%harmonic oscillator, hw = 1eV 

N = 200; dx = 1; h = 1; m = 1; %Set constants
w = 1/h; %hw = 1eV condition

%second derivative finite differences operator in 1D
D2x=1/dx^2*(spdiags(-2*ones(N,1),0,sparse(N,N))+spdiags(ones(N,1),1,sparse(N,N))+spdiags(ones(N,1),-1,sparse(N,N)));
Vx = 0.5*m*(w^2)*(1/N-1/2:1/N:1/2).^2; %V(x) = .5m(w^2)(x^2)
V = sparse(diag(Vx)); %potential operator
H = -(h^2/2/m)*D2x + V; %Hamiltonian FD operator
[F,D] = eigs(H,10,'SM'); %calculate eigenvalues and eigenvectors
En = fliplr(diag(D)'); %eigenvalues in order from n=1,...,10

%plot results
figure(1); plot(N*En,'o'); xlabel('quantum number'); ylabel('Energy'); 
title('10 lowest eigenvalues of 1D harmonic oscillator potential'); 

X = 1:1:N;
Y1 = F(:,end).^2; Y2 = F(:,end-1).^2;
figure(2); plot(X,Y1,X,Y2,X,Vx./sqrt(N)); 
legend('mode 1', 'mode 2', 'Vx - scaled'); ylabel('|psi|^2'); xlabel('x'); 
title('lowest two eigenstates superimposed on V(x)')

%% lowest eigenvalues of a 1D potential for an electron
%finite quantum well, 1eV deep and 1nm wide
close all; clear all;

N = 200; dx = .1; h = 1; m = 1; %Set constants
a = .25*N; V0 = -1; 

%second derivative finite differences operator in 1D
D2x=1/dx^2*(spdiags(-2*ones(N,1),0,sparse(N,N))+spdiags(ones(N,1),1,sparse(N,N))+spdiags(ones(N,1),-1,sparse(N,N)));
Vx = [zeros(1,.375*N),V0*ones(1,a),zeros(1,.375*N)]; %V(x) = [0 V0 0]
V = sparse(diag(Vx)); %potential operator
H = -(h^2/2/m)*D2x + V; %Hamiltonian FD operator
[F,D] = eigs(H,10,'SM'); %calculate eigenvalues and eigenvectors
En = fliplr(diag(D)'); %eigenvalues in order from n=1,...,10

%plot results
figure(1); plot(1:10,En,'o',1:10,zeros(1,10),'-r'); xlabel('quantum number'); ylabel('Energy'); 
title('10 lowest eigenvalues of 1D finite quantum well potential'); 

X = 1:1:N;
Y1 = F(:,end).^2; Y2 = F(:,end-1).^2;
figure(2); plot(X,Y1,X,Y2,X,Vx./sqrt(N));
legend('mode 1', 'mode 2', 'Vx - scaled'); ylabel('|psi|^2'); xlabel('x'); 
title('lowest two eigenstates superimposed on V(x)')

%%  lowest eigenvalues of a 1D potential for an electron
%double quantum well, 1nm separation, 1eV deep
close all; clear all;

N = 200; dx = .1; h = 1; m = 1; %Set constants
a = .075*N; V0 = -1; 

%second derivative finite differences operator in 1D
D2x=1/dx^2*(spdiags(-2*ones(N,1),0,sparse(N,N))+spdiags(ones(N,1),1,sparse(N,N))+spdiags(ones(N,1),-1,sparse(N,N)));
Vx = [zeros(1,.375*N),V0*ones(1,a),zeros(1,.1*N),V0*ones(1,a),zeros(1,.375*N)]; %V(x) = [0 V0 0 V0 0]
V = sparse(diag(Vx)); %potential operator
H = -(h^2/2/m)*D2x + V; %Hamiltonian FD operator
[F,D] = eigs(H,10,'SM'); %calculate eigenvalues and eigenvectors
En = fliplr(diag(D)'); %eigenvalues in order from n=1,...,10

%plot results
figure(1); plot(1:10,En,'o',1:10,zeros(1,10),'-r'); xlabel('quantum number'); ylabel('Energy'); 
title('10 lowest eigenvalues of 1D double quantum well potential'); 

X = 1:1:N;
Y1 = F(:,end).^2; Y2 = F(:,end-1).^2;
figure(2); plot(X,Y1,X,Y2,X,Vx./sqrt(N));
legend('mode 1', 'mode 2', 'Vx - scaled'); ylabel('|psi|^2'); xlabel('x'); 
title('lowest two eigenstates superimposed on V(x)')

%%  lowest eigenvalues of a 1D potential for an electron
%exponential well, ?V0 exp(?|x|/x0): V0 =2 eV, x0 =1 nm.
close all; clear all;

N = 200; dx = .1; h = 1; m = 1; %Set constants
V0 = 2; x0 = 1; 

%second derivative finite differences operator in 1D
D2x=1/dx^2*(spdiags(-2*ones(N,1),0,sparse(N,N))+spdiags(ones(N,1),1,sparse(N,N))+spdiags(ones(N,1),-1,sparse(N,N)));
Vx = -V0*exp(-abs(.1*(1-N/2:1:N/2))./x0); %V(x) = ?V0 exp(?|x|/x0)
V = sparse(diag(Vx)); %potential operator
H = -(h^2/2/m)*D2x + V; %Hamiltonian FD operator
[F,D] = eigs(H,10,'SM'); %calculate eigenvalues and eigenvectors
En = fliplr(diag(D)'); %eigenvalues in order from n=1,...,10

%plot results
figure(1); plot(1:10,En,'o',1:10,zeros(1,10),'-r'); xlabel('quantum number'); ylabel('Energy'); 
title('10 lowest eigenvalues of 1D exponential well potential'); 

X = 1:1:N;
Y1 = F(:,end).^2; Y2 = F(:,end-1).^2;
figure(2); plot(X,Y1,X,Y2,X,Vx./sqrt(N));
legend('mode 1', 'mode 2', 'Vx - scaled'); ylabel('|psi|^2'); xlabel('x'); 
title('lowest two eigenstates superimposed on V(x)')

%% Variation in the ground state energy of a finite quantum well
close all; clear all;

h = 1; m = 1; dx = 20/500; N = 500; x = 0:20/499:20; W = 20:10:160; %set initial parameters

%second derivative FD operator
D2x=1/dx^2*(spdiags(-2*ones(N,1),0,sparse(N,N))+spdiags(ones(N,1),1,sparse(N,N))+spdiags(ones(N,1),-1,sparse(N,N)));

En1 = []; En2 = []; %initialize solution arrays

for k = 1:length(W); %loop over different values of well width     
a = W(k); V0 = -2; b = (N - a)/2;
Vx = [zeros(1,b),V0*ones(1,a),zeros(1,b)]; %V(x) = [0 V0 0]
V = sparse(diag(Vx));
H = -(h^2/2/m)*D2x + V; 
[F,D] = eig(full(H)); 
En(k) = -D(1,1)-2; 
end

%plot results
plot(W./max(W),En,'-o'); xlabel('well width [nm]'); ylabel('ground state energy');
title('Ground state energy of finite quantum well as a function of width');

%%  Variation in energies of the two lowest bound states of a finite double quantum wel
close all; clear all;

h = 1; m=1; dx = 20/500; N = 500; x = 0:20/499:20; W = 20:10:160; %set initial parameters

%second derivative FD operator
D2x=1/dx^2*(spdiags(-2*ones(N,1),0,sparse(N,N))+spdiags(ones(N,1),1,sparse(N,N))+spdiags(ones(N,1),-1,sparse(N,N)));

En1 = []; En2 = []; %initialize solution arrays

for k = 1:length(W); %loop over different values of barrier width
a = W(k); b = (N-a)/2-25; 
Vx = [zeros(1,b), -ones(1,25), zeros(1,a), -ones(1,25), zeros(1,b)]; 
V = sparse(diag(Vx)); %potential operator
H = -(h^2/2/m)*D2x + V; %full operator
[F,D] = eig(full(H)); %eigenvalue solutions of initial operator
En1(k) = D(1,1); %store lowest two eigenvalues
En2(k) = D(2,2);
end

%plot results
plot(W./100,En1); hold on; plot(W./100,En2); legend('symmetric', 'antisymmetric');
xlabel('barrier width [nm]'); ylabel('lowest energy eigenvalues');
title('Variations in energy of lowest two bound states double quantum well');

%%  1 nm-wide infinite quantum well in the range 0.5 < E < 2 eV
%Numerov integration and shooting method
close all; clear all;

%use the bisection algorithm to find E values that minimize deviation of
%numerov wavefunction from boundary conditions, code is commented out since
%value of 0.5 seems to provide the best results

%F = @(e) numerov(e); 
%E1 = bisection(F,.5,2,1e-2); 

%plot results 
EE = .5:.01:.55;

for i = 1:length(EE); %loop over different values of energy
h = 1; m = 1; dx = .001; h = dx; N = 200; E = EE(i); %set parameters

W = zeros(1,N); psi = zeros(1,N);
V = 0; f = (2*m/h)*(V-E); 
psi(1) = 0; psi(2) = 0.1;
W(1) = (1-(h^2/12)*f)*psi(1); 
 
for n = 2:N-1; 
W(n) = (1-(h^2/12)*f)*psi(n); 
W(n+1) = h^2*f*psi(n) + 2*W(n) - W(n-1);
psi(n+1) = W(n+1)/(1-(h^2/12)*f);
end

%plot results
x = 1/N:1/N:1;
plot(x,-psi); hold on;
xlabel('position [nm]'); ylabel('psi');
title('waveforms of infinite quantum well at different energy values'); 
end


%%  the exponential potential V (x) = ?V0 exp(?|x|/x0)
close all; clear all;

h = 1; m = 1; dx = .006*1e-3; h = dx; N = 650; E = 2; %set parameters
x0 = 1; V0 = 1; 

W = zeros(1,N); psi = zeros(1,N);
V = -V0*exp(-abs(.1*(1-N/2:1:N/2))./x0); 
f(1) = (2*m/h)*(V(1)-E); 
psi(1) = 0; psi(2) = 0.1;
W(1) = (1-(h^2/12)*f(1))*psi(1); 
 
for n = 2:N-1; 
f(n) = (2*m/h)*(V(n)-E); 
W(n) = (1-(h^2/12)*f(n))*psi(n); 
W(n+1) = h^2*f(n)*psi(n) + 2*W(n) - W(n-1);
psi(n+1) = W(n+1)/(1-(h^2/12)*f(n));
end

%plot results
x = 1/N:1/N:1;
plot(x,abs(psi).^2/N); 
xlabel('distance'); ylabel('|psi|^2');
title('waveforms of exponential at different energy values');

%%  the energy and  wavefunctions of 1st 10 bound states of infinite QW
close all; clear all; 

N = 16; a = 1; h = 1; m0 = 1;
x = -a:2*a/200:a; 

f = []; 
for i = 1:N;
    n = i-1;
    f(i,:) = (x.^(n+2)-(a^2*x.^n))/(a^(n+2.5));
end
H = zeros(N,N); S = zeros(N,N);
for i = 1:N;
    n = i-1;
    for j = 1:N;
        m = j-1;
        if mod(n+m,2) == 0;
            H(i,j) = (-h^2/2/m0/a^2).*8.*(1-m-n-2.*m.*n)./(n+m+3)./(n+m+1)./(n+m-1);
            S(i,j) = 2./(n+m+5)-4./(n+m+3)+2./(n+m+1);
        elseif mod(n+m,2) == 1;
            H(i,j) = 0; S(i,j) = 0; 
        end
    end
end

[F,D] = eig(H,S); 
Y = (F'*f)'; Y = abs(Y(:,1:3));

%plot results
figure(1); plot(x,Y.^2); xlabel('position'); ylabel('|psi|^2');
title('bound states of an infinite QW');
figure(2); En = diag(D)./N; plot(En(1:end),'o'); 
title('bound energy of infinite QW');
xlabel('quantum number'); ylabel('energy'); 








