%% 02/21, Victor Sellemi

%% 2d Laplace's equation with edge boundary conditions

%We implement the relaxation method to estimate the solution of del^2f = 0.
%We first construct a Nx x Ny matrix (M) of values that f(x,y) could take
%evaluated at the point (i,j) of the grid. We then apply the boundary 
%conditions and implement an initial guess of the solution by averaging 
%over the initial conditions. We run iterations of the relaxation method
%that are given by M = AM + MA' where A is a Nx x Ny matrix of zeros with
%0.25 at every entry on the first superdiagonal and first subdiagonal 

Nx = 20; Ny = 20; %set the size of an Nx by Ny grid
M = zeros(Nx,Ny); %initialize the Nx by Ny array
M(:,1) = 0; M(1,:) = 1; M(:,end) = 1; M(end,:) = 0; %initial conditions
M(:,2:end-1) = mean(M(:)); %initial guess given by the average over B.C's

%iterations of relaxation method are given by M = AM + MA' 
%define A below 
A = 0.25*(diag(ones(1,Ny-1),-1) + diag(ones(1,Ny-1),1));

diff = 1; %first value of variation threshold
while max(diff(:)) > 1e-6; %iterate until 10^-6 variation threshold is met
    M = A*M + M*A'; %relaxation method iterations
    %reapply boundary conditions each iteration
    M(:,1) = 0; M(:,end) = 1; M(1,:) = 1;  M(end,:) = 0; 
    %calculate absolute deviation between iterations
    diff1 = abs(M-(A*M + M*A')); diff = diff1(2:end-1,2:end-1);     
end

%plot solution
surf(M); xlabel('x'); ylabel('y'); zlabel('f(x,y)');
title('Solution to 2D Laplace equation with edge boundary conditions');

%%  2d Laplace's equation with interior boundary condition
clf;
%we implement the same relaxation method from above except with different
%initial conditions. In this case, we hold the value of f(x,y) = 0 at all
%of the boundaries and f(x,y) = 1 in the center of the array

Nx = 20; Ny = 20; %set the size of an Nx by Ny grid
M = zeros(Nx,Ny); %initialize the Nx by Ny array
%initial conditions
M(:,1) = 0; M(1,:) = 0; M(:,end) = 0; M(end,:) = 0; 
M(round(Nx/2),round(Ny/2))=1;
M(:,2:end-1) = mean(M(:)); %initial guess given by the average over B.C's

%iterations of relaxation method are given by M = AM + MA' 
%define A below 
A = 0.25*(diag(ones(1,Ny-1),-1) + diag(ones(1,Ny-1),1));

diff = 1; %first value of variation threshold
while max(diff(:)) > 1e-6; %iterate until 10^-6 variation threshold is met
    M = A*M + M*A'; %relaxation method iterations
    %reapply initial conditions each iteration
    M(:,1) = 0; M(:,end) = 0; M(1,:) = 0;  M(end,:) = 0; M(Nx/2,Ny/2)=1;
    %calculate absolute deviation between iterations
    diff1 = abs(M-(A*M + M*A')); diff = diff1(2:Nx/2-1,2:Ny/2-1); 
end

%plot solution
surf(M); xlabel('x'); ylabel('y'); zlabel('f(x,y)');
title('Solution to 2D Laplace equation with interior boundary conditions');

%%  Processing time of relaxation method as a function of N=Nx=Ny
clf;

%here, we attempt to quantify the relationship between the dimension of the
%grid N = Nx = Ny and the time it takes to produce a solution using the
%relaxation method with a 10^-6 maximum absolute variation threshold. Here
%we will test for the case where the boundary conditions are at the edges

time = []; N = [];

for i = 1:70;
tic
N(i) = i; Nx = i; Ny = i;    
M = zeros(Nx,Ny); %initialize the Nx by Ny array
M(:,1) = 0; M(1,:) = 1; M(:,end) = 1; M(end,:) = 0; %initial conditions
M(:,2:end-1) = mean(M(:)); %initial guess given by the average over B.C's

%iterations of relaxation method are given by M = AM + MA' 
%define A below 
A = 0.25*(diag(ones(1,Ny-1),-1) + diag(ones(1,Ny-1),1));

diff = 1; %first value of variation threshold
while max(diff(:)) > 1e-6; %iterate until 10^-6 variation threshold is met
    M = A*M + M*A'; %relaxation method iterations
    %reapply boundary conditions each iteration
    M(:,1) = 0; M(:,end) = 1; M(1,:) = 1;  M(end,:) = 0; 
    %calculate absolute deviation between iterations
    diff1 = abs(M-(A*M + M*A')); diff = diff1(2:end-1,2:end-1);     
end

time(i) = toc;

end

plot(N,time,'-*'); xlabel('N=Nx=Ny'); ylabel('processing time');
title('Processing time of relaxation method');

%the processing time appears to increase quadratically over time for N
%values between 1 and 70,

%%  Laplacian matrix operator to solve 2D laplace's equation
clf;
time = []; stime = []; N = [];
for i = 1:70;

N(i) = i; Nx = i; Ny = i; dx = 1; dy = 1; 

tic
%full matrices
L1x = (1/dx^2)*(-2*eye(Nx)+diag(ones(Nx-1,1),1)+diag(ones(Nx-1,1),-1));
L1y = (1/dy^2)*(-2*eye(Ny)+diag(ones(Ny-1,1),1)+diag(ones(Ny-1,1),-1));
L2 = kron(eye(Ny),L1x)+kron(L1y,eye(Nx));

M = zeros(Nx,Ny);
M(round(Nx/2),round(Ny/2)) = 1;
Y = L2\reshape(-M,Nx*Ny,1);

time(i) = toc;

tic
%sparse matrices
sL1x=1/dx^2*(spdiags(-2*ones(Nx,1),0,sparse(Nx,Nx))+spdiags(ones(Nx,1),1,sparse(Nx,Nx))+spdiags(ones(Nx,1),-1,sparse(Nx,Nx)));
sL1y=1/dy^2*(spdiags(-2*ones(Ny,1),0,sparse(Ny,Ny))+spdiags(ones(Ny,1),1,sparse(Ny,Ny))+spdiags(ones(Ny,1),-1,sparse(Ny,Ny))); 
sL2=kron(speye(Ny),sL1x)+kron(sL1y,speye(Nx));

M = zeros(Nx,Ny);
M(round(Nx/2),round(Ny/2)) = 1;
Y = sL2\reshape(-M,Nx*Ny,1);

stime(i) = toc;

end

plot(N,time); hold on; plot(N,stime); xlabel('N = Nx = Ny'); 
ylabel('processing time'); legend('full matrices','sparce matrices');
title('scaling behavior of laplaces equation using sparse and full matrices');

%%  2d Laplace's equation on circular disk (zero at perimeter)
clf;
%we attempt to solve 2d laplace's equation on a circular disk with zero
%on the perimeter and unity in the center. We use matrix operator form of
%d/dr and d^2/dr^2 and assume approximate rotational symmetry so we can
%ignore the differential term with respect to theta. 

Nr = 80; R = 1; dr = R/Nr;

%analytic solution
r = 0+1/Nr:1/Nr:1;
Z = @(r) 1 + log(dr./r)/log(R/dr);
Z_analytic = Z(r);

%calculated solution using sparse matrix operator form of d/dr and d2/dr2
%assume d/dth is zero due to rotational symmetry
L2r = sparse((1/dr^2)*(-2*eye(Nr)+diag(ones(Nr-1,1),1)+diag(ones(Nr-1,1),-1)));
L1r = sparse((1/(2*dr))*(diag(-ones(Nr-1,1),-1)+diag(ones(Nr-1,1),1)));
L1r(1,1) = -2*1/(2*dr); L1r(1,2) = 2*1/(2*dr); L2r(1,1) = -1*1/dr^2;
L2 = L2r + diag(1./linspace(2*dr,R,Nr))*L1r;

M = zeros(Nr,1);
M(1,1) = -3/4*(1/dr^2);
Z_calc = L2\-M;

%plot the results
subplot(2,1,1); plot(r,-Z_calc,'o'); hold on; plot(r,Z_analytic); 
xlabel('position r'); ylabel('f(r)'); legend('calc','analytic');

subplot(2,1,2);
Nths=30; 
Z=repmat(Z_calc,1,Nths); 
ths=linspace(0,2*pi,Nths);
[THS,RS]=meshgrid(ths,r);
[X,Y]=pol2cart(THS,RS); 
surf(X,Y,-Z);


%%  2d Laplace's equation on a circular disk
clf;
Nr = 80; R = 1; dr = R/Nr; Nth = 80; dth = 2*pi/Nth;

%calculated solution using sparse matrix operator form of d/dr and d2/dr2
L2r = sparse((1/dr^2)*(-2*eye(Nr)+diag(ones(Nr-1,1),1)+diag(ones(Nr-1,1),-1)));
L1r = sparse((1/(2*dr))*(diag(-ones(Nr-1,1),-1)+diag(ones(Nr-1,1),1)));
L2th = sparse((1/dth^2)*(-2*eye(Nth)+diag(ones(Nth-1,1),1)+diag(ones(Nth-1,1),-1)));
L1r(1,1) = -2*1/(2*dr); L1r(1,2) = 2*1/(2*dr); L2r(1,1) = -1*1/dr^2; 
L2th(1,end) = 1/dth^2; L2th(end,1) = 1/dth^2;
r = spdiags(ones(Nr,1),0,sparse(Nr,Nr));
r2 = spdiags(ones(Nr,1),0,sparse(Nr,Nr));

%diag(1./linspace(2*dr,R,Nr))
L2 = kron(speye(Nth),L2r) + kron(speye(Nth),r*L1r) + kron(speye(Nth),r2)*kron(L2th,speye(Nr));

M = zeros(Nr*Nth,1);
M(1,1) = -3/4*(1/dr^2)*dr;
M(round(Nr/2),1) = 1;
B = L2\M;
Z = reshape(B,Nr,Nth);
r = 0+1/Nr:1/Nr:1;
ths=linspace(0,2*pi,Nth);
[THS,RS]=meshgrid(ths,r);
[X,Y]=pol2cart(THS,RS);
surf(X,Y,Z);



%%  2d Laplace's equation on a circular disk
clf;
Nr = 80; R = 1; dr = R/Nr; Nth = 80; dth = 2*pi/Nth;

%calculated solution using sparse matrix operator form of d/dr and d2/dr2
L2r = sparse((1/dr^2)*(-2*eye(Nr)+diag(ones(Nr-1,1),1)+diag(ones(Nr-1,1),-1)));
L1r = sparse((1/(2*dr))*(diag(-ones(Nr-1,1),-1)+diag(ones(Nr-1,1),1)));
L2th = sparse((1/dth^2)*(-2*eye(Nth)+diag(ones(Nth-1,1),1)+diag(ones(Nth-1,1),-1)));
L1r(1,1) = -2*1/(2*dr); L1r(1,2) = 2*1/(2*dr); L2r(1,1) = -1; 
L2th(1,end)=1; L2th(end,1)=1;
L2 = kron(speye(Nth),L2r) + kron(speye(Nth),diag(1./linspace(2*dr,R,Nr))*L1r) + kron(speye(Nth),diag(1./(linspace(2*dr,R,Nr).^2))) + kron(L2th,speye(Nr));

M = zeros(Nr,Nth);
theta = 0+1/Nth:1/Nth:1;
M(end,:) = sin(2*theta);
Z = L2\reshape(-M,Nr*Nth,1);
r = 0+1/Nr:1/Nr:1;
ths=linspace(0,2*pi,Nth);
[THS,RS]=meshgrid(ths,r);
[X,Y]=pol2cart(THS,RS);
surf(X,Y,reshape(Z,Nr,Nth))









