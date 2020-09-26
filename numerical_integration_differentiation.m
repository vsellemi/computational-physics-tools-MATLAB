%%  02/07, Victor Sellemi
%% catastrophic cancellation in function evaluations
%(1)(i): f(x) = 1 - cosx
x = [-4e-08:0.01e-08:4e-08]; %set the domain of the plot
y = 1 - cos(x); %evaluate 1-cosx at the domain
z = 2*(sin(x/2)).^2; %evaluate 2sin^2(x/2) at the domain
plot(x,y); hold on; %plot 1-cosx and keep the plot open
plot(x,z); %plot 2sin^2(x/2) on the same plot
legend('1-cos(x)', '2sin(x/2)^2','Location', 'NorthEast'); %add a legend
title('Catastrophic cancellation of 1-cosx'); xlabel('x'); %add titles
hold off;clear all;

%%  f(x) = 1 - expx 
x = [-1e-15:1e-17:1e-15]; %set the domain of the function
y = exp(x) - 1; %expx - 1
w = tanh(x/2).*(exp(x) + 1); %rewritten to avoid catastrophic cancellation
z = expm1(x); %z = exmp1(x), which is equivalent to above function w(x)
plot(x,y); hold on %plot of all three functions evaluated on domain
plot(x,w); hold on
plot(x,z);
axis([-1e-15 1e-15 -1.5e-15 1.5e-15]); %define axis bounds and labels
legend('exp(x)-1', 'tanh(x/2)(exp(x)+1)', 'expm1(x)', 'Location', 'NorthEast');
title('Catastrophic cancellation of expx-1'); xlabel('x') %add titles 
hold off; clear all;

%%  quadratic formula
a = 1; c = 0.5; %fix a and c values
b = round(10.^(6:0.25:8)); %set range of b values
x1 = (-b + sqrt(b.^2 - 4*a*c))/(2*a); %x1 solution based on quadratic form
x2 = (2*c)./(-b-sqrt(b.^2-4*a*c)); %avoid catastrophic cancellation
x = a*c*(b.^(-2)); %set x = ac/b^2
y = abs((x2 - x1)./(x2)); %relative error

%build figure with two plots
%Figure 1 (3-D): X~a/b, Y~b/c, Z = relative error
%Figure 2 (2-D): x = ac/b^2, y = relative error
X = logspace(-8,-6,10); Y = logspace(-8,-6,10); [X,Y] = meshgrid(X,Y);
X2 = (-2*Y)./(1+sqrt(1-(4*X.*Y)));
X1 = (-1./X)*0.5*(1-sqrt(1-(4*X.*Y)));
Z = abs((X1-X2)./X2);
figure %generate figure for the two plots
subplot(1,2,1); %add first plot
S = surf(X,Y,Z); xlabel('a/b'); ylabel('c/b'); zlabel('relative error');
title('Relative error of solution to quadratic formula');
subplot(1,2,2); %add second plot
plot(log(x),y); xlabel('ac/b^2'); ylabel('relative error');
title('Catastrophic cancellation in quadratic formula');
clear all; 

%% instability of iterative schemes to calculate powers of phi
close all;
phi = ((5^0.5) - 1)/2; %initial value of phi
i = 1:80; %vector of integers
phi_i(i) = phi.^i; %generate vector of powers of phi by multiplication
phi_j = []; %initialize powers of phi vector for iterative method
phi_j(1) = ((5^0.5) - 1)/2; %initialize first term phi^1
phi_j(2) = phi_j(1)^2; %initialize second term phi^2
for j = 2:79; %iterate for powers of phi from 3-80
    phi_j(j+1) = phi_j(j-1) - phi_j(j); %iteration scheme using previous powers of phi
end
plot(i,log(phi_i)); hold on %plot results using multiplication
plot(i,log(phi_j)); %plot results using iterative scheme
legend('multiplication', 'linear recursion'); %add legend
title('Instability of iterative methods for powers of phi'); %add title
xlabel('n'); ylabel('log phi^n'); % add axis labels
hold off;clear all;

%%  scaling of absolute error
%(i)forward and centered difference formulae for df/dx|x=0 with f(x)=exp(x)
gf_h = @(x) abs(((exp(x)-1)./x) - 1); %forward difference error function
gc_h = @(y) abs(((exp(y)-exp(-y))./(2.*y)) - 1); %centered difference error function
h = 1 * 10.^-(fliplr(0:0.5:10)); %vector of h values
error_fd = gf_h(h); %forward difference error evaluated on h
error_cd = gc_h(h); %centered difference error evaluated on h
plot(log(h),log(error_fd)); hold on %plot error on log-log plot
plot(log(h),log(error_cd));
%plot labels
legend('forward', 'centered');
xlabel('log of h'); ylabel('log of abs error d/dx|x=0 expx')
title('Absolute error of forward and centered difference formulas for d/dx expx at x=0');
hold off; clear all;

%% centered difference for second derivative at zero with f(x)=exp(x)
f = @(x) exp(1).^x; % f(x) = exp(x)
gc_h = @(h) abs(((f(h)-2*f(0)+f(-h))./(h.^2)) - 1); %abs error formula 
h = 10.^-(fliplr(1:8)); %h vector
y = gc_h(h); %error evaluated at h
plot(log(h),log(y)); %plot absolute error on log-log plot
%plot labels
ylabel('log of absolute error in d^2/dx^2 exp(x)');
xlabel('log of h');
title('Absolute error for centered difference formula for d^2/dx^2 exp(x)');
clear all;
%% trapezoidal & simpson's rules for integrating exp(x) from 0 to 1
f = @(x) exp(x); %f(x) = exp(x)
a = 0; b = 1; %set bounds of integration
%initialize arrays for loop
int_trapezoidal = []; int_simpsons = []; int_actual = []; abs_error =[];
X = [];

for j = 1:3; %loop over different values of N
    N = 10^(j); %set value of N for desired accuracy
    X(j) = log(N);
    h = 1/N; i = 1:(N-1); %set value of h based on N 
    %calculate integral estimates using trapezoid rule and simpsons rule
    int_trapezoidal(j) = (h/2)*(f(a) + 2*sum(f(a + i*h)) + f(b));
    int_simpsons(j) = (h/3)*(f(a) + 4*sum(f(a + h*i(1:2:end))) + 2*sum(f(a+h*i(2:2:end))) + f(b));
    int_actual = exp(b) - exp(a); %analytic result for the integral using FTC
    abs_error(j,1) = abs(int_trapezoidal(j) - int_actual); %absolute error array
    abs_error(j,2) = abs(int_simpsons(j) - int_actual);
end
Y = log(abs_error);
plot(X,Y(:,1)); hold on
plot(X,Y(:,2)); legend('trapezoid rule','simpsons rule');
xlabel('N'); ylabel('absolute error'); 
title('integrating exp(x) using trapezoid and simpsons rules');
hold off; clear all;

%% monte carlo method for integrating expx from 0 to 1
f = @(x) exp(x); %f(x) = exp(x)
abs_error = []; X = []; %initialize absolute error vector
for i = 1:13; %loop over different values of N
N = round(10^(2 + i*.25)); 
X(i) = N;
x = rand(1,N); %randomize the x component 
y = exp(1)*rand(1,N); %randomize the y component 
n_inside = sum(y < f(x)); %count the number of points that fall under the curve
estimate = (n_inside./N)*exp(1); %estimate the area underneath the curve
actual = exp(1) - 1; %evaluate the actual integral using FTC
abs_error(i) = abs(estimate - actual); %calculate the absolute error
Z(i) = 1/(N^0.5); %vector for f = 1/(N^0.5) to compare
end
Y = abs_error; 
plot(log(X),log(Y)); hold on; %plot the absolute error for various N
plot(log(X),log(Z),'--');
legend('monte carlo', 'N^-0.5');
title('Absolute error of monte carlo integration of exp(x) from 0 to 1');
ylabel('Absolute error'); xlabel('log(N)');
hold off; clear all;

%% volume of the 2,4,and 6-d unit hypersphere using Monte-Carlo
dim = [2,4,6]; %vector of dimensions
abs_error = []; %initialize array of absolute error
X = []; %initialize vector of X values for later plotting
for i = 1:8; %loop over N values
   N = round(10^(.5*i)); %set N to be an integer value ~ powers of 10
   X(i) = log(N);
   for k = 1:length(dim); %loop over the three dimension values
       x = -1 + 2*rand(N,dim(k)); %generate random dim dimension points
       r = sum((x.*x)'); %calculate the magnitudes of the random points
       n_inside = sum(r < 1); %count the number of points with r < 1
       V_estimate(k) = (n_inside/N)*(2^dim(k)); %generate volume estimate
   end
   V_actual = [pi, 0.5*pi^2, (1/6)*pi^3]; %actual volume values
   abs_error(i,1:3) = abs(V_estimate - V_actual); %calc. absolute error
end
Y = log(abs_error); %take natural log of abs_error values for plot
plot(X/2,Y(:,1)); hold on; %plot volume absolute error vs. samples/dim.
plot(X/4,Y(:,2)); hold on;
plot(X/6,Y(:,3)); hold on;
%plot labels
legend('2D', '4D', '6D');
xlabel('Monte carlo samples per dimension'); 
ylabel('Unit hypersphere volume absolute error');
title('Monte carlo method for calculating volume of 2,4,and 6-d unit hypersphere');
hold off; clear all;