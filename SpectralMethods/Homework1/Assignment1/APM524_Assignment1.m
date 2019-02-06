close all; clear variables

%% Problem 4
N = 20;
x = chebfun('x',[0 2*pi]);
f = sin(x./2).^3;
figure
plot(f)

A = exp(1i*x*(-N:N));
lambda = (1/2*pi)*A'*f;

% figure
% plot(real(lambda),imag(lambda),'*')

figure
loglog(-N:N,lambda,'*')

%%
clear variables;
N = 20:10:100;
x = chebfun('x',[0 2*pi]);
f = sin(x./2).^3;
for j = 1:length(N)
    A = exp(1i*x*(-N(j):N(j)));
    lambda = (1/2*pi)*A'*f;
    fn = A*lambda;
    err(j) = norm(f-fn,2);
end

figure
loglog(err,'*')

%% Problem 5

N = 500;

x = chebfun('x',[0 2*pi]);
f = 3./(5-4*cos(x));
figure
plot(f)

A = exp(1i*x*(-N:N));
lambda = (1/2*pi)*A'*f;

a = acos(5/4)
fa = 3./(5-4*cos(x-a*1i));
M = max(fa);
errbound = 4*pi*M/(a*(exp(-a*N)))

figure
loglog(-N:N,lambda,'*')
hold on
loglog(-N:N,errbound*ones(length(-N:N)),'r')

%% 
N = 20:10:100;
for i = 1:length(N)
x = chebfun('x',[0 2*pi]);
f = 1./(2+sin(x).^2);
A = exp(1i*x*(-N(i):N(i)));
lambda = (1/(2*pi)*A'*f);
fn = A*lambda;

err(i) = norm(f - fn,2);
end
%%
figure
loglog(N,err,'*')

% Z = asin(sqrt(2)*1i);

% figure
% semilogy(-N:N,lambda,'*')

%%
x = chebfun('x',[0 2*pi]);
f = 1./(2+sin(x).^2);
A = exp(1i*x*(-10:10));
lambda = (1/(2*pi)*A'*f);


a = asin(sqrt(2)*1i);
fa = 1./(2+sin(x-a*1i).^2);
M = max(fa)
errbound = 4*pi*M/a*(exp(-a*10))

figure
semilogy(-10:10,lambda,'*')
hold on
semilogy(-10:10,errbound*ones(length(-10:10)),'r')


