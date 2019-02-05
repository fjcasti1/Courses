clear all
close all
clc

x = chebfun('x',[0 2*pi]);
f = 1./(2+sin(x).^2);
plot(f)

 
% N=10;
% A = exp(1i*x*(-N:N));
% I = A*A^(-1)
% I(1,1)
% cccc
N = 20:10:150;
for i=1:length(N)
    A = exp(1i*x*(-N(i):N(i)));
    fhat_k = (1/2*pi)*A'*f;
    fN = A*fhat_k;
    err(i) = max(abs(f-fN));
    i
end
semilogy(N,err,'*')
%%

close all; clear variables

%% Problem 4
close all
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
semilogy(err,'*')

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
lambdabound = 2*pi*M*exp(-a*N)

figure
loglog(-N:N,lambda,'*')
hold on
loglog(-N:N,errbound*ones(length(-N:N)),'r')
loglog(-N:N,lambdabound*ones(length(-N:N)),'g')

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
close all
clear all
clc

x = chebfun('x',[0 2*pi]);
f = 1./(2+sin(x).^2);
A = exp(1i*x*(-10:10));
lambda = (1/(2*pi)*A'*f);

fn = A*lambda;
err = norm(f-fn,2)


zroot = asin(sqrt(2)*1i);
a = imag(zroot)
fa = 1./(2+sin(x-a).^2);
M = max(fa)
errbound = 4*pi*M/a*(exp(-a*abs((-10:10))))';
lambdabound = 2*pi*M*exp(-a*abs((-10:10)))';

figure
semilogy(-10:10,lambda,'*')
hold on
semilogy(-10:10,errbound,'r')
semilogy(-10:10,lambdabound,'g')

