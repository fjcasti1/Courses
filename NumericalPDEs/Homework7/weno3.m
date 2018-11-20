function weno3(rhol,vl,pl,rhor,vr,pr,gamma,tf,N,cfl)
%weno3(rhol,vl,pl,rhor,vr,pr,gamma,tf,N,cfl) solves 1D gas dynamical 
% Riemann Problems using 3rd-order WENO with Lax-Friedrichs flux splitting. 
%   BCs are through-flow.  The exact RP solution is superimposed.
%       rhol,vl,pl =  left RP parameters: density, velocity, pressure
%       rhor,vr,pr = right RP parameters: density, velocity, pressure
%       gamma = polytropic gas constant
%       tf = final time
%       N = number of dx
%       cfl = CFL factor
%   Examples:
%       1. weno3(5,30,0.4127,0.5,0,0.4127,5/3,0.015,500,0.5) Jet RP
%           Ha & Gardner, J Sci Comp 34 (2008) 247
%       2. weno3(1,0,1,0.125,0,0.1,1.4,0.2,200,0.9) Sod RP
%   From Liska & Wendroff, SIAM J Sci Comp 25 (2003) 995:
%       3. weno3(1,0.75,1,0.125,0,0.1,1.4,0.2,200,0.9) Toro/Sod RP
%       4. weno3(1,-19.59745,1000,1,-19.59745,0.01,1.4,0.01,200,0.9)
%           Stationary Contact RP (Toro)  
%       5. weno3(1,-2,0.4,1,2,0.4,1.4,0.1,400,0.5) Near-Vacuum RP (Toro)
%       6. weno3(1,1,10^-6,1,-1,10^-6,5/3,1,200,0.75) Noh RP
%       7. weno3(0.1261192,8.9047029,782.92899,6.591493,...
%           2.2654207,3.1544874,1.4,0.0055,800,0.9) Peak RP (Kucharik)
%       8. weno3(5.99924,19.5975,460.894,5.99242,-6.19633,46.095,1.4,...
%           0.02,400,0.9) Two Strong Shocks RP (Toro)

%   Developed by Jeremiah Jones and Carl Gardner, Arizona State University
%       Email comments to carl.gardner@asu.edu
%   WENO implementation based on a Fortran code of Guan-Shan Jiang and
%   Chi-Wang Shu

%   NOTES:
%       0. To speed the code up significantly, comment out the movie!
%       1. The MATLAB code is optimized for readability
%       2. 5th-order WENO gives better resolution of the RPs (good student
%          project!)
%       3. If transformation to characteristic variables space is not used
%          with WENO fluxes, the Jet, Stationary Contact, and Noh RP 
%          examples above fail, while the Peak RP example gives p and v 
%          overshoots and oscillations (good HW problem!)
%       4. For difficult problems, a positivity-preserving version of WENO
%          should be used
%       5. The computed internal energy is not as accurate as rho, v, p

%   Reverse RP Examples for Testing Code:
%       weno3(0.125,0,0.1,1,0,1,1.4,0.2,200,0.9) reverse Sod RP
%       weno3(0.125,0,0.1,1,-0.75,1,1.4,0.2,200,0.9) reverse Toro/Sod RP
%       weno3(0.5,0,0.4127,5,-30,0.4127,5/3,0.015,500,0.5) reverse Jet RP

%   The Euler equations of 1D gas dynamics are written as conservation laws
%       q_t + f_x = 0
%   which are then approximated by 
%       (dq/dt)(j) = -(F(j+1/2)-F(j-1/2))/dx
%   Conserved quantities solution array
%       q(j,eq#) 
%   is an (N+1)x3 matrix:
%       rho = q(:,1) = mass density
%         m = q(:,2) = momentum density
%         E = q(:,3) = energy density

% if N is odd, reset N to N+1
if mod(N,2) 
    N = N + 1;
    fprintf('N reset = %g to be even\n',N);
end

dx = 1/N;
x = (-0.5:dx:0.5)'; % grid points

% exact RP solution [rho,v,p]:
qexact = exact_RP_sol(x,rhol,vl,pl,rhor,vr,pr,gamma,tf);

tic;

% allocate memory for solution array
q = zeros(N+1,3); % 3 conserved quantities for 1D gas dynamics

% left RP parameters
ml = rhol*vl;
El = pl/(gamma - 1) + 0.5*rhol*vl^2;

% right RP parameters
mr = rhor*vr;
Er = pr/(gamma - 1) + 0.5*rhor*vr^2;

% set up initial conditions
q(1:N/2,1) = rhol;
q(1:N/2,2) = ml;
q(1:N/2,3) = El;

q(N/2+1:N+1,1) = (rhol + rhor)/2;
q(N/2+1:N+1,2) = (ml + mr)/2;
q(N/2+1:N+1,3) = (El + Er)/2;

q(N/2+2:N+1,1) = rhor;
q(N/2+2:N+1,2) = mr;
q(N/2+2:N+1,3) = Er;

% time loop
figure % for movie
n = 0;
t = 0;
while t < tf
    n = n + 1;
    p = pressure(q,gamma);
    % compute dt 
    alpha = max_speed(q,p,gamma);
    dt = cfl*dx/alpha;
    if 1.01*dt >= tf - t
        dt = tf - t;
    end
    t = t + dt;
    
    % (TVD) RK3 timestep for dq/dt = rhs = -df/dx
    q1 = q + dt*rhs(dx,q,p,gamma,alpha);
    p1 = pressure(q1,gamma);
    alpha1 = max_speed(q1,p1,gamma);
    q2 = 0.75*q + 0.25*(q1 + dt*rhs(dx,q1,p1,gamma,alpha1));
    p2 = pressure(q2,gamma);
    alpha2 = max_speed(q2,p2,gamma);
    q = (q + 2*(q2 + dt*rhs(dx,q2,p2,gamma,alpha2)))/3;
    
    % movie
    plot(x,q(:,1),'b-','LineWidth',2); title('Density');
    % plot(x,q(:,2)./q(:,1),'c-','LineWidth',2); title('Velocity'); 
    % p = pressure(q,gamma); plot(x,p,'m-','LineWidth',2); title('Pressure');
    set(gca,'fontsize',24);
    xticks([-0.5 -0.25 0 0.25 0.5]);
    getframe;
    % M(n) = getframe;
end

figure;
plot(x,qexact(:,1),'r-',x,q(:,1),'b-','LineWidth',1);
set(gca,'fontsize',24);
xticks([-0.5 -0.25 0 0.25 0.5]);
title('Density');
% legend('exact','WENO','Location','North');

figure;
plot(x,qexact(:,2),'r-',x,q(:,2)./q(:,1),'b-','LineWidth',1);
set(gca,'fontsize',24);
xticks([-0.5 -0.25 0 0.25 0.5]);
title('Velocity'); 
% legend('exact','WENO','Location','South');

figure; % since ~ p/rho, computational errors are compounded
p = pressure(q,gamma);
plot(x,qexact(:,3)./((gamma - 1)*qexact(:,1)),'r-',...
    x,p./((gamma - 1)*q(:,1)),'b-','LineWidth',1);
set(gca,'fontsize',24);
xticks([-0.5 -0.25 0 0.25 0.5]);
title('Internal Energy = T/(\gamma-1)');
% legend('exact','WENO','Location','South');

figure;
p = pressure(q,gamma);
plot(x,qexact(:,3),'r-',x,p,'b-','LineWidth',1);
set(gca,'fontsize',24);
xticks([-0.5 -0.25 0 0.25 0.5]);
title('Pressure');
legend('exact','WENO','Location','South');

toc;

end

function f = flux(q,p) 
% exact flux for 1D gas dynamics
v = q(:,2)./q(:,1);
f(:,1) = q(:,2);
f(:,2) = q(:,1).*v.^2 + p;
f(:,3) = v.*(q(:,3) + p);   
end

function p = pressure(q,gamma)
p = (gamma - 1)*(q(:,3) - 0.5*q(:,2).^2./q(:,1));
if (numel(find(p < 0)) > 0)
    error('pressure < 0 computed in WENO: try decreasing cfl and/or dx');
end 
end

function alpha = max_speed(q,p,gamma)
% computes maximum characteristic speed
v = q(:,2)./q(:,1);
c = sqrt(gamma*p./q(:,1));
alpha = max(abs(v) + c);
end

function RHS = rhs(dx,q,p,gamma,alpha)
N = size(q,1) - 1;
F = num_flux(q,p,gamma,alpha);
RHS = -(F(2:N+2,:) - F(1:N+1,:))/dx;  
end

function F = num_flux(q,p,gamma,alpha)
% WENO method
N = size(q,1) - 1;

% through-flow BCs with ghost points
qghost = [q(1,:);q(1,:);q;q(N+1,:);q(N+1,:)];
pghost = [p(1);p(1);p;p(N+1);p(N+1)];

f = flux(qghost,pghost);

% Lax-Friedrichs flux splitting
% there are N+5 elements of fplus and fminus
fplus = 0.5*(f + alpha*qghost); 
fminus = 0.5*(f - alpha*qghost);

% compute left and right Roe matrices
L = zeros(3,3,N+2);
R = zeros(3,3,N+2);
for j = 1:N+2
    [LL,RR] = roe_matrix(qghost(j+1:j+2,:),pghost(j+1:j+2),gamma);
    L(:,:,j) = LL;
    R(:,:,j) = RR;
end

% compute positive WENO flux 
                        % F(j+1/2) and F(j-1/2) stencils
fjm1 = fplus(1:N+2,:);  % f(j-1)       f(j-2)
fj = fplus(2:N+3,:);    % f(j)         f(j-1)
fjp1 = fplus(3:N+4,:);  % f(j+1)       f(j)
% transform to characteristic variables space
for j = 1:N+2
   fj(j,:) = fj(j,:)*L(:,:,j)';
   fjp1(j,:) = fjp1(j,:)*L(:,:,j)';
   fjm1(j,:) = fjm1(j,:)*L(:,:,j)';
end
Fplus = weno3_flux(fjm1,fj,fjp1);

% compute negative WENO flux
                        % F(j+1/2) and F(j-1/2) stencils
fj = fminus(2:N+3,:);   % f(j)         f(j-1)
fjp1 = fminus(3:N+4,:); % f(j+1)       f(j)
fjp2 = fminus(4:N+5,:); % f(j+2)       f(j+1)
% transform to characteristic variables space
for j = 1:N+2
   fjp2(j,:) = fjp2(j,:)*L(:,:,j)';
   fjp1(j,:) = fjp1(j,:)*L(:,:,j)';
   fj(j,:) = fj(j,:)*L(:,:,j)';
end
Fminus = weno3_flux(fjp2,fjp1,fj);

% compute WENO flux in characteristic variables space
F = Fplus + Fminus;

% transform WENO flux back to physical space
for j = 1:N+2
   F(j,:) = F(j,:)*R(:,:,j)';
end
end

function [L,R] = roe_matrix(q,p,gamma)
% q and p have only 2 rows j and j+1: 2x3 matrix
% computes left and right Roe matrices at x(j+1/2) based on Roe average of
% two states q(j) and q(j+1)                  

v = q(:,2)./q(:,1);
c = sqrt(gamma*p./q(:,1));
h = (q(:,3) + p)./q(:,1); % enthalpy
w = sqrt(q(:,1));

% mean values
t = w(1)/(w(1) + w(2));
vm = t*v(1) + (1 - t)*v(2);
cm = t*c(1) + (1 - t)*c(2);
hm = t*h(1) + (1 - t)*h(2);

% define Roe matrices with mean values
R = [1, 1, 1; vm - cm, vm, vm + cm; hm - cm*vm, 0.5*vm^2, hm + cm*vm];
% L = R\eye(3);
b = (gamma - 1)/cm^2;
L = [0.5*(b*0.5*vm^2 + vm/cm), -0.5*(b*vm + 1/cm), 0.5*b;
     1 - b*0.5*vm^2, b*vm, -b;
     0.5*(b*0.5*vm^2 - vm/cm), -0.5*(b*vm - 1/cm) ,0.5*b];
end

function F = weno3_flux(fjm1,fj,fjp1)
epsilon = 10^-6;
c1 = 1/3; 
c2 = 2/3; 

% compute WENO flux 
beta1 = (fj - fjm1).^2;        % left smoothness indicator
beta2 = (fj - fjp1).^2;        % right smoothness indicator
w1 = c1./(epsilon + beta1).^2; % non-scaled left weight
w2 = c2./(epsilon + beta2).^2; % non-scaled right weight 
sum = w1 + w2;                 % sum of non-scaled weights
w1 = w1./sum;                  % scaled left weight
w2 = w2./sum;                  % scaled right weight
f1 = 0.5*(3*fj - fjm1);        % numerical flux on left stencil 
f2 = 0.5*(fj + fjp1);          % numerical flux on right stencil 
F = w1.*f1 + w2.*f2;           % WENO flux    
end

function q = exact_RP_sol(x,rhol,vl,pl,rhor,vr,pr,gamma,t)
% q = [rho,v,p](j)
SHOCK = 1;
RAREFACTION = 2;
N = size(x,1) - 1;
q = zeros(N+1,3);

[pstar,vstar,Ml,Mr,lwave,rwave] = ...
    exact_RP(rhol,vl,pl,rhor,vr,pr,gamma,SHOCK,RAREFACTION);
fprintf('pstar = %g, vstar = %g\n',pstar,vstar);
if lwave == SHOCK
    fprintf('left wave = SHOCK\n');
elseif lwave == RAREFACTION
    fprintf('left wave = RAREFACTION\n');
end
if rwave == SHOCK
    fprintf('right wave = SHOCK\n');
elseif rwave == RAREFACTION
    fprintf('right wave = RAREFACTION\n');
end

xc = vstar*t; % location of contact
for j = 1:N+1
    if x(j) >= xc
        if rwave == SHOCK
            vs = vr + Mr/rhor;
            if x(j) > vs*t
                q(j,:) = [rhor,vr,pr];
            else
                rhostar = Mr/(vs - vstar);
                q(j,:) = [rhostar,vstar,pstar];
            end
        else % rwave == RAREFACTION
            cr = sqrt(gamma*pr/rhor);
            cstar = cr + 0.5*(vstar - vr)*(gamma - 1);
            x1 = (vstar + cstar)*t; % left boundary of rarefaction
            x2 = (vr + cr)*t; % right boundary of rarefaction
            if x(j) < x1
                rhostar = gamma*pstar/cstar^2;
                q(j,:) = [rhostar,vstar,pstar];
            elseif x(j) > x2
                q(j,:) = [rhor,vr,pr];
            else % inside rarefaction
                vj = x(j)/t;
                cbar = (2*cr + (gamma - 1)*(vj - vr))/(1 + gamma);
                vbar = vj - cbar;
                rhobar = (cbar^2*rhor^gamma/(gamma*pr))^(1/(gamma - 1));
                pbar = cbar^2*rhobar/gamma;
                q(j,:) = [rhobar,vbar,pbar];
            end
        end
    else % x(j) < xc
        if lwave == SHOCK
            vs = vl - Ml/rhol;
            if x(j) < vs*t
                q(j,:) = [rhol,vl,pl];
            else
                rhostar = Ml/(vstar - vs);
                q(j,:) = [rhostar,vstar,pstar];
            end
        else % lwave == RAREFACTION
            cl = sqrt(gamma*pl/rhol);
            cstar = cl + 0.5*(vl - vstar)*(gamma - 1);
            x1 = (vl - cl)*t; % left boundary of rarefaction
            x2 = (vstar - cstar)*t; % right boundary of rarefaction 
            if x(j) < x1
                 q(j,:) = [rhol,vl,pl];
            elseif x(j) > x2
                rhostar = gamma*pstar/cstar^2;
                q(j,:) = [rhostar,vstar,pstar];
            else % inside rarefaction
                vj = x(j)/t;
                cbar = (2*cl + (gamma - 1)*(vl - vj))/(1 + gamma);
                vbar = vj + cbar;
                rhobar = (cbar^2*rhol^gamma/(gamma*pl))^(1/(gamma - 1));
                pbar = cbar^2*rhobar/gamma;
                q(j,:) = [rhobar,vbar,pbar];
            end            
        end
    end
end

end

function [pstar,vstar,Ml,Mr,lwave,rwave] = ...
    exact_RP(rhol,vl,pl,rhor,vr,pr,gamma,SHOCK,RAREFACTION)
% based on Chorin, J Comp Phys 22 (1976) 517, modification of Godunov
MAX_ITER_alpha = 50;
MAX_ITER = 50;
EPS = 8*eps;
MIN_PRESSURE = realmin;

% initial guess
pstar = 0.5*(pl + pr);
coefl = sqrt(pl*rhol);
coefr = sqrt(pr*rhor);
Ml = coefl*phi(pstar/pl,gamma);
Mr = coefr*phi(pstar/pr,gamma);

% iterate at least twice to avoid spurious convergence when pl = pr    
pstar = ((vl - vr)*Ml*Mr + pr*Ml + pl*Mr)/(Ml + Mr);
if pstar < MIN_PRESSURE 
    pstar = MIN_PRESSURE;
end
Ml = coefl*phi(pstar/pl,gamma);
Mr = coefr*phi(pstar/pr,gamma);

% iterative loop using Godunov's method, as adapted by Chorin
n = 1;
alpha = 2;
for n1 = 1:MAX_ITER_alpha
alpha = alpha/2;
for n2 = 1:MAX_ITER
    n = n + 1;
    Ml_old = Ml;
    Mr_old = Mr;  
    ptilde = ((vl - vr)*Ml*Mr + pr*Ml + pl*Mr)/(Ml + Mr);
    if ptilde <= MIN_PRESSURE 
        ptilde = MIN_PRESSURE;
    end
    pstar = alpha*ptilde + (1 - alpha)*(pstar);
    Ml = coefl*phi(pstar/pl,gamma);
    Mr = coefr*phi(pstar/pr,gamma);
    if abs(Ml - Ml_old) <= EPS*Ml && abs(Mr - Mr_old) <= EPS*Mr
    % if abs(Ml - Ml_old) <= EPS && abs(Mr - Mr_old) <= EPS
        vstar = (pl - pr + Mr*vr + Ml*vl)/(Ml + Mr);
        if pstar >= pl      
            lwave = SHOCK;
        else
            lwave = RAREFACTION;
        end
        if pstar >= pr       
            rwave = SHOCK;
        else
            rwave = RAREFACTION;
        end
        fprintf('exact RP solution converged in %g iterations\n',n);
        if (pstar == MIN_PRESSURE)
            fprintf('pstar = MIN_PRESSURE: likely error in exact RP\n');
        end
        return;
    end
end
end

error('exact RP solution failed to converge in %g iterations\n',n);
end

function PHI = phi(r,gamma)
EPS = 8*eps;

if abs(1 - r) <= EPS
    PHI = sqrt(gamma);
elseif r < 1
    beta = 0.5*(gamma - 1)/gamma;
    PHI = 0.5*(gamma - 1)*(1 - r)/(sqrt(gamma)*(1 - r^beta));
else % r > 1
    PHI = sqrt(0.5*(gamma + 1)*r + 0.5*(gamma - 1));
end
end