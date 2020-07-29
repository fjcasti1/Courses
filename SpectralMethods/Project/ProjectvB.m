%% Spectral Methods Project - Francisco Castillo
clear all; close all; clc;
format long;

%% TO DO %% 
% - RHS of the equations
% - Decide and implement observables
% - Check Poisson equation in reserch, missing a minus??? 
% - Probably much more
% ----------------------------------------------------------------------------%
%% TO IMPROVE %% 
% - Implement different Nr, Nz
% ----------------------------------------------------------------------------%
%% TO ASK %% 
% - Sense of the vertical direction. Seems weird, show rr zz for N=10.
% - Correct approach about the 1/r factor?
% ----------------------------------------------------------------------------%

%%
% Global variables
global mgraph;
global m;
global Nsteps;
global N;

% Parameters
N = 30;  % Number of Grid Points
Re = 3000;
dt = 1e-10;
 
% Grid and diff matrices
[D,xch] = cheb(N-1);  
r = -xch; D = -D;   % To trasnlate the domain to [-1,1]
[rr,zz] = meshgrid(r);

Dp = D';
D2 = D^2; D2p = D2';
indb = find(rr==-1|rr==1|zz==0|zz==1); % For what?? Impose BCs in velocities

% laplacian
A = D2(2:end-1,2:end-1);
B = D2p(2:end-1,2:end-1)-Dp(2:end-1,2:end-1)*diag(1./r(2:end-1));
L = kron(eye(N-2),A)+kron(B',eye(N-2));  %Laplacian of the interior
%Linv = inv(L);
[lo,up,per] = lu(L,'vector');  % LU factorization

% Initial velocity & pre-allocate memory
psi = 0*rr;
g = psi;
w = psi;
gp = g;
wp = w;
ur = 0*rr;
uw = ur;
uz = uw;
% Everything initialized to zero
count = 0;
TU = 1;
Nsteps = TU/dt;
mgraph = 1;

for m=1:Nsteps
    
    % Calculate vorcity and angular momentum
    w = D*ur-uz*Dp;
    g = uw*diag(r);
    
%     g
%     w
%     psi
%     keyboard

    check(w,g,psi)
    
    % Prediction, interior
    [rhsg,rhsw] = RHS(r,Re,g,w,D,D2,Dp,D2p,ur,uz); % Right hand side of PDE
    gp(2:end-1,2:end-1) = g(2:end-1,2:end-1) + dt*rhsg(2:end-1,2:end-1);
    wp(2:end-1,2:end-1) = w(2:end-1,2:end-1) + dt*rhsw(2:end-1,2:end-1);

    % Solve Streamfunction, interior
    psi(2:end-1,2:end-1) = sylvester(D2(2:end-1,2:end-1),D2p(2:end-1,2:end-1)-Dp(2:end-1,2:end-1)*diag(1./r(2:end-1)),-wp(2:end-1,2:end-1))*diag(r(2:end-1));
%     C = -wp(2:end-1,2:end-1)*diag(r(2:end-1)); C=C(:);
%     psi(2:end-1,2:end-1) = reshape(up\(lo\(C(per))),N-2,N-2);

    check(wp,gp,psi)
    
    % Calculate velocities
    ur = -(D*psi)*diag(1./r);
    uw = gp*diag(1./r);
    uz = (psi*Dp)*diag(1./r);
    % Boundary Conditions imposed directly on the velocities (so far updated the interior)
    [ur,uw,uz] = BCS(ur,uw,uz);
    
    % Calculate vorcity and angular momentum
    wp = D*ur-uz*Dp;
    gp = uw*diag(r);
    
% %     gp
% %     wp
% %     psi
% %     keyboard
    
    % Correction, interior
    [rhsg,rhsw] = RHS(r,Re,gp,wp,D,D2,Dp,D2p,ur,uz); % Right hand side of PDE
    g(2:end-1,2:end-1) = 0.5*(g(2:end-1,2:end-1)+gp(2:end-1,2:end-1)+dt*rhsg(2:end-1,2:end-1));
    w(2:end-1,2:end-1) = 0.5*(w(2:end-1,2:end-1)+wp(2:end-1,2:end-1)+dt*rhsw(2:end-1,2:end-1));
    
    % Solve Streamfunction, interior
    psi(2:end-1,2:end-1) = sylvester(D2(2:end-1,2:end-1),D2p(2:end-1,2:end-1)-Dp(2:end-1,2:end-1)*diag(1./r(2:end-1)),-w(2:end-1,2:end-1))*diag(r(2:end-1));
%     C = -w(2:end-1,2:end-1)*diag(r(2:end-1)); C=C(:);
%     psi(2:end-1,2:end-1) = reshape(up\(lo\(C(per))),N-2,N-2);

    check(w,g,psi)    
    
    % Calculate velocities
    ur = -(D*psi)*diag(1./r);
    uw = g*diag(1./r);
    uz = (psi*Dp)*diag(1./r);
    % Boundary Conditions imposed directly on the velocities (so far updated the interior)
    [ur,uw,uz] = BCS(ur,uw,uz);
    
    % Observables
    if (mod(m,mgraph)==0)
        plotting(rr,zz,w,g,psi,ur,uw,uz)
        
    end
    end

function [rhsg,rhsw] = RHS(r,Re,g,w,D,D2,Dp,D2p,ur,uz)
    
% % %     rhsg = ((D*psi).*(g*Dp)-(psi*Dp).*(D*g))*diag(1./r)+(1/Re)*(D2*g+g*D2p-g*Dp*diag(1./r));
% % %     rhsw = ((D*psi).*(w*Dp)-(psi*Dp).*(D*w))*diag(1./r)-w.*(D*psi)*diag(1./r.^2)...
% % %         +2*g.*(D*g)*diag(1./r.^3)+(1/Re)*(D2*w+w*D2p+w*Dp*diag(1./r)-w*diag(1./r.^2));
    
    rhsg = -ur.*(g*Dp)-uz.*(D*g)+(1/Re)*(D2*g+g*D2p-g*Dp*diag(1./r));
    rhsw = -ur.*(w*Dp)-uz.*(D*w)+ur.*w*diag(1./r)...
        +2*g.*(D*g*diag(1./abs(r.^3)))+(1/Re)*(D2*w+w*D2p+w*Dp*diag(1./r)-w*diag(1./r.^2));
%     keyboard
%     if m == 10
%         keyboard
%     end
% 	if (mod(m,mgraph)==0)
%         rhsg
%         rhsw
%         m
%         keyboard
%     end
%     if (mod(m,100)==0)
%         keyboard
%     end
end

function [ur,uw,uz] = BCS(ur,uw,uz) % Boundary Conditions imposed directly on the velocities
    global N;
    % Outer Radius (Left & Right): no slip boundary condition
        ur(:,[1 end]) = 0;
        uw(:,[1 end]) = 0;
        uz(:,[1 end]) = 0;
    % Bottom: no slip boundary condition
        ur(1,:) = 0;
        uw(1,:) = 0;
        uz(1,:) = 0;
    % Top: angular velocity = 1
        ur(end,:) = 0;
        uw(end,1:N/2) = -1;
        uw(end,N/2+1:end) = 1;
        uz(end,:) = 0;
end
function check(w,g,psi)
    global m;
    if(max(max(isnan(w))) || max(max(isinf(w))))
        m
        keyboard
        ccccc
    elseif(max(max(isnan(g))) || max(max(isinf(g))))
        m
        keyboard
        ccccc
    elseif(max(max(isnan(psi))) || max(max(isinf(psi))))
        m
        keyboard
        ccccc
    end
end

function plotting(rr,zz,w,g,psi,ur,uw,uz)
    global m;   
    global Nsteps;
    fprintf('Progress: %d %% \n',m*100/Nsteps)

    if max(max(isnan(w)))==1
        m
        ccc
    elseif max(max(isnan(w)))==1
        m
        ccc
    elseif max(max(isnan(w)))==1
        m
        ccc
    end

    subplot(2,3,1)
    contourf(rr,zz,psi)
    axis([-1 1 -1 1]), axis square
    title('Streamlines','fontsize',16)

    subplot(2,3,2)
    contourf(rr,zz,w)
    axis([-1 1 -1 1]), axis square
    title('Vorticity','fontsize',16)

    subplot(2,3,3)
    contourf(rr,zz,g)
    axis([-1 1 -1 1]), axis square
    title('Ang Momentum','fontsize',16)
    
    speed = sqrt(ur.^2+uw.^2+uz.^2);
       
	subplot(2,3,4)
	quiver3(rr,0*zz,zz,ur./speed,uw./speed,uz./speed)
    axis([-1 1 -1 1]), axis square
    title('Velocity','fontsize',16)
    view(-35,45)
    hold on
    [X,Y,Z] = cylinder;
    Z(1,:) = -1;
    testsubject = surf(X,Y,Z); 
    set(testsubject,'FaceAlpha',0.25)
    set(testsubject,'FaceAlpha',0.25,'EdgeColor','red','EdgeAlpha',0.4,...
    'DiffuseStrength',1,'AmbientStrength',1)
    hold off
    
    subplot(2,3,5)
	quiver3(rr,0*zz,zz,ur./speed,uw./speed,uz./speed)
    axis([-1 1 -1 1]), axis square
    title('Velocity','fontsize',16)
    view(0,90)
    hold on
    [X,Y,Z] = cylinder;
    Z(1,:) = -1;
    testsubject = surf(X,Y,Z); 
    set(testsubject,'FaceAlpha',0.25)
    set(testsubject,'FaceAlpha',0.25,'EdgeColor','red','EdgeAlpha',0.4,...
    'DiffuseStrength',1,'AmbientStrength',1)
    hold off
    
    subplot(2,3,6)
	quiver3(rr,0*zz,zz,ur./speed,uw./speed,uz./speed)
    axis([-1 1 -1 1]), axis square
    title('Velocity','fontsize',16)
    view(0,0)
%     hold on
%     [X,Y,Z] = cylinder;
%     Z(1,:) = -1;
%     testsubject = surf(X,Y,Z); 
%     set(testsubject,'FaceAlpha',0.25)
%     set(testsubject,'FaceAlpha',0.25,'EdgeColor','red','EdgeAlpha',0.4,...
%     'DiffuseStrength',1,'AmbientStrength',1)
%     hold off

    drawnow
    pause
end