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
% Parameters
N = 30;  %40
Re = 1000;
dt = 1e-12;
 
% Grid and diff matrices
[D,xch] = cheb(N-1); 
r = (-xch+1)/2; D = -2*D;   % To trasnlate the domain to [0,1]
[rr,zz] = meshgrid(r);

Dp = D';
D2 = D^2; D2p = D2';
indb = find(rr==0|rr==1|zz==0|zz==1); % For what?? Impose BCs in velocities

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
Nsteps = TU/dt
global mgraph;
mgraph = 1000;

for m=1:Nsteps
    
    % Calculate vorcity and angular momentum
    w = D*ur-uz*Dp;
    g = uw*diag(r);

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
%     if m==24
%         keyboard
%     end
    % Prediction, interior
    [rhsg,rhsw] = RHS(m,r,Re,psi,g,w,D,D2,Dp,D2p,ur,uw,uz); % Right hand side of PDE
    gp(2:end-1,2:end-1) = g(2:end-1,2:end-1) + dt*rhsg(2:end-1,2:end-1);
    wp(2:end-1,2:end-1) = w(2:end-1,2:end-1) + dt*rhsw(2:end-1,2:end-1);

    % Solve Streamfunction, interior
    psi(2:end-1,2:end-1) = sylvester(D2(2:end-1,2:end-1),D2p(2:end-1,2:end-1)-Dp(2:end-1,2:end-1)*diag(1./r(2:end-1)),-wp(2:end-1,2:end-1))*diag(r(2:end-1));

    
    if(max(max(isnan(wp))) || max(max(isinf(wp))))
        m
        keyboard
        ccccc
    elseif(max(max(isnan(gp))) || max(max(isinf(gp))))
        m
        keyboard
        ccccc
    elseif(max(max(isnan(psi))) || max(max(isinf(psi))))
        m
        keyboard
        ccccc
    end
    
    
    
    % Calculate velocities
    ur = -(D*psi)*diag(1./r);
    uw = gp*diag(1./r);
    uz = (psi*Dp)*diag(1./r);

    % Boundary Conditions imposed directly on the velocities (so far updated the interior)
    % Axis and Outer Radius: no slip boundary condition (the axis has no flow)
        ur(:,[1 end]) = 0;
        uw(:,[1 end]) = 0;
        uz(:,[1 end]) = 0;
    % Bottom: no slip boundary condition
        ur(1,:) = 0;
        uw(1,:) = 0;
        uz(1,:) = 0;
    % Top: angular velocity = 1
        ur(end,:) = 0;
        uw(end,:) = 1;
        uz(end,:) = 0;
    
    % Calculate vorcity and angular momentum
    wp = D*ur-uz*Dp;
    gp = uw*diag(r);

    if(max(max(isnan(wp))) || max(max(isinf(wp))))
        m
        keyboard
        ccccc
    elseif(max(max(isnan(gp))) || max(max(isinf(gp))))
        m
        keyboard
        ccccc
    elseif(max(max(isnan(psi))) || max(max(isinf(psi))))
        m
        keyboard
        ccccc
    end
    
    
    % Correction, interior
    [rhsg,rhsw] = RHS(m,r,Re,psi,gp,wp,D,D2,Dp,D2p,ur,uw,uz); % Right hand side of PDE
    g(2:end-1,2:end-1) = 0.5*(g(2:end-1,2:end-1)+gp(2:end-1,2:end-1)+dt*rhsg(2:end-1,2:end-1));
    w(2:end-1,2:end-1) = 0.5*(w(2:end-1,2:end-1)+wp(2:end-1,2:end-1)+dt*rhsw(2:end-1,2:end-1));
    
    % Solve Streamfunction, interior
    psi(2:end-1,2:end-1) = sylvester(D2(2:end-1,2:end-1),D2p(2:end-1,2:end-1)-Dp(2:end-1,2:end-1)*diag(1./r(2:end-1)),-w(2:end-1,2:end-1))*diag(r(2:end-1));

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
    
    % Calculate velocities
    ur = -(D*psi)*diag(1./r);
    uw = g*diag(1./r);
    uz = (psi*Dp)*diag(1./r);

    % Boundary Conditions imposed directly on the velocities (so far updated the interior)
    % Axis and Outer Radius: no slip boundary condition (the axis has no flow)
        ur(:,[1 end]) = 0;
        uw(:,[1 end]) = 0;
        uz(:,[1 end]) = 0;
    % Bottom: no slip boundary condition
        ur(1,:) = 0;
        uw(1,:) = 0;
        uz(1,:) = 0;
    % Top: angular velocity = 1
        ur(end,:) = 0;
        uw(end,:) = 1;
        uz(end,:) = 0;

    % Observables
    %     m
    if (mod(m,mgraph)==0)
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
        
        subplot(1,3,1)
        contourf(rr,zz,psi)
        axis([0 1 0 1]), axis square
        title('Streamlines','fontsize',16)

        subplot(1,3,2)
        contourf(rr,zz,w)
        axis([0 1 0 1]), axis square
        title('Vorticity','fontsize',16)

        subplot(1,3,3)
        contourf(rr,zz,g)
        axis([0 1 0 1]), axis square
        title('Ang Momentum','fontsize',16)

        drawnow

        %         psi
        %         g
        %         w
        %         keyboard
    end
end
function [rhsg,rhsw] = RHS(m,r,Re,psi,g,w,D,D2,Dp,D2p,ur,uw,uz)
    global mgraph;
    
% % %     rhsg = ((D*psi).*(g*Dp)-(psi*Dp).*(D*g))*diag(1./r)+(1/Re)*(D2*g+g*D2p-g*Dp*diag(1./r));
% % %     rhsw = ((D*psi).*(w*Dp)-(psi*Dp).*(D*w))*diag(1./r)-w.*(D*psi)*diag(1./r.^2)...
% % %         +2*g.*(D*g)*diag(1./r.^3)+(1/Re)*(D2*w+w*D2p+w*Dp*diag(1./r)-w*diag(1./r.^2));
    
    rhsg = -ur.*(g*Dp)-uz.*(D*g)+(1/Re)*(D2*g+g*D2p-g*Dp*diag(1./r));
    rhsw = -ur.*(w*Dp)-uz.*(D*w)+ur.*w*diag(1./r)...
        +2*g.*(D*g*diag(1./r.^3))+(1/Re)*(D2*w+w*D2p+w*Dp*diag(1./r)-w*diag(1./r.^2));
    
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