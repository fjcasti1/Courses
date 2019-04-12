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
N = 40;  %40
Re = 1200;
dt = 1e-6;
 
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
% Everything initialized to zero
count = 0;
TU = 10;
Nsteps = TU/dt;

for m=1:Nsteps
    % Prediction, interior
%     [rhsg,rhsw] = RHS(m,r(2:end-1),Re,psi(2:end-1,2:end-1),g(2:end-1,2:end-1),...
%         w(2:end-1,2:end-1),D(2:end-1,2:end-1),D2(2:end-1,2:end-1),...
%         Dp(2:end-1,2:end-1),D2p(2:end-1,2:end-1));
    [rhsg,rhsw] = RHS(m,r,Re,psi,g,w,D,D2,Dp,D2p);
%     keyboard
    gp(2:end-1,2:end-1) = g(2:end-1,2:end-1) + dt*rhsg(2:end-1,2:end-1);
    wp(2:end-1,2:end-1) = w(2:end-1,2:end-1) + dt*rhsw(2:end-1,2:end-1);
    % Solve Streamfunction, interior
    psi(2:end-1,2:end-1) = sylvester(D2(2:end-1,2:end-1),D2p(2:end-1,2:end-1)-Dp(2:end-1,2:end-1)*diag(1./r(2:end-1)),-wp(2:end-1,2:end-1))*diag(r(2:end-1));

    % Boundary Conditions (so far updated the interior)
        % Axis: psi, g, w were initialized to zero and remain unchanged
	    % Outer Radius: psi, g were initialized to zero and remain unchanged
    wp(:,end) = -(1/r(end))*psi*D2p(:,end);
        % Bottom: psi, g were initialized to zero and remain unchanged
	wp(1,2:end-1) = -(D2(1,:)*psi(:,2:end-1))./r(2:end-1)';
        % Top: psi was initialized to zero and remain unchanged
    wp(end,2:end-1) = -(D2(end,:)*psi(:,2:end-1))./r(2:end-1)';
    gp(end,:) = r;
%     wp
    
    
    % Correction, interior
%     [rhsg,rhsw] = RHS(m,r(2:end-1),Re,psi(2:end-1,2:end-1),gp(2:end-1,2:end-1),...
%         wp(2:end-1,2:end-1),D(2:end-1,2:end-1),D2(2:end-1,2:end-1),...
%         Dp(2:end-1,2:end-1),D2p(2:end-1,2:end-1));
    [rhsg,rhsw] = RHS(m,r,Re,psi,gp,wp,D,D2,Dp,D2p);
    g(2:end-1,2:end-1) = 0.5*(g(2:end-1,2:end-1)+gp(2:end-1,2:end-1)+dt*rhsg(2:end-1,2:end-1));
    w(2:end-1,2:end-1) = 0.5*(w(2:end-1,2:end-1)+wp(2:end-1,2:end-1)+dt*rhsw(2:end-1,2:end-1));
    % Solve Streamfunction, interior
    psi(2:end-1,2:end-1) = sylvester(D2(2:end-1,2:end-1),D2p(2:end-1,2:end-1)-Dp(2:end-1,2:end-1)*diag(1./r(2:end-1)),-w(2:end-1,2:end-1))*diag(r(2:end-1));
    % Boundary Conditions (so far updated the interior)
        % Axis: psi, g, w were initialized to zero and remain unchanged
	    % Outer Radius: psi, g were initialized to zero and remain unchanged
    w(:,end) = -(1/r(end))*psi*D2p(:,end);
        % Bottom: psi, g were initialized to zero and remain unchanged
	w(1,2:end-1) = -(D2(1,:)*psi(:,2:end-1))./r(2:end-1)';
        % Top: psi was initialized to zero and remain unchanged
    w(end,2:end-1) = -(D2(end,:)*psi(:,2:end-1))./r(2:end-1)';
    g(end,:) = r;
    % Observables
    m
    if (mod(m,0000)==0)
        psi
        g
        w
        keyboard
    end
end
function [rhsg,rhsw] = RHS(m,r,Re,psi,g,w,D,D2,Dp,D2p)
    rhsg = ((D*psi).*(g*Dp)-(psi*Dp).*(D*g))*diag(1./r)+(1/Re)*(D2*g+g*D2p-g*Dp*diag(1./r));
    rhsw = ((D*psi).*(w*Dp)-(psi*Dp).*(D*w))*diag(1./r)-w.*(D*psi)*diag(1./r.^2)...
        -2*g.*(D*g)*diag(1./r.^3)+(1/Re)*(D2*w+w*D2p+w*Dp*diag(1./r)-w*diag(1./r.^2));

%     if (mod(m,100)==0)
%         keyboard
%     end
end