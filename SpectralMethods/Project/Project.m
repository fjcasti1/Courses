%% Spectral Methods Project - Francisco Castillo
clear all; close all; clc;
format long;

%% TO DO %% 
% - RHS of the equations
% - Decide and implement observables
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
N = 10;  %40
Re = 1200;
dt = 1e-4;
 
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
% Everything initialized to zero
count = 0;
TU = 10;
Nsteps = TU/dt;
for m=1:Nsteps
    % Prediction, interior
    [rhsg,rhsw] = RHS(r,Re,psi,g,w,D,D2,Dp,D2p);
    gp = g + dt*(1);
    wp = w + dt*(1);
    % Solve Streamfunction, interior
    % Boundary Conditions

    % Correction, interior
    g = 0.5*(g+gp+dt*(2));
    w = 0.5*(w+wp+dt*(2));
    % Solve Streamfunction, interior
    % Boundary Conditions
    
    % Observables
    m
end
function [rhsg,rhsw] = RHS(r,Re,psi,g,w,D,D2,Dp,D2p)
    rhsg = ((D*psi).*(g*Dp)-(psi*Dp).*(D*g))*diag(1./r)+(1/Re)*(D2*g+g*D2p-g*Dp*diag(1./r));
    rhsw = ((D*psi).*(w*Dp)-(psi*Dp).*(D*w))*diag(1./r)-w.*(D*psi)*diag(1./r.^2)...
        -2*g.*(D*g)*diag(1./r.^3)+(1/Re)*(D2*w+w*D2p+w*Dp*diag(1./r)-w*diag(1./r.^2));
end