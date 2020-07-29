function phi = updateGhostCells(phi,t,M)
% Left boundary, Dirichlet BC
phi(3)=2*phi_left(t)-phi(4); 
phi(2)=2*phi_left(t)-phi(5);
phi(1)=2*phi_left(t)-phi(6);
% Right boundary, zero Neumann BC
phi(M+4)=phi(M+3);
phi(M+5)=phi(M+2);
phi(M+6)=phi(M+1);
end