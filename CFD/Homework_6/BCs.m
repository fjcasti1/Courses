function BC=BCs(i,j,hx,hy,Z)
xi=(i-1)*hx;
yj=(j-1)*hy;
opt=inputname(5);
switch opt
    case 'u'
        M=size(Z,1)-1;
        N=size(Z,2)-2;
        if (i==1 && yj>0.5 && yj<=1)% Left Boundary
            BC=2;
        elseif (i==M+1 && yj>1 && yj<=1.5)% Right Boundary
            BC=-1;
        elseif (j==N+2)% Top Boundary
            BC=-Z(i,N+1);
        elseif (j==1)% Bottom Boundary
            BC=-Z(i,2);
        else % No slid Boundary Condition for the walls
            BC=0;
        end
        
    case 'v'
        if (xi==0 && yj>=0.5 && yj<=1) %v_in1
            BC = 0;
        elseif (xi==4 && yj>=1 && yj<=1.5) %v_in2
            BC = 0;
        elseif (yj==2 && xi>=0.5 && xi<=1) %v_in3
            BC = -1;
        elseif (yj==0 && xi>=1.5 && xi<=2) %v_out
            BC = 0;
            Neumann_vout
        else % no-slip BC
            BC = 0;
        end
    case 'Y'
        if (xi==0 && yj>=0.5 && yj<=1) %Y_in1
            BC = 1;
        elseif (xi==4 && yj>=1 && yj<=1.5) %Y_in2
            BC = 0.25;
        elseif (yj==2 && xi>=0.5 && xi<=1) %Y_in3
            BC = 0;
        elseif (yj==0 && xi>=1.5 && xi<=2) %Y_out
            BC = 0;
            Neumann_Yout
        else % no-slip BC
            BC = 0;
            Neumann_Yin
        end
end
end