function phi = applyBCs(phi,M,N,hx,hy,opt)
switch opt
    case 'u'
        y=linspace(-hy/2,2+hy/2,N+2);
        %---------- Boundary Conditions for u----------%
        % Left and right
        for j=1:N+2
            if (y(j)>=0.5 && y(j)<=1) % Inlet 1
                phi(1,j)=-48*y(j)^2+72*y(j)-24;
            elseif (y(j)>=1 && y(j)<=1.5) % Inlet 2
                phi(M+1,j)=24*y(j)^2-60*y(j)+36;
            end
        end
        % Top and bottom
        phi(:,1)=-phi(:,2); 
        phi(:,N+2)=-phi(:,N+1);
	case 'v'
        x=linspace(-hx/2,4+hx/2,M+2);
        %---------- Boundary Conditions ----------%
        % Top and bottom
        for i=1:M+2
            if (x(i)>=0.5 && x(i)<=1) % Inlet 3
                phi(i,N+1)=24*x(i)^2-36*x(i)+12;
            elseif (x(i)>=1.5 && x(i)<=2.5) % Outlet
                phi(i,1)=(4*phi(i,2)-phi(i,3))/3;
            end
        end
        % Left and right
        phi(1,:)=-phi(2,:); 
        phi(M+2,:)=-phi(M+1,:);
    case 'Y'
        x=linspace(-5*hx/2,4+5*hx/2,M+6);
        y=linspace(-5*hy/2,2+5*hy/2,N+6);
        %---------- Boundary Conditions ----------%
        % Zero Neumann for walls and outlet
        %Left
        phi(3,4:N+3)=phi(4,4:N+3); 
        phi(2,4:N+3)=phi(5,4:N+3); 
        phi(1,4:N+3)=phi(6,4:N+3); 
        %Right
        phi(M+4,4:N+3)=phi(M+3,4:N+3); 
        phi(M+5,4:N+3)=phi(M+2,4:N+3); 
        phi(M+6,4:N+3)=phi(M+1,4:N+3); 
        %Bottom
        phi(4:M+3,3)=phi(4:M+3,4);
        phi(4:M+3,2)=phi(4:M+3,5);
        phi(4:M+3,1)=phi(4:M+3,6);
        %Top
        phi(4:M+3,N+4)=phi(4:M+3,N+3);
        phi(4:M+3,N+5)=phi(4:M+3,N+2);
        phi(4:M+3,N+6)=phi(4:M+3,N+1);
        % Dirichlet for the Inlets
        %Inlets 1 and 2
        for j=3:N+3
            if (y(j)>=0.5 && y(j)<=1) % Inlet 1
                phi(3,j)=2-phi(4,j);
                phi(2,j)=2-phi(5,j);
                phi(1,j)=2-phi(6,j);
            elseif (y(j)>=1 && y(j)<=1.5) % Inlet 2
                phi(M+4,j)=2*0.25-phi(M+3,j);
                phi(M+5,j)=2*0.25-phi(M+2,j);
                phi(M+6,j)=2*0.25-phi(M+1,j);
            end
        end
        %Inlet 3
        for i=3:M+3
            if (x(i)>=0.5 && x(i)<=1) % Inlet 3
                phi(i,N+4)=-phi(i,N+3);
                phi(i,N+5)=-phi(i,N+2);
                phi(i,N+6)=-phi(i,N+1);
            end
        end
end
end