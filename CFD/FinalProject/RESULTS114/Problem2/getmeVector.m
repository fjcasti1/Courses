function Vector=getmeVector(M,N,z)
opt=inputname(3);
switch opt
    case 'u'
        Vector=zeros((M+1)*(N+2),1);
        for l=1:M+1
            for k=1:N+2
                Vector((k-1)*(M+1)+l)=z(l,k);
            end
        end
	case 'v'
        Vector=zeros((M+2)*(N+1),1);
        for l=1:M+2
            for k=1:N+1
                Vector((k-1)*(M+2)+l)=z(l,k);
            end
        end
	case 'Y'
        Vector=zeros((M+6)*(N+6),1);
        for l=1:M+6
            for k=1:N+6
                Vector((k-1)*(M+6)+l)=z(l,k);
            end
        end
end