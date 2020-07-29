function L=ErrorNorm2D(v,v_app,k)
    [M,N]=size(v);
    if k==inf
        L=max(max(abs(v-v_app)));
    elseif k==1
        L=sum(sum(abs(v-v_app)))/(M*N);
    elseif k==2
        L=sqrt(sum(sum((v-v_app).^2))/(M*N));
    end
end