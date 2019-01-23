function L=ErrorNorm(v,v_app,k)
    N=length(v);
    if k==inf
        L=max(abs(v-v_app));
    elseif k==1
        L=sum(abs(v-v_app))/N;
    elseif k==2
        L=sqrt(sum((v-v_app).^2)/N);
    end
end