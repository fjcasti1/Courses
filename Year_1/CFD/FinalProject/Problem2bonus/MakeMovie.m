% Make movie

if boolMovie
MOV_u=struct('cdata',[],'colormap',[]);
MOV_v=struct('cdata',[],'colormap',[]);
MOV_Y=struct('cdata',[],'colormap',[]);
MOV_mix=struct('cdata',[],'colormap',[]);
end
if boolMovie
[MOV_u,MOV_v,MOV_Y,MOV_mix]=getMovie(u,v,Y,M,N,hx,hy,MOV_u,MOV_v,MOV_Y,MOV_mix,tstep);
end
if boolMovie
[MOV_u,MOV_v,MOV_Y,MOV_mix]=getMovie(u,v,Y,M,N,hx,hy,MOV_u,MOV_v,MOV_Y,MOV_mix,tstep);
end

%% Save movies
    if boolMovie
    %Movie for u
    myVideo=VideoWriter([pathMOVIES,'uMovie.avi']);
    myVideo.FrameRate=12;
    open(myVideo)
    writeVideo(myVideo,MOV_u);
    close(myVideo)
    %Movie for u
    myVideo=VideoWriter([pathMOVIES,'vMovie.avi']);
    myVideo.FrameRate=12;
    open(myVideo)
    writeVideo(myVideo,MOV_v);
    close(myVideo)
    %Movie for Y
    myVideo=VideoWriter([pathMOVIES,'YMovie.avi']);
    myVideo.FrameRate=12;
    open(myVideo)
    writeVideo(myVideo,MOV_Y);
    close(myVideo)
    %Movie for mix Y(1-Y)
    myVideo=VideoWriter([pathMOVIES,'YmixMovie.avi']);
    myVideo.FrameRate=12;
    open(myVideo)
    writeVideo(myVideo,MOV_mix);
    close(myVideo)
    end