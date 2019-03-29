clear all
%%% matlab code to generate a spatio-temporal sources along a fault of the overthrust model


fid=fopen('win_overt.bin','r');
A=fread(fid,[160 450],'real*4');
fclose(fid);


%%% create 16 sources between (231,45) to (216, 60)

faultx=linspace(216,231,16);
faultz=linspace(60,45,16);

tpeak = linspace(0.2,0.8,16);
dt=0.002; ns=3000;
t=(0:1:ns-1).*dt;
fc=30.0;
t0=1/fc;

for i=1:1:16
    gn = wgn(ns,1,10); %% white gaussian noise
    gn = gn./max(abs(gn)); %% normalize
    
    %%% filter between 1 to 80 Hz
    [N,Fo,Ao,W] = firpmord([0 1 60 80],[0 1 0],[0.1 0.01 0.1],1/dt);
    fp=firpm(N,Fo,Ao,W);
    snois(:,i) = 0.5.*gn;
    
    ts=tpeak(i);
    tmp = pi*fc.*(t-ts-t0).^2;
	wt = (1.0-2.0.*tmp).*exp(-tmp);
    swlt(:,i) = wt.'+0.5.*gn;
    
end

nx=450; nz=160;
src3d=zeros(ns*nx*nz,1);

for is=1:1:length(faultx)
    wlt=swlt(:,is);
    
    for it=1:1:ns
        ind=it+ns*faultx(is)+ns*nx*faultz(is);
        src3d(ind)=wlt(it);
    end
end

fid=fopen('source3d.bin','w');
fwrite(fid,src3d,'real*4');
fclose(fid);
