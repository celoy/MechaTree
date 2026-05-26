%%
data = load ('ZAllocation.dat');

t  = data(:,1);
Nb = data(:,2);
Nl = data(:,3);
W  = data(:,4);
Nt = data(:,5);
Ns = data(:,6);
Np  = data(:,7);
R  = data(:,8)/(0.25*pi*0.01);
% write(16,*) iout, t%n_branches, t%n_leaves, wind, Ntwigs,Nseeds,Npruned, t%Reserve

%%
figure(1);clf
set(gcf,'Color','white','Position',[700 800 1402 702])

semilogy(t,Nb,'k-'); hold on
semilogy(t,W,'r-')
semilogy(t,Ns,'-m')
semilogy(t,Np,'-oc')
semilogy(t,R,'--k')
xlabel('t')
legend('N branches','Wind','New seeds','N pruned','Reserve')

% semilogy(t,Nb,t,Nl,t,W,t,Nt,'-o',t,Ns,t,Np,'-o',t,R)
% xlabel('t')
% legend('N branches','N leaves','Wind','New twigs','New seeds','N pruned','Reserve')
