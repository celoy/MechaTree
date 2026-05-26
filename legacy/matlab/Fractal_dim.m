%%
file = 'ZZZ1r0001gen';
figure(1);clf
set(gcf,'Color','white','Position',[700 800 1402 702])

istep = 1000;

%%
filestep = strcat('ZStat_',file,num2str(istep,'%#08d'),'.dat');
DATA = load (filestep);
Nbranches = size(DATA,1);
Age = max(DATA(:,1));

figure(1);cla
plot([0 0],scale*[-2 20],'w-')
plot(scale*[-10 10],[0 0], 'k-')
hold on
C=colormap('jet');

locationX = zeros(Nbranches,1);
locationY = zeros(Nbranches,1);
for i = 1:Nbranches
    locationX(i) = DATA(i,10);
    locationY(i) = DATA(i,10);
end
[~, iX] = sort(locationX);
[~, iY] = sort(locationY,'descend');

for j = 1:Nbranches
    i = iY(j);
    length      = DATA(i,3);
    diameter    = DATA(i,2);
    location    = DATA(i,10:12);
    reserve     = DATA(i,13);
    t           = DATA(i,4:6);
    Strahler    = DATA(i,1);
    Xmin        = location(1);
    Ymin        = location(2);
    Zmin        = location(3);
    Xmax        = location(1) + length * t(1);
    Ymax        = location(2) + length * t(2);
    Zmax        = location(3) + length * t(3);
    plot([Xmin Xmax], [Zmin Zmax], '-k','LineWidth', ...
        max(diameter,.001)*10,'Color',C(ceil(64*Strahler/Age),:))
end

axis equal
axis off
title(strcat('Generation: ', num2str(istep,    '%#03d'), ...
    '     Nb branches: ',num2str(Nbranches,'%#04d')))

%%
Strahlermax = max(DATA(:,1));
area = 0.25*pi*DATA(:,2).^2;
area_tot = zeros(1,Strahlermax);
for i = 1:Nbranches
    j = ceil(DATA(i,1));
    area_tot(j) = area_tot(j) + area(i);
end

%%
filelength  = strcat('Z_length_',file,'.dat');
DATA = load(filelength);
gen = DATA(:,1);
genmax = istep;
lengths = DATA(:,2:Strahlermax + 1);
figure(4); clf
semilogy(lengths(genmax,:),'ok'); hold on

filesegment = strcat('Z_Nsegments_',file,'.dat');
DATA = load(filesegment);
Nsegments = DATA(:,2:Strahlermax + 1);
semilogy(Nsegments(genmax,:),'+k')



semilogy(area_tot./Nsegments(genmax,1:Strahlermax),'vk')


%%

xlabel('Strahler order')
legend('Mean length','Number of branches','Mean area','Location','NorthEast')