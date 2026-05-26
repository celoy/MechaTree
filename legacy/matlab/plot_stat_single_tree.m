%%
file = 'ZZZ1r0001gen';
figure(11);clf
set(gcf,'Color','white','Position',[700 800 1402 702])
scale= 1.0

%%
filelength  = strcat('Z_length_',file,'.dat');
DATA = load(filelength);
gen = DATA(:,1);
genmax = size(DATA,1);
genmax = 190;
Strahlermax=sum(sum(DATA > 0) > 0) - 1;
lengths = DATA(:,2:Strahlermax + 1);
figure(4); clf
%semilogy(lengths(1:genmax,:)','o-')
semilogy(lengths(genmax,:),'ok')
xlabel('Strahler order')
ylabel('Mean length')

filesegment = strcat('Z_Nsegments_',file,'.dat');
DATA = load(filesegment);
Nsegments = DATA(:,2:Strahlermax + 1);
figure(5); clf
semilogy(Nsegments(genmax,:),'ok')
xlabel('Strahler order')
ylabel('Number of branches')


fileTokuna  = strcat('Z_Tokunaga_',file,'.dat');
DATA = load(fileTokuna);
sizeT = size(DATA,2);
Tokunaga = DATA(genmax,2:sizeT);
Tokunaga = reshape(Tokunaga,sqrt(sizeT-1),sqrt(sizeT-1));
Tokunaga = Tokunaga(1:Strahlermax-1,1:Strahlermax-1)
figure(6); clf
rgbcolor=colormap(cool);
semilogy(0,1,'-y')
hold on
for i=2:Strahlermax-1;
    plot((2:i)-1,fliplr(Tokunaga(i,1:i-1)),'o-k','Color',rgbcolor(ceil(64*(i-1)/(Strahlermax-1)),:))
end
xlabel('Strahler order')
ylabel('Number of branches')

%%
for istep = 0:100:1000
    filestep = strcat('ZStat_',file,num2str(istep,'%#08d'),'.dat');
    DATA = load (filestep);
    Nbranches = size(DATA,1);
    Age = max(DATA(:,1));
    
    figure(11);cla
    plot([0 0],scale*[-2 20],'w-')
    plot(scale*[-24 24],[0 0], 'k-')
    hold on
    C=colormap('copper');
    C=colormap('jet');
    
    locationX = zeros(Nbranches,1);
    locationY = zeros(Nbranches,1);    
    for i = 1:Nbranches
        locationX(i) = DATA(i,10);
        locationY(i) = DATA(i,10);
    end
    [toto, iX] = sort(locationX);
    [toto, iY] = sort(locationY,'descend');
       
    for j = 1:Nbranches
        i = iY(j);
        length      = DATA(i,3);
        diameter    = DATA(i,2);
        location    = DATA(i,10:12);
        reserve     = DATA(i,13);
        t           = DATA(i,4:6);
        generation  = DATA(i,1);
        Xmin        = location(1);
        Ymin        = location(2);
        Zmin        = location(3);
        Xmax        = location(1) + length * t(1);
        Ymax        = location(2) + length * t(2);
        Zmax        = location(3) + length * t(3);
        plot([Xmin Xmax]-scale*12, [Zmin Zmax], '-k','LineWidth', ...
            max(diameter,.001)*10,'Color',C(ceil(64*generation/Age),:))
    end
    
    axis equal
    axis off
    plot(log(reserve)*[-1 1],[-1 -1],'b-')
    d_trunk    = DATA(1,2);
    title(strcat('Generation: ', num2str(istep,    '%#03d'), ...
            '     Nb branches: ',num2str(Nbranches,'%#04d')))
    pause(0.1)
end