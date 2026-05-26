%%
file = 'ZZZ1r0001gen';
figure(11);clf
set(gcf,'Color','white','Position',[700 800 1402 702])
scale= 1.0

%%
for istep = 0:1000:1000
    filestep = strcat('ZStat2_',file,num2str(istep,'%#08d'),'.dat');
    DATA = load (filestep);
    Nbranches = size(DATA,1);
    Age = max(DATA(:,1));
    
    figure(11);cla
    plot([0 0],scale*[-2 20],'w-')
    plot(scale*[-24 24],[0 0], 'k-')
    hold on
    C=colormap('copper');
%     C=colormap('jet');
    
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
        distanceG   = DATA(i,13);
        t           = DATA(i,4:6);
        distanceL   = DATA(i,1);
        Xmin        = location(1);
        Ymin        = location(2);
        Zmin        = location(3);
        Xmax        = location(1) + length * t(1);
        Ymax        = location(2) + length * t(2);
        Zmax        = location(3) + length * t(3);
        plot([Xmin Xmax]-scale*12, [Zmin Zmax], '-k','LineWidth', ...
            max(diameter,0.01)*10,'Color',C(ceil(64*(1.0001-distanceL/Age)),:))
    end
       
    figure(11)
    axis equal
    axis off
    title(strcat('Generation: ', num2str(istep,    '%#03d'), ...
        '     Nb branches: ',num2str(Nbranches,'%#04d')))
    
    %%
    distanceL   = DATA(:,1);
    distanceG   = DATA(:,13);
    area        = 0.25*pi*DATA(:,2).^2;
    
    maxN = ceil(max(distanceL));
    N = zeros(1,maxN);
    area_tot = N;
    for i = 1:Nbranches
        j = ceil(distanceL(i));
        N(j) = N(j) + 1;
        area_tot(j) = area_tot(j) + area(i);
    end
    figure(12)
    subplot(2,1,1)
    loglog((1:maxN),N,'ro',(1:maxN),area_tot./N,'ko')
    xlabel('mean distance to leaves')
    legend('Number of branches','mean section')
    subplot(2,1,2)
    plot((1:maxN),area_tot,'ko')
    xlabel('mean distance to leaves')
    ylabel('Total area')
    
    %%
    filestep = strcat('ZAreaRatio_',file,num2str(istep,'%#08d'),'.dat');
    DATA = load (filestep);
    Nnodes = size(DATA,1);
    if Nnodes>0
        distanceL   = DATA(:,1);
        distanceG   = DATA(:,2);
        area0       = DATA(:,3);
        area1       = DATA(:,4);
        area2       = DATA(:,5);
        area_ratio  = (area1+area2)./area0;
        
        figure(13)
        semilogy(distanceL,area_ratio,'ko')
        xlabel('mean distance to leaves')
        ylabel('Area ratio: r')
        disp(['mean area ratio: ', ...
        num2str(mean(area_ratio(1:sum(distanceL>1.55)))), ...
        '  // STD area ratio:',...
        num2str(std(area_ratio(1:sum(distanceL>1.55)))) ])
    end
    %%
    pause(0.1)
    
end