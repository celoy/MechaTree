
SizeForest = 100;
skip = 5;

%%
figure(10);clf
set(gcf,'Color','white','Position',[950 100 802 802])
figure(11);clf
set(gcf,'Color','white','Position',[10 80 902 802])

%%
for generation = 00:500:10000,
    file = 'ZZZ';
    DATA = load( strcat(file,num2str(generation,'%#08d'),'.dat') );
    Ntrees = size(DATA,1);
    genes = zeros(Ntrees,3);
    DATA = sortrows(DATA,-6);
    datas = DATA(1:skip:Ntrees,7:46);
    N2 = size(datas,1);
    NNbranch = datas(:,4:13);
    NNreserve= datas(:,14:31);
    
    %     genes(:,1) = max(0.0, min(1.0, .5*DATA(:,7) - .5*DATA(:,8) + .5 ));
    %     genes(:,2:3) = DATA(:,9:10);
    genes = DATA(:,41:1:43);      % 41:46 junk DNA
    %     genes(:,1) = genes(:,1)  / min(1.0, max(genes(:,1))*2) ;
    %     genes(:,1) = ( genes(:,1) - min(genes(:,1)) ) / ( max(genes(:,1)) - min(genes(:,1)) );
    %     genes(:,2) = ( genes(:,2) - min(genes(:,2)) ) / ( max(genes(:,2)) - min(genes(:,2)) );
    %     genes(:,3) = ( genes(:,3) - min(genes(:,3)) ) / ( max(genes(:,3)) - min(genes(:,3)) );
    
    %%
    nb_leaves   = 5;
    max_stress  = 1;
    Nleaves     = 10;
    vol_relative= 0.5; 
    Safety   = zeros(N2,1);
    Phototro = zeros(N2,1);
    Pleaves  = zeros(N2,1);
    Pseeds   = zeros(N2,1);
    
    figure(11);
    subplot(2,2,1);
    scatter(90*DATA(1:skip:Ntrees,7),-90*DATA(1:skip:Ntrees,8),(DATA(1:skip:Ntrees,6) + 1)*4,genes(1:skip:Ntrees,:),'fill');
    xlabel('\theta_1');
    ylabel('\theta_2');
    axis([0 90 -90 0]);
    title (strcat('time: ', num2str(generation) ) );
    
    subplot(2,2,2);
    scatter(90*DATA(1:skip:Ntrees,7),360*DATA(1:skip:Ntrees,9),(DATA(1:skip:Ntrees,6) + 1)*4,genes(1:skip:Ntrees,:),'fill');
    xlabel('\theta_1');
    ylabel('\gamma');
    axis([0 90 0 360]);
    
    for i = 1:N2
        Safety(i) = neural_branch(nb_leaves, max_stress, NNbranch(i,:));
        [Pleaves(i) Pseeds(i) Phototro(i)] = neural_reserve(Nleaves,vol_relative, NNreserve(i,:));
    end
    
    subplot(2,2,3)
    scatter(Phototro,Safety,(DATA(1:skip:Ntrees,6) + 1)*4,genes(1:skip:Ntrees,:),'fill');
    xlabel('Phototropism')
    xlim([0 1])
    ylabel('Safety')
    ylim([1 5])
    title (strcat('generation: ', num2str(generation) ) )
    
    subplot(2,2,4)
    scatter(Pleaves,Pseeds,(DATA(1:skip:Ntrees,6) + 1)*4,genes(1:skip:Ntrees,:),'fill');
    xlabel('Pleaves')
    xlim([0 1.0])
    ylabel('Pseeds')
    ylim([-0.1 1.1])
    
    %%
    figure(10);cla;
    viscircles([0 0],SizeForest,'EdgeColor',[0 0 0]);
    hold on;
    viscircles([0 0],(SizeForest+20),'EdgeColor',[1 1 1]);
    
    axis equal;
    xlabel('x'); ylabel('y');
    axis off;
    
    %%
    figure(10);
    scatter(DATA(1:skip:Ntrees,3),DATA(1:skip:Ntrees,4),(DATA(1:skip:Ntrees,6) + 0.5)*60,genes(1:skip:Ntrees,:),'fill');
    %     for j = 1:skip:Ntrees
    %         i = iZ(j);
    %         x = DATA(i,1);
    %         y = DATA(i,2);
    %         X = DATA(i,3);
    %         Y = DATA(i,4);
    %         R = DATA(i,5) + .1;
    %         H = DATA(i,6) + 1;
    %         G = genes(i,:);
    %         plot([x X],[y Y],'-k','LineWidth',20*H/SizeForest,'Color',G(1:3))
    %         viscircles([X Y],R,'LineWidth',20*H/SizeForest,'EdgeColor',G(1:3))
    %     end
    title (strcat('time: ', num2str(generation) ) );
    
    %%
    pause(0.1)
end