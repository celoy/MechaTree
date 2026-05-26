%%
file = 'YYY_arch00021000.dat';
datas = load(file);

angles   = datas(:,1:3);
NNbranch = datas(:,4:13);
NNreserve= datas(:,14:31);
data     = datas(:,41:42);
data(:,3)= datas(:,40);
Nindividuals = size(data,1);

[maxx, ix] = max(data(:,1));
[maxy, iy] = max(data(:,2));
[maxz, iz] = max(data(:,3));

disp(strcat('>> max Moment leaves:', num2str(ix),      ' / max: ', num2str(maxx)))
disp(strcat('max N seeds:', num2str(iy), ' / max: ', num2str(maxy) ))

%%
N=128;
nb_leaves  = zeros(N,N);
max_stress = zeros(N,N);
Nleaves    = zeros(N,N);
vol_relat  = zeros(N,N);

Safety     = zeros(N,N);
Pleaves    = zeros(N,N);
Pseeds     = zeros(N,N);
Phototro   = zeros(N,N);

for k = 1:2
    figure(k)
    clf
    set(gcf,'Color','white','Position',[k*100 k*100+650 800 800])
    colormap(winter)
    
    if k==1
        ii = ix;
        disp('Best Moment leaves')
    elseif k==2
        ii = iy;
        disp('Best N seeds')
    end
    disp(strcat('all perfs: ',num2str(data(ii,:))))
    
    for i = 1:N
        for j = 1:N
            nb_leaves(i,j)  = 50*(i-1)/(N-1);
            max_stress(i,j) = 2*(j-1)/(N-1);
            Nleaves(i,j)    = 100*(i-1)/(N-1);
            vol_relat(i,j)  = 0.0079 * (j-1)/(N-1);
            
            Safety(i,j) = neural_branch(nb_leaves(i,j), max_stress(i,j), NNbranch(ii,:));
            [Pseeds(i,j), Pleaves(i,j), Phototro(i,j)] = neural_reserve(Nleaves(i,j),vol_relat(i,j), NNreserve(ii,:));
        end
    end
    
    subplot(2,2,1)
    %     [C,h] = contour(nb_leaves, max_stress, Safety, (1:.2:3));
    [C,h] = contourf(nb_leaves, max_stress, Safety); colorbar
    clabel(C,h);
    xlabel('N_{leaves}')
    ylabel('max_(stress)')
    title('Safety')
    
    subplot(2,2,2)
    [C,h] = contourf(Nleaves, vol_relat, Pseeds); colorbar
    clabel(C,h);
    xlabel('N_{leaves}')
    ylabel('Vol relative')
    title('% seeds')
    
    subplot(2,2,3)
    [C,h] = contourf(Nleaves, vol_relat, Pleaves); colorbar
    clabel(C,h);
    xlabel('N_{leaves}')
    ylabel('Vol relative')
    title('% leaves')
    
    subplot(2,2,4)
    [C,h] = contourf(Nleaves, vol_relat, Phototro); colorbar
    clabel(C,h);
    xlabel('N_{leaves}')
    ylabel('Vol relative')
    title('Phototropism')
    Phototro
end

