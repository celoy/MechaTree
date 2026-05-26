%%
file = 'S2.dat';
datas = load(file);
theta1 = 90*mean(datas(1:100,7))
theta2 =-90*mean(datas(1:100,8))
gamma  =360*mean(datas(1:100,9))
datas=datas(4,7:46);

angles   = datas(:,1:3);
NNbranch = datas(:,4:13);
NNreserve= datas(:,14:31);
ii=1;

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

figure(1)
clf
set(gcf,'Color','white','Position',[300 300+650 800 800])

for i = 1:N
    for j = 1:N
        nb_leaves(i,j)  = 100*(i-1)/(N-1);
        max_stress(i,j) = 0.5*(j-1)/(N-1);
        Nleaves(i,j)    = 100*(i-1)/(N-1);
        vol_relat(i,j)  = 0.5*(j-1)/(N-1);
        
        Safety(i,j) = neural_branch(nb_leaves(i,j), max_stress(i,j), NNbranch(ii,:));
        [Pseeds(i,j), Pleaves(i,j), Phototro(i,j)] = neural_reserve(Nleaves(i,j),vol_relat(i,j), NNreserve(ii,:));
    end
end

%% Figure
figure(1); clf
set(gcf,'Color','white','Position',[600 600 267*3.75 440])
colormap(parula)

subplot(2,3,1)
Safety(1,1)=0;
Safety(1,N)=4;
[C,h] = contour(nb_leaves, max_stress, Safety, (0:.2:4)); colorbar
%imagesc(nb_leaves(:,1)', max_stress(1,N:-1:1), Safety'); colorbar
%clabel(C,h);
% xlabel('# leaves')
% ylabel('max stress')
% title('Safety')

subplot(2,3,2)
loglog(nb_leaves(:,1),Safety(:,39),nb_leaves(:,1),Safety(:,52))
hold on
loglog((2:10),3.3*(2:10).^(2/30))
legend('relative stress: 0.15','relative stress: 0.20')
xlim([1 100])
ylim([2.5 5])
% xlabel('# leaves')
% ylabel('Safety')

subplot(2,3,4)
Pseeds(1,1)=0;
Pseeds(1,N)=1;
[C,h] = contour(Nleaves, vol_relat, Pseeds, (0:.1:1)); colorbar
%clabel(C,h);
% xlabel('Total # leaves')
% ylabel('Vol relative')
% title('% seeds')

subplot(2,3,5)
Pleaves(1,1)=0;
Pleaves(1,N)=1;
[C,h] = contour(Nleaves, vol_relat, Pleaves, (0:.1:1)); colorbar
%clabel(C,h);
% xlabel('N leaves')
% ylabel('Vol relative')
% title('% segments')

subplot(2,3,6)
Phototro(1,1)=0;
Phototro(1,N)=1;
[C,h] = contour(Nleaves, vol_relat, Phototro, (0:.1:1)); colorbar
%clabel(C,h);
%xlabel('N leaves')
%ylabel('Vol relative')
% title('Photosensitivity')
