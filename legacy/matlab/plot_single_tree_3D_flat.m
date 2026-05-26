%%
clear all
file = 'ZZZ1r0001gen';
figure(11);clf
set(gcf,'Color','white','Position',[200 500 902 802])
imovie = 1;
clear M

rgbcolor = [
 14  74 106 
103 167 208
185 226 251
248 208 190
224 128 105
173  12   1
]/255;

%rgbcolor = flipud(rgbcolor);
%rgbcolor = .5*ones(6,3);

%%
for istep = 00:1:250
    filestep = strcat('ZStat_',file,num2str(istep,'%#08d'),'.dat');
    DATA = load (filestep);
    Nbranches = size(DATA,1);
    Age = max(DATA(:,1));
    
    figure(11);cla
    mesh(30*[-1 1],30*[-1 1]',10*[0 0; 0 0], ...
        'FaceColor',.9*[1 1 1],'EdgeColor','w');hold on
    %plot3(30,30,18,'.w')
    
    
    ct=cos(.1*[0 1 2 3 4 5 6 7 8 9 10]*2*pi);
    st=sin(.1*[0 1 2 3 4 5 6 7 8 9 10]*2*pi);
    [xs,ys,zs] = sphere;
    
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
        generation  = DATA(i,1);
        diameter    = DATA(i,2);
        length      = DATA(i,3);
        t           = DATA(i,4:6);
        n           = DATA(i,7:9);
        location    = DATA(i,10:12);
        
        b = cross(t,n);
        b = b / norm(b);
        X = (location(1) + length*[0;t(1)])*ones(1,11) + ...
            0.5*diameter*[1;1]*n(1)*ct + ...
            0.5*diameter*[1;1]*b(1)*st;
        Y = (location(2) + length*[0;t(2)])*ones(1,11) + ...
            0.5*diameter*[1;1]*n(2)*ct + ...
            0.5*diameter*[1;1]*b(2)*st;
        Z = (location(3) + length*[0;t(3)])*ones(1,11) + ...
            0.5*diameter*[1;1]*n(3)*ct + ...
            0.5*diameter*[1;1]*b(3)*st;
        surf(X,Y,.01*ones(size(Z)),'EdgeColor','none','FaceLighting','flat', ...
            'FaceColor',.1*[1 1 1] );
    end
    
    
    for j = 1:Nbranches
        i = iY(j);
        generation  = DATA(i,1);
        diameter    = DATA(i,2);
        length      = DATA(i,3);
        t           = DATA(i,4:6);
        n           = DATA(i,7:9);
        location    = DATA(i,10:12);
        
        b = cross(t,n);
        b = b / norm(b);
        X = (location(1) + length*[0;t(1)])*ones(1,11) + ...
            0.5*diameter*[1;1]*n(1)*ct + ...
            0.5*diameter*[1;1]*b(1)*st;
        Y = (location(2) + length*[0;t(2)])*ones(1,11) + ...
            0.5*diameter*[1;1]*n(2)*ct + ...
            0.5*diameter*[1;1]*b(2)*st;
        Z = (location(3) + length*[0;t(3)])*ones(1,11) + ...
            0.5*diameter*[1;1]*n(3)*ct + ...
            0.5*diameter*[1;1]*b(3)*st;
        surf(X,Y,Z,'EdgeColor','none','FaceLighting','phong', ...
            'FaceColor',rgbcolor(generation,:) );
    end
    
    material dull
    %     camlight(30,80)
    %     camlight(10,30)
    axis equal
    xlabel('x'); ylabel('y')
    axis off
    
    view([180 14])
    h=camlight(0,76,'infinite')
    set(h,'Color',0.75*[1 1 1])
    h=camlight('headlight','infinite')
    %set(h,'Color',0.75*[1 1 1])
    
    title(strcat('Generation: ', num2str(istep,'%#03d'), ...
        '     Nb branches: ',num2str(Nbranches,'%#04d')))

end