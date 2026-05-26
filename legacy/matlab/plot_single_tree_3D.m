%%
clear all
file = 'ZZZ1r0001gen';
figure(11);clf
set(gcf,'Color','white','Position',[200 500 902 802])
imovie = 1;
clear M

%%
for istep = 0:1:250
    filestep = strcat(file,num2str(istep,'%#08d'),'.dat');
    DATA = load (filestep);
    Nbranches = size(DATA,1);
    Age = max(DATA(:,1));
    
    figure(11);cla
    rgbcolor=colormap(copper);
    mesh(10*[-1 1],10*[-1 1]',10*[0 0; 0 0], ...
        'FaceColor',[0.7 1 0.7],'EdgeColor','k');hold on
    plot3(30,30,18,'.w')
    axis equal
    xlabel('x'); ylabel('y')
    axis off
    set(gca,'Projection','orthographic', ...
        'CameraViewAngle',35, ...
        'CameraPosition',[-35 -35 30], ...
        'CameraTarget', [0 0 14])
    
    ct=cos(.1*[0 1 2 3 4 5 6 7 8 9 10]*2*pi);
    st=sin(.1*[0 1 2 3 4 5 6 7 8 9 10]*2*pi);
    [xs,ys,zs] = sphere;
    
    for i = 1:Nbranches,
        generation  = DATA(i,1);
        diameter    = DATA(i,2);
        length      = DATA(i,3);
        t           = DATA(i,4:6);
        n           = DATA(i,7:9);
        location    = DATA(i,10:12);
        light       = DATA(i,13);
        reserve     = DATA(i,14);
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
            'FaceColor',rgbcolor(ceil(16 + 32*generation/Age),:));
%         if light > 1E-8
%             surf (location(1) + length*t(1) + 0.5*xs, ...
%                 location(2) + length*t(2) + 0.5*ys, ...
%                 location(3) + length*t(3) + 0.5*zs, ...
%                 'EdgeColor','none','FaceLighting','phong', ...
%                 'FaceColor',[.1 1 0]*(.9*max(min(light,1),0)+.1))
%         end
    end
    
    %material dull
%     camlight(30,80)
%     camlight(10,30)
    title(strcat('Generation: ', num2str(istep,'%#03d'), ...
        '     Nb branches: ',num2str(Nbranches,'%#04d')))
    pause(0.1)
    h=figure(11);
    M(imovie) = getframe(h,[50 -410 800 700]);
    imovie = imovie + 1;
end