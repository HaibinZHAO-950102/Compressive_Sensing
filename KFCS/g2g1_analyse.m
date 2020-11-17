clc
clear
close all

printfigure = 1;
z = 0 : 0.001 : 10;

y = 0.01 * z;

y11 = 0.0001 *(z-3).^2;
y12 = 0.0001 *(z - 10).^2;
y1 = y11 + y12;
plot(z,y,'r-','linewidth',5)
xlim([0 10])
ylim([0 0.5])
setplt('Erste Term','$\Delta z$','','term1',printfigure)

figure
plot(z,y1,'b-','linewidth',5)
xlim([0 10])
ylim([0 0.1])
setplt('Zweiter Term','$\Delta z$','','term2',printfigure)


figure
for g2g1 = 1 : 0.5 : 50
    clf
    y2 = g2g1 * y1;
    Y = y + y2;
    [a,b] = min(Y);
    plot(z, Y,'k-','linewidth',5)
    hold on
    plot([z(b),z(b)],[0,a],'k--','linewidth',2);
    if z(b) < 3.5
        xticks([ z(b),3.5])
        xticklabels({'$\Delta \hat{z}$','$\Delta z ^*$'})
        yticks([])
        yticklabels({})
    elseif z(b) == 3.5
        xticks([3.5])
        xticklabels({'$\Delta \hat{z} = \Delta z ^*$'})
                yticks([])
        yticklabels({})

    else
        xticks([ 3.5,z(b)])
        xticklabels({'$\Delta z ^*$', '$\Delta \hat{z}$'})
                yticks([])
        yticklabels({})

    end
    f1 = gca;
    f1.TickLabelInterpreter = 'latex';



    xlim([0 10])
    ylim([0 0.5])
    setplt('Werte','$\Delta z$','','term3',printfigure)
    drawnow
    frame=getframe(gcf);
    imind=frame2im(frame);
    [imind,cm] = rgb2ind(imind,256);
    name = ['g2g1_analyse.gif'];
    if printfigure == 1
        if g2g1 == 1
             imwrite(imind,cm,name,'gif', 'Loopcount',inf,'DelayTime',1e-6);
        else
             imwrite(imind,cm,name,'gif','WriteMode','append','DelayTime',1e-6);
        end
    end
end

close all