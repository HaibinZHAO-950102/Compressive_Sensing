function setplt(figurename,xname,yname,filename,printfigure)
    set(gca,'Fontsize',30)
    set(gca,'fontname','times new Roman')
    T = title(figurename,'fontsize',80);
    set(T,'Interpreter','latex')
    T = xlabel(xname,'fontsize',40);
    set(T,'Interpreter','latex')
    T = ylabel(yname,'fontsize',40);
    set(T,'Interpreter','latex')
    set(gcf,'outerposition',get(0,'screensize'));
    
    if printfigure == 1
        name = [filename,'.png'];
        exportgraphics(gcf, name,'Resolution',200)
    end
    drawnow
end
