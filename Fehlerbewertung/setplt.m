function setplt(figurename,xname,yname,filename,printfigure)
    set(gca,'Fontsize',60)
    set(gca,'fontname','times new Roman')
    T = title(figurename,'fontsize',60);
    set(T,'Interpreter','latex')
    T = xlabel(xname,'fontsize',60);
    set(T,'Interpreter','latex')
    T = ylabel(yname,'fontsize',60);
    set(T,'Interpreter','latex')
    set(gcf,'outerposition',get(0,'screensize'));
    if printfigure == 1
        print(filename,'-dpng','-r600')
    end
    drawnow
end
