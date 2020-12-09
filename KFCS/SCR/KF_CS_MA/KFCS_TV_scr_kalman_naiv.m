% wavelet + kalman integrated

clc
clear
close all

printfigure = 1;

load('Messwerte')

x = linspace(0,10,1024);
nx = length(x);
p = floor(linspace(1,1024,64));
m = f(p,1:10:end);

nt = size(m,2);








for kkk = [1,7,10]
    
    % random select of sensors
number = 12;
S = zeros(nt,number);
S(:,1) = ceil(rand(nt,1)*64);
for t = 1 : nt
    k = 2;
    while k <= number
        temp = ceil(rand()*64);
        if abs(S(t,:) - temp) ~= 0
            S(t,k) = temp;
            k = k + 1;
        end
    end
    S(t,:) = sort(S(t,:));
end


% CS Basis
load basis_ml
WT = THETA^-1;

Phi_cs = zeros(nt, number, 64);
for t = 1 : nt
    for i = 1 : number
        Phi_cs(t,i,(S(t,i))) = 1;
    end
end

z = zeros(64, nt);
z(:,1) = WT * f(p,1);


% Kalman-modell

lambda = 0 : 50;
lambda = lambda * pi / 10;

Psi = zeros(nx, 51);  % Eigenfunktionen

Psi(:, 1) = sqrt(1 / 10);
for i = 2 : 51
    for n = 1 : nx
        Psi(n, i) = sqrt(2 / 10) * cos(lambda(i) * x(n));
    end
end

k = 0.1;

A = zeros(51);
for i = 1 : 51
    A(i,i) = 1 - 0.01 * k * lambda(i)^2;
end

Phi_KF = zeros(length(p),nx);
for i = 1 : length(p)
    Phi_KF(i,p(i)) = 1;
end

    

    T = double(zeros(51, nt));
    T(:,1) = double((Psi' * Phi_KF' * Phi_KF * Psi)^-1 * Psi' * Phi_KF'* THETA * z(:, 1));

    G2G1 = 10^(kkk-1);
    f_e_wt_kf_naiv = zeros(nx,nt);
    f_e_wt_kf_naiv(:,1) = Psi * T(:,1);



% KFCS
for t = 2 : nt
    t
    h = squeeze(Phi_cs(t,:,:)) * THETA;
    y = m(S(t,:),t);
    T(:,t-1) = double((Psi' * Phi_KF' * Phi_KF * Psi)^-1 * Psi' * Phi_KF' * THETA * z(:, t-1));
    e = 0.05;   
    cvx_begin quiet;
        variable a(64,1);
        minimize(norm(a,1) + G2G1 * ([y - h * (z(:,t-1) + a) ; A * T(:,t-1) - (Psi' * Phi_KF' * Phi_KF * Psi)^-1 * Psi' * Phi_KF' * THETA * (z(:,t-1) + a)])'*([y - h * (z(:,t-1) + a) ; A * T(:,t-1) - (Psi' * Phi_KF' * Phi_KF * Psi)^-1 * Psi' * Phi_KF' * THETA * (z(:,t-1) + a)]));
    cvx_end;
    
    z(:, t) = z(:,t-1) + a;
    T(:,t) = double((Psi' * Phi_KF' * Phi_KF * Psi)^-1 * Psi' * Phi_KF'  * THETA * z(:, t));
end

% for t = 2 : nt
%     f_e_wt_kf_naiv(:,t) = THETA * z(:,t);
% end
for t = 2 : nt
    f_e_wt_kf_naiv(:,t) = Psi * T(:,t);
end



figure
for n = 1 : 5 : nt
    clf
    nf = (n-1)*10 + 1;
    plot(x, f(:,nf),'k-','LineWidth',5)
    hold on
    plot(x, f_e_wt_kf_naiv(:,n),'c-','LineWidth',5)
    hold on
    plot(x(p(S(n,:))),f(p(S(n,:)),nf)','r.','Markersize',40)
    legend('Signal','geschätztes Signal','Messungen')
    xlim([0 10])
    ylim([-0.5 2.5])
    txt = ['$t = ',num2str((n-1)*0.1),'$'];
    TEXT = text(0.8,0.2,txt,'FontSize',60);
    set(TEXT,'Interpreter','latex')
    txt = ['$M = 12$'];
    TEXT = text(0.8,0.6,txt,'FontSize',60);
    set(TEXT,'Interpreter','latex')
    txt = ['$G2/G1 = 10^{',num2str(kkk-1),'}$'];
    TEXT = text(0.8,1,txt,'FontSize',60);
    set(TEXT,'Interpreter','latex')
    setplt('','$x$','$f$','',0)
    drawnow
    frame=getframe(gcf);
    imind=frame2im(frame);
    [imind,cm] = rgb2ind(imind,256);
    name = ['CS_wavelet_KF_naiv_',num2str(G2G1),'.gif'];
    if printfigure == 1
        if n==1
             imwrite(imind,cm,name,'gif', 'Loopcount',inf,'DelayTime',1e-6);
        else
             imwrite(imind,cm,name,'gif','WriteMode','append','DelayTime',1e-6);
        end
    end
end


t = 0 : 0.1:20;
for n = 0 : 5
    figure
    timef = 2000 / 5 * n + 1
    timefe = 200 / 5 * n + 1
    plot(x, f(:,timef),'k-','LineWidth',5)
    hold on
    plot(x, f_e_wt_kf_naiv(:,timefe),'c-','LineWidth',5)
    hold on
    plot(x(p_index(S(timefe,:))),m(S(timefe,:),timefe),'r.','Markersize',40)
    legend('Signal','geschätztes Signal','Messungen')

    xlim([0 10])
    ylim([-0.5 2.5])
    txt = ['$t = ',num2str(4*n),'$'];
    T = text(0.8,0.2,txt,'FontSize',60);
    set(T,'Interpreter','latex')
    txt = ['$M = ',num2str(number),'$'];
    T = text(0.8,0.6,txt,'FontSize',60);
    set(T,'Interpreter','latex')
        txt = ['$G2/G1 = 10^{',num2str(kkk-1),'}$'];
    TEXT = text(0.8,1,txt,'FontSize',60);
    set(TEXT,'Interpreter','latex')
    drawnow
    frame=getframe(gcf);
    imind=frame2im(frame);
    [imind,cm] = rgb2ind(imind,256);
    name = ['KFCS_naiv_',num2str(G2G1),'_shot_',num2str(n+1)];
    setplt('','$x$','$f$',name,printfigure)
end

% eval(['f_e_wt_kf_naiv_kkk',num2str(G2G1),' = f_e_wt_kf_naiv;']);
% 
% save(['f_e_wt_kf_naiv_kkk',num2str(G2G1),'.mat'],['f_e_wt_kf_naiv_kkk',num2str(G2G1)]);
% 
close all

end
