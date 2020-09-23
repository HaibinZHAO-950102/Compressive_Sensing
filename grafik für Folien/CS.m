% meine approximation

clc
clear
close all

printfigure = 1;

load('Messwerte')
load basis_ml
load f_e_scr_kf_fdm

x = 0 : dx : 10;
nx = length(x);
Dt = 0.1;
m = m(:,1:Dt/dt:end);
nt = size(m,2);

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



Phi = zeros(nt, number, 64);
for t = 1 : nt
    for i = 1 : number
        Phi(t,i,S(t,i)) = 1;
    end
end

z = zeros(64, nt);
z(:,1) = MT * f(p_index,1);

% for t = 2 : nt
%     t
%     h = squeeze(Phi(t,:,:)) * THETA;
%     y = m(S(t,:),t);
%     
%     e = 0.03;
%     cvx_begin quiet;
%         variable a(64,1);
%         minimize(norm(a,1));
%         subject to;
%             norm(h * (z(:,t-1) + a) - y) <= e;
%     cvx_end;
%     z(:, t) = z(:,t-1) + a;
% end
% 
% yp = zeros(64,nt);
% for n = 1 : nt
%     yp(:,n) = THETA * z(:,n);
% end
% 

n = 136
t = (n-1) * Dt;
t_f = t/dt + 1;
clf
EM = plot(x(p_index(S(n,:))),m(S(n,:),n)','r.','Markersize',40)
legend('real Measure Points')
xlim([0 10])
ylim([-0.5 2.5])
set(gcf,'outerposition',get(0,'screensize'));
txt = ['$t = ',num2str((n-1)*Dt),'$'];
z = text(0.8,0.6,txt,'FontSize',30);
set(z,'Interpreter','latex')
txt = ['$M = ',num2str(number),'$'];
z = text(0.8,0.8,txt,'FontSize',30);
set(z,'Interpreter','latex')
drawnow
setplt('Temperature Distribution','$x$','$f$','CS-idea 1',printfigure)

figure
n = 136
t = (n-1) * Dt;
t_f = t/dt + 1;
clf
reSignal = plot(x, f_e_scr_kf_fdm(:,n),'-','LineWidth',5)
set(reSignal,'Color',[115 121 244]/255)
hold on
EM = plot(x(p_index(S(n,:))),m(S(n,:),n)','r.','Markersize',40)
legend('estimated Signal by CS','real Measure Points')
xlim([0 10])
ylim([-0.5 2.5])
set(gcf,'outerposition',get(0,'screensize'));
txt = ['$t = ',num2str((n-1)*Dt),'$'];
z = text(0.8,0.6,txt,'FontSize',30);
set(z,'Interpreter','latex')
txt = ['$M = ',num2str(number),'$'];
z = text(0.8,0.8,txt,'FontSize',30);
set(z,'Interpreter','latex')
drawnow
setplt('Temperature Distribution','$x$','$f$','CS-idea 2',printfigure)



Phi_1024 = zeros(nt, number, 1024);
for t = 1 : nt
    for i = 1 : number
        Phi_1024(t,i,p_index(S(t,i))) = 1;
    end
end

THETA_1024 = idct(eye(1024));
MT_1024 = THETA_1024^-1;

z_1024 = zeros(1024, nt);
z_1024(:,1) = MT_1024 * f(:,1);

t = 136;
h_1024 = squeeze(Phi_1024(t,:,:)) * THETA_1024;
y = m(S(t,:),t);
e = 0.03;
cvx_begin quiet;
    variable z_1024(1024,1);
    minimize(norm(z_1024,1));
    subject to;
        norm(h_1024 * z_1024 - y) <= e;
cvx_end;

f_1024 = THETA_1024 * z_1024;

figure
n = 136
t = (n-1) * Dt;
t_f = t/dt + 1;
clf
Signal = plot(x,f(:,t_f),'-','LineWidth',5)
set(Signal,'Color',[0.9,0.9,0.9])
hold on
reSignal = plot(x, f_1024,'-','LineWidth',5)
set(reSignal,'Color',[115 121 244]/255)
hold on
EM = plot(x(p_index(S(n,:))),m(S(n,:),n)','r.','Markersize',40)
legend('Signal','estimated Signal by CS','real Measure Points')
xlim([0 10])
ylim([-0.5 2.5])
set(gcf,'outerposition',get(0,'screensize'));
txt = ['$t = ',num2str((n-1)*Dt),'$'];
z = text(0.8,0.6,txt,'FontSize',30);
set(z,'Interpreter','latex')
txt = ['$M = ',num2str(number),'$'];
z = text(0.8,0.8,txt,'FontSize',30);
set(z,'Interpreter','latex')
drawnow
setplt('Temperature Distribution','$x$','$f$','CS-idea 3',printfigure)


