clc
clear

read_kf = textread('error_kf.txt');
read_kfcs = textread('error_kfcs.txt');

e_kf = read_kf(:,1:4);
e_kfcs = read_kfcs(:,1:4);

E_KF = zeros(10,10);
E_KFCS = zeros(10,10);

for i = 1 : 10
  for j = 1 : 10
    
    [a, b] = find(e_kf(:,2) == i);
    temp = e_kf(a,:);
    [a, b] = find(temp(:,3) == j);
    temp = temp(a,:);
    E_KF(i,j) = mean(temp(:,4));

    [a, b] = find(e_kfcs(:,2) == i);
    temp = e_kfcs(a,:);
    [a, b] = find(temp(:,3) == j);
    temp = temp(a,:);
    E_KFCS(i,j) = mean(temp(:,4));
  end
end

surf(E_KF)
xticks(2:2:10)
yticks(2:2:10)
yticklabels({'0.002','0.004','0.006','0.008','0.01'})
xticklabels({'0.02','0.04','0.06','0.08','0.1'})
setmesh('Fehler','$\sigma_M$','$\sigma_S$','Fehler','error_kf',1)


figure

surf(E_KFCS)
xticks(2:2:10)
yticks(2:2:10)
yticklabels({'0.002','0.004','0.006','0.008','0.01'})
xticklabels({'0.02','0.04','0.06','0.08','0.1'})
setmesh('Fehler','$\sigma_M$','$\sigma_S$','Fehler','error_kfcs',1)


E = E_KF - E_KFCS;

figure

surf(E)
xticks(2:2:10)
yticks(2:2:10)
yticklabels({'0.002','0.004','0.006','0.008','0.01'})
xticklabels({'0.02','0.04','0.06','0.08','0.1'})
setmesh('Differenzfehler','$\sigma_M$','$\sigma_S$','Fehler','error_difference',1)



figure
plot(E(1,2:end),'linewidth',5)
hold on
plot(E(3,2:end),'linewidth',5)
hold on
plot(E(5,2:end),'linewidth',5)
hold on
plot(E(7,2:end),'linewidth',5)
hold on
plot(zeros(1,9),'k--','linewidth',3)
hold on

xticks(2:2:10)
xticklabels({'0.01','0.03','0.05','0.07','0.09'})
legend('\sigma_S = 0.001','\sigma_S = 0.003','\sigma_S = 0.005','\sigma_S = 0.007')
setplt('Differenzfehler','$\sigma_M$','Fehler','error_line',1)


