clc
clear
close all

Error_KF = zeros(4,10,10);

Sigma_sr = 0: 0.001 : 0.01% Systemrauschen
Sigma_mu = 0: 0.01 : 0.1 % Messunsicherheit

for kkk = 1 : 4
    for iii = 2 : 11
        for jjj = 2 : 11
            [kkk,iii,jjj]
            LOAD = ['load f_e_kfcs_',num2str(kkk),'_',num2str(iii),'_',num2str(jjj)];
            eval(LOAD)
            LOAD = ['load ',num2str(kkk),'_',num2str(iii),'_',num2str(jjj)];
            eval(LOAD)
            f = f_sr(:,1:10:end);
            
            Error_KF(kkk,iii-1,jjj-1) = bewertung(f_e_kfcs,f);
        end
    end
end

Meanerror_KFCS = zeros(10,10);

for iii = 1 : 10
    for jjj = 1 : 10
        Meanerror_KFCS(iii,jjj) = mean(Error_KF(:,iii,jjj));
    end
end

surf(Meanerror_KFCS)

xlswrite('Meanerror_KFCS.xlsx',Meanerror_KFCS);

            