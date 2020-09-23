clc
clear
close all

Error_KF = zeros(10,10,10);

Sigma_sr = 0: 0.001 : 0.01% Systemrauschen
Sigma_mu = 0: 0.01 : 0.1 % Messunsicherheit

for kkk = 1 : 10
    for iii = 2 : 11
        for jjj = 2 : 11
            [kkk,iii,jjj]
            LOAD = ['load f_e_',num2str(kkk),'_',num2str(iii),'_',num2str(jjj)];
            eval(LOAD)
            LOAD = ['load ',num2str(kkk),'_',num2str(iii),'_',num2str(jjj)];
            eval(LOAD)
            f = f_sr(:,1:10:end);
            
            Error_KF(kkk,iii-1,jjj-1) = bewertung(f_e,f);
        end
    end
end

Meanerror_KF = zeros(10,10);

for iii = 1 : 10
    for jjj = 1 : 10
        Meanerror_KF(iii,jjj) = mean(Error_KF(:,iii,jjj));
    end
end

surf(Meanerror_KF)

xlswrite('Meanerror_KF.xlsx',Meanerror_KF);

            