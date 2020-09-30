clc
clear
close all

Sigma_sr = 0.001: 0.001 : 0.01;% Systemrauschen
Sigma_mu = 0.01: 0.01 : 0.1; % Messunsicherheit


for kkk = 1 : 1000
    for iii = 1 : 10
        for jjj = 1 : 10
            clc
            
            sigma_mu = Sigma_mu(iii);
            sigma_sr = Sigma_sr(jjj);
            name = [num2str(kkk),'_',num2str(iii),'_',num2str(jjj)]
            
            processing = ['Generating Temperature Distribution']
            maketemperature(sigma_mu, sigma_sr, name);
            
            processing = ['Kalman Filter']
            kalman(sigma_mu, sigma_sr, name);
            
            processing = ['Kalman Filtered Compressive Sensing']
            kfcs(sigma_mu, sigma_sr, name);

            LOAD = ['load ',name];
            eval(LOAD);
            LOAD = ['load f_e_',name];
            eval(LOAD);
            LOAD = ['load f_e_kfcs_',name];
            eval(LOAD);
            
            e_kf = bewertung(f_e, f_sr);
            e_kfcs = bewertung(f_e_kfcs, f_sr);
            
            read_kf = xlsread('error_kf.xlsx');
            read_kfcs = xlsread('error_kfcs.xlsx');
            
            write_kf = [read_kf;kkk, iii, jjj, e_kf];
            write_kfcs = [read_kfcs;kkk, iii, jjj, e_kfcs];
            
            xlswrite('error_kf.xlsx', write_kf);
            xlswrite('error_kfcs.xlsx', write_kfcs);
            
            DELET = ['delete ',name,'.mat'];
            eval(DELET);
            DELET = ['delete f_e_',name,'.mat'];
            eval(DELET);
            DELET = ['delete f_e_kfcs_',name,'.mat'];
            eval(DELET);
        end
    end
end