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
            [f_sr,f_mu,k,Length,dt,dx,p,m_rh,p_index] = maketemperature(sigma_mu, sigma_sr);
            
            processing = ['Kalman Filter']
            f_e = kalman(sigma_mu, sigma_sr, f_sr,f_mu,k,Length,dt,dx,p,m_rh,p_index);
            
            processing = ['Kalman Filtered Compressive Sensing']
            f_e_kfcs = kfcs(sigma_mu, sigma_sr, f_sr,f_mu,k,Length,dt,dx,p,m_rh,p_index);

            
            e_kf = bewertung(f_e, f_sr);
            e_kfcs = bewertung(f_e_kfcs, f_sr);
            
            read_kf = xlsread('error_kf.xlsx');
            read_kfcs = xlsread('error_kfcs.xlsx');
            
            write_kf = [read_kf;kkk, iii, jjj, e_kf];
            write_kfcs = [read_kfcs;kkk, iii, jjj, e_kfcs];
            
            xlswrite('error_kf.xlsx', write_kf);
            xlswrite('error_kfcs.xlsx', write_kfcs);
            
        end
    end
end