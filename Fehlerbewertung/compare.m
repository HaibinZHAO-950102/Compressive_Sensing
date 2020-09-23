clc
clear
close all

E_KF = xlsread('Meanerror_KF.xlsx');
E_KFCS = xlsread('Meanerror_KFCS.xlsx');

surf(E_KF)

hold on
surf(E_KFCS)

% figure
% 
% dE = E_KF - E_KFCS;
% surf(dE)