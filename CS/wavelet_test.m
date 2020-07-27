clc
clear
close all

printfigure = 0;

N = 8192;
T = 20;
t = 0 : T/(N-1) : T;
x = sin(2*pi*t);

plot(t,x,'k-','LineWidth',5)
setplt('Function','$t$','$x$','function',printfigure)

Phi_wt = DWT(N,'haar');
Phi_iwt = Phi_wt^-1;

y_wt = Phi_wt * x';
[y2,y2l] = wavedec(x,log2(N),'haar');

figure
plot(y_wt,'r.')
xlim([0 8192])
setplt('Wavelet Matrix','','','Wavelet Matrix',printfigure)
figure
plot(y2,'b.')
xlim([0 8192])
setplt('Wavelet Transformation','','','Wavelet Transformation',printfigure)
