clc
clear
close all

printfigure = 0;

dt = 0.0001;
t = 0 : dt : 3;
N = length(t) - 1;
y = 2 * sin(2 * pi * 10 * t) + sin(2 * pi * 25 * t);

plot(t,y,'k','LineWidth',2)
setplt('Signal','$t$','$y(t)$','Signal',printfigure)

figure
plot(t(1:10001),y(1:10001),'k','LineWidth',2)
setplt('Signal scaled','$t$','$y(t)$','Signal_scaled',printfigure)

A = fft(y);
df = 1 / dt / N;
f = 0 : df : 1 / dt;
f = f - ceil(max(f)/2);
A = [A(ceil((N+1)/2):N+1),A(1:ceil((N+1)/2)-1)];
figure
plot(f,abs(A),'k','LineWidth',2)
xlim([-50 50])
setplt('Signal_Frequence','$f$','$|A(f)|$','Signal_Frequence',printfigure)

sampling_frequence = 50;
sampling_number = 3 * 50;
CS_sampling_number = ceil(sampling_number * 0.1);
sample_index = sort(ceil(rand(1, CS_sampling_number) * N));
t_sample = t(sample_index);
y_sample = y(sample_index);

plot(t,y,'k','LineWidth',2)
hold on
plot(t_sample,y_sample,'b-','LineWidth',2)
hold on
plot(t_sample,y_sample,'r.','markersize',30)
setplt('Sampling','$t$','$y(t)$','Signal_sampling',printfigure)

Fourier_matrix = zeros(CS_sampling_number, 31);
for i = 1 : 31
    Fourier_matrix(:,i) = sin(2*pi*(i-1)*t_sample');
end
e = 0.1;

fun = @(a) sum(abs(a));
c = @(G,a,e) [norm(G * a) - e,0];
con = @(a)c(Fourier_matrix,a,e);
a = fmincon(fun,zeros(31,1),[],[],[],[],[],[],con);

    
