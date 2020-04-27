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
setplt('Signal_Spectrum','$f$','$|A(f)|$','Signal_Spectrum',printfigure)

sampling_frequence = 50;
sampling_number = 3 * sampling_frequence;
sample_index_shannon = ceil(linspace(1, N + 1, sampling_number));
t_shannon = t(sample_index_shannon);
y_shannon = y(sample_index_shannon);

figure
plot(t,y,'k','LineWidth',2)
hold on
plot(t_shannon,y_shannon,'b-','LineWidth',2)
hold on
plot(t_shannon,y_shannon,'r.','markersize',30)
setplt('Traditional Sampling','$t$','$y(t)$','Signal_traditional_sampling',printfigure)


CS_sampling_number = ceil(sampling_number * 0.2);
CS_sample_index = sort(ceil(rand(1, CS_sampling_number) * N));
t_sample = t(CS_sample_index);
y_sample = y(CS_sample_index);

figure
plot(t,y,'k','LineWidth',2)
hold on
plot(t_sample,y_sample,'b-','LineWidth',2)
hold on
plot(t_sample,y_sample,'r.','markersize',30)
setplt('CS Sampling','$t$','$y(t)$','Signal_CS_sampling',printfigure)

Fourier_matrix = zeros(CS_sampling_number, 31);
for i = 1 : 31
    Fourier_matrix(:,i) = sin(2*pi*(i-1)*t_sample');
end

e = 0.0001;
fun = @(a) sum(abs(a));
con = @(a)constraint_1(Fourier_matrix,a,y_sample',e);
a = fmincon(fun,ones(31,1),[],[],[],[],[],[],con);
threshold = floor(log(1 / (sum(abs(a)) * 0.01)) / log(10));
a = round(a,threshold);

y_re = zeros(1,length(t));
for i = 1 : 31
    y_re = y_re + a(i) * sin(2 * pi * (i-1) * t);
end

figure
plot(t(1:10001),y(1:10001),'k-','LineWidth',2)
hold on
plot(t(1:10001),y_re(1:10001),'r--','LineWidth',2)
legend('Original Signal','Reconstructed Signal')
setplt('Signal Reconsruction','$t$','$y(t)$','Signal_reconsruction',printfigure)

