clc
clear
close all

printfigure = 1;

dt = 0.001;
t = 0 : dt : 3;
N = length(t) - 1;
f1 = 10;
f2 = 25;
fm = max(f1,f2);
y = 2 * sin(2 * pi * f1 * t) + sin(2 * pi * f2 * t);

plot(t,y,'k','LineWidth',2)
setplt('','$t$','$y(t)$','Signal',printfigure)

figure
plot(t(1001:1501),y(1001:1501),'k','LineWidth',2)
setplt('','$t$','$y(t)$','Signal_scaled',printfigure)

A = fft(y);
df = 1 / dt / N;
f = 0 : df : 1 / dt;
f = f - ceil(max(f)/2);
A = [A(ceil((N+1)/2):N+1),A(1:ceil((N+1)/2)-1)];
figure
plot(f,abs(A),'k','LineWidth',3)
xlim([-1.5*fm 1.5*fm])
setplt('','$f$','$|A(f)|$','Signal_Spectrum',printfigure)

sampling_frequence = 2*fm;
sampling_number = 3 * sampling_frequence;
sample_index_shannon = ceil(linspace(1, N + 1, sampling_number));
t_shannon = t(sample_index_shannon);
y_shannon = y(sample_index_shannon);

figure
plot(t,y,'k','LineWidth',3)
hold on
L = plot(t_shannon,y_shannon,'--','LineWidth',1)
set(L,'Color',[0.5 0.5 0.5])
hold on
plot(t_shannon,y_shannon,'r.','markersize',40)
setplt('','$t$','$y(t)$','Signal_traditional_sampling',printfigure)


CS_sampling_number = 30;
CS_sample_index = sort(ceil(rand(1, CS_sampling_number) * N));
t_sample = t(CS_sample_index);
y_sample = y(CS_sample_index);

figure
plot(t,y,'k','LineWidth',3)
hold on
L2 = plot(t_sample,y_sample,'--','LineWidth',1)
set(L2,'Color',[0.5 0.5 0.5])
hold on
plot(t_sample,y_sample,'r.','markersize',40)
setplt('','$t$','$y(t)$','Signal_CS_sampling',printfigure)

Psi = zeros(N + 1, N + 1);
for i = 1 : N + 1
    for j = 1 : N + 1
        Psi(i,j) = sin(2*pi*(j-1)* t(i));
    end
end
Phi = zeros(CS_sampling_number, N + 1);
for i = 1 : CS_sampling_number
    Phi(i,CS_sample_index(i)) = 1;
end

A = Phi * Psi;

e = 0.00001;
cvx_begin
    variable a(N + 1,1)
    minimize(norm(a,1))
    subject to
        norm(A * a - y_sample') <= e
cvx_end

y_re = Psi * a;

figure
plot(t(1:1000),y(1:1000),'k-','LineWidth',5)
hold on
plot(t(1:1000),y_re(1:1000),'r--','LineWidth',5)
legend('Originales Signal','Rekonstruiertes Signal')
setplt('','$t$','$y(t)$','Signal_reconsruction',printfigure)
close all
