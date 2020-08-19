function ww = DWT(N, wavelet)

% 'db1' or 'haar', 'db2', ... ,'db10', ... , 'db45'
% 'coif1', ... , 'coif5'
% 'sym2', ... , 'sym8', ... ,'sym45'
% 'bior1.1', 'bior1.3', 'bior1.5'
% 'bior2.2', 'bior2.4', 'bior2.6', 'bior2.8'
% 'bior3.1', 'bior3.3', 'bior3.5', 'bior3.7'
% 'bior3.9', 'bior4.4', 'bior5.5', 'bior6.8'
% 'rbio1.1', 'rbio1.3', 'rbio1.5'
% 'rbio2.2', 'rbio2.4', 'rbio2.6', 'rbio2.8'
% 'rbio3.1', 'rbio3.3', 'rbio3.5', 'rbio3.7'
% 'rbio3.9', 'rbio4.4', 'rbio5.5', 'rbio6.8'

[hp,tp] = wfilters(wavelet,'d');

L = length(hp);
Z_max = log2(N);
Z_min = floor(log2(L))+1;
ww=1;

for schichte = Z_min : Z_max
    size = 2 ^ schichte;

    p_hp = sparse([hp,zeros(1,size-L)]);
    p_tp = sparse([tp,zeros(1,size-L)]);

    for i=1:size/2
        p1(i,:) = circshift(p_hp',2*(i-1))';
        p2(i,:) = circshift(p_tp',2*(i-1))';
    end
    
    w1 = [p1;p2];
    mm = N - length(w1);
    w = sparse([w1,zeros(length(w1),mm);zeros(mm,length(w1)),eye(mm,mm)]);
    ww = ww * w;
    clear p1;clear p2;
end
