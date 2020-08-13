%==========================================================================
% Author: Sujit Kumar Sahoo, School of Electrical and Electronic
% Engineering, Nanyang Technological University, Singapore.
% -----------------------------------------------------------------------
% Demo 1 of the following article.

% [1] Sahoo, S.K.; Makur, A., "Signal recovery from random measurements via
% extended orthogonal matching pursuit", Transactions on Signal Processing,
% IEEE , vol.63, no.10, pp.2572-2581, May 2015.

%-------------------Set the Parameterss-------------------
d = 1024; % Signal length
P = 1000; % Number of test sample signals

fac = floor(2*d/256);
M = fac*(9:13); % sparsity

N = 64;

A = 0:0.02:1.2;

%% --- Initializing Sparse Vectors
Sr = zeros(d,P);
%--- Initializing Results
ResultOMP = zeros(length(A),length(M));

%-------------------Initialize Dictionaries---------------
Phi = randn(N,d);
Phi = Phi./repmat(sqrt(sum(Phi.^2,1)),N,1);

for k = 1:length(M) % sparsity of the signal
    tic;
    
    m = M(k);
    disp(['starting m = ',num2str(m)]);
    disp(' ');
    S = GetSparseData(m,d,P,2);
    
    %-------------------Sensing-------------------------------
    V = Phi*S;
    
    %-------------------Recovery------------------------------
    %         for i = 1:length(A)
    %         alpha = A(i);
    %         tic;
    %         Sr = OMPa(Phi,V,m,alpha);
    %         toc;
    %         ResultOMP(i,k)  = 100*mean(sum((S-Sr).^2)<10^-6);
    %         end
    
    
    Sr = OMPvar(Phi,V,m,A);
    toc;
    for i = 1:length(A)
        ResultOMP(i,k)  = 100*mean(sum((S-Sr(:,:,i)).^2)<10^-6);
    end
    
    toc;
end

%% Plotting % of Success
Line = {':','-.','-',':','-.','-',':','-.','-'};
Marker = ['o','+','*','x','s','d'];
Color = ['b','k','r','m','g','c','y','b','g','r','m','k','c','y','b','g','r','m','k','c','y','b','g','r','m','k','c','y'];


figure,hold;
xlabel('\alpha'), ylabel('% of signal recovered');
title('Performance Comparison')

Percentage = A;
for k = 1:length(M)
    Percentage(:) = ResultOMP(:,k);
    plot(A,Percentage,['-o',Color(k)]);
end
hold;

str = {};
for k = 1:length(M)
    str = [str, {['m = ', num2str(M(k))]}];
end
legend(str);
