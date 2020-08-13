%==========================================================================
% Author: Sujit Kumar Sahoo, School of Electrical and Electronic
% Engineering, Nanyang Technological University, Singapore.
% -----------------------------------------------------------------------
% Demo 2 of the following article.

% [1] Sahoo, S.K.; Makur, A., "Signal recovery from random measurements via 
% extended orthogonal matching pursuit", Transactions on Signal Processing, 
% IEEE , vol.63, no.10, pp.2572-2581, May 2015.
  
%-------------------Set the Parameterss-------------------
  d = 256; % Signal length
  P = 1000; %Number of sample signals
  
  fac = floor(4*d/256);
  M = fac*(1:8);
  
  step = floor(2*d/256);
  NoM = step:step:d;
  
  A = [0,1/16,1/8,1/4,1/2,1,inf];
  
  % --- Initializing Sparse Vectors
  Sr = zeros(d,P);
  %--- Initializing Results
  ResultOMP = zeros(length(A),length(M),length(NoM));
  
  %% ----------Recovery-------------
  % matlabpool(2)
  for k = 1:length(M) % sparsity of the signal
      tic;
      
      m = M(k);
      disp(['starting m = ',num2str(m)]);
      disp(' ');
      S = GetSparseData(m,d,P,1);
      
      Nmin = max(1.70*m*log(d/m)-60, 2*m);
      Nmax = min(1.94*m*log(d)+30, d);
      
      for n = floor(Nmin/step):ceil(Nmax/step)
          N = NoM(n);
          %-------------------Initialize Dictionaries---------------
          Phi = randn(N,d);
          Phi = Phi./repmat(sqrt(sum(Phi.^2,1)),N,1);

          %-------------------Sensing-------------------------------
          V = Phi*S;
          
          %-------------------Recovery------------------------------
          %         for i = 1:length(A)
          %         alpha = A(i);
          %         Sr = OMPa(Phi,V,m,alpha);
          %         ResultOMP(i,k,n)  = 100*mean(sum((S-Sr).^2)<10^-6);
          %         end
          
          Sr = OMPvar(Phi,V,m,A);
          for i = 1:length(A)
              ResultOMP(i,k,n)  = 100*mean(sum((S-Sr(:,:,i)).^2)<10^-6);
          end
      end
      ResultOMP(:,k,(n+1):end) = 100;
      
      toc;
  end
  % matlabpool close;
  
%% Plotting 
  %-- Plot Specs---
  Line = {':','-.','-',':','-.','-',':','-.','-'};
  Marker = ['o','*','+','x','s','d','o','*','+','x','s','d'];
  Color = ['r','g','b','r','k','m','c','y','b','g','r','m','k','c','y','b','g','r','m','k','c','y','b','g','r','m','k','c','y'];
  
%--- Probabilitypercentage of Success--
  figure,hold;  
  c=0;
  for k = 2:4:length(M)
      c = c+1;   s=0;
      for i = 1:length(A)
          s = s+1;
          Percentage(:) = ResultOMP(i,k,:);
          plot(NoM,Percentage,['-',Marker(c),Color(s)]);
      end
  end
hold;
  
str = {};
for k = 2:4:length(M)
    for i = 1:length(A)
    str = [str, {['m = ', num2str(M(k)),' \alpha = ', num2str(A(i))]}];
    end
end
legend(str);

axis([NoM(1) 1.1*NoM(end) 0 100]);
xlabel('No. of measurements (N)'), ylabel('% of signal recovered');
title('Performance Comparison')

%% --- N vs m for 95% Success probability ---- 
Max = 7;
PC = 95; 
OverSamplingOMP = zeros(length(A),length(M));

figure,hold;
for i = 1:length(A)
   for k = 1:length(M)
       % sparsity of the signal
       for n = 1:length(NoM)
           if ResultOMP(i,k,n) >= PC;
               OverSamplingOMP(i,k)=NoM(n);
               break;
           end
       end
   end
   Max = sum(OverSamplingOMP(i,:)>0);
   plot(M(1:Max),OverSamplingOMP(i,1:Max),[':',Color(i),Marker(i)])
end
hold;

hold;
for i = 1:length(A)
   a = A(i);
   if a == inf
       f = fittype({['x*log(',num2str(d),'/x)'],'1'})
   else
       f = fittype({['x*log(',num2str(d),'/(x*',num2str(a),'+1))'],'1'})
   end
   Max = sum(OverSamplingOMP(i,:)>0);
   [Curve,gof2] = fit(M(1:Max)',OverSamplingOMP(i,1:Max)',f)
   plot(Curve,Color(i));
end
hold;

str = {};
for i = 1:length(A)
    str = [str, {['\alpha = ', num2str(A(i))]}];
end
legend(str);

axis([M(1) 1.1*M(end) 1 d ]);
xlabel('Sparsity (m)'), ylabel('No of measurement (N)');
title(['Oversampling Factor needed for ', num2str(PC),'% success'])

