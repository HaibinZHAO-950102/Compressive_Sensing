function X = GetSparseData(m,d,Ns,type)
%==========================================================================
% Author: Sujit Kumar Sahoo, School of Electrical and Electronic
% Engineering, Nanyang Technological University, Singapore.
% -----------------------------------------------------------------------
%-------------------Generate Sparse Signals----------------

x = zeros(d,1);
X = zeros(d,Ns);

if type ==1
    for i = 1:Ns
        idx = randperm(d);
        %     x(idx(1:m)) = 2*randn(m,1);
        x(idx(1:m)) = 2;
        
        %     idx = (binornd(1,m/d,d,1))>0;
        %     x(idx) = 5;
        %
        %     id = randi(d);
        %     if id > d-m
        %         id = id -m;
        %     end
        %         x(id+(1:m)) = 10*randn(m,1);
        %     x(id+(1:m)) = 5;
        
        X(:,i) = x;
        x(:) = 0;
    end
end

if type == 2
    for i = 1:Ns
        idx = randperm(d);
        x(idx(1:m)) = 2*randn(m,1);
        X(:,i) = x;
        x(:) = 0;
    end
end

if type ==3
    for i = 1:Ns
        id = randi(d-m+1);
        
        x(id-1+(1:m)) = 2*randn(m,1);
        
        X(:,i) = x;
        x(:) = 0;
    end
end

if type ==4
    for i = 1:Ns
        id = randi(d-m+1);
        
        x(id-1+(1:m)) = 2;
        
        X(:,i) = x;
        x(:) = 0;
    end
end

return;