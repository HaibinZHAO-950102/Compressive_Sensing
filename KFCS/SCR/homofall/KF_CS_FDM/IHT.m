function [ z ] = IHT( y,theta,K,mu,epsilon,loopmax )    
    if nargin < 6    
        loopmax = 3000;    
    end    
    if nargin < 5      
        epsilon = 1e-3;      
    end     
    if nargin < 4      
        mu = 1;      
    end     
    [x_rows,x_columns] = size(y);      
    if x_rows<x_columns      
        y = y';%x should be a column vector      
    end    
    n = size(theta,2);    
    z = zeros(n,1);%Initialize y=0    
    loop = 0;    
    while(norm(y-theta*z)>epsilon && loop < loopmax)    
        z = z + theta'*(y-theta*z)*mu;%update y    
        %the following two lines of code realize functionality of H_M(.)    
        %1st: permute absolute value of y in descending order    
        [ysorted inds] = sort(abs(z), 'descend');    
        %2nd: set all but M largest coordinates to zeros    
        z(inds(K+1:n)) = 0;    
        loop = loop + 1;    
    end    
end   