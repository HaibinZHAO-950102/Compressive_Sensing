function H = coherence(A)
n = size(A,2);
h = zeros(n*(n-1)/2,1);
k = 1;
for i = 1 : n-1
    for j = i+1 : n
        p1 = A(:,i);
        p2 = A(:,j);
        h(k) = abs(p1'*p2)/norm(p1)/norm(p2);
        k = k + 1;
    end
end
H = max(h);