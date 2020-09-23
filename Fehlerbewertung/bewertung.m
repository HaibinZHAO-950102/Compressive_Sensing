function e = bewertung(f1,f)

error_1 = f - f1;

error1 = zeros(size(f,2),1);

for i = 1 : size(f,2)
    error1(i) = norm(error_1(:,i)) / norm(f(:,i));
end

e = mean(error1);
