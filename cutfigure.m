function Inew = cutfigure(filename)
['Processing "',filename,'" ....']
I = imread(filename);
Ig = rgb2gray(I);
[m,n] = size(Ig);
cutedge = [1 1 m n];
for i = 1 : m
    if sum(Ig(i,:)) < n * 255
        cutedge(1) = i;
        break
    end
end
for i = 1 : n
    if sum(Ig(:,i)) < m * 255
        cutedge(2) = i;
        break
    end
end
for i = m : -1 : 1
    if sum(Ig(i,:)) < n * 255
        cutedge(3) = i;
        break
    end
end
for i = n : -1 : 1
    if sum(Ig(:,i)) < m * 255
        cutedge(4) = i;
        break
    end
end
Cutedge = cutedge + [-100 -100 100 100];
Cutedge(1) = max(Cutedge(1),1);
Cutedge(2) = max(Cutedge(2),1);
Cutedge(3) = min(Cutedge(3),m);
Cutedge(4) = min(Cutedge(4),n);

Inew = I(Cutedge(1):Cutedge(3),Cutedge(2):Cutedge(4),:);

