function Inew = cutfigure(filename)
['Processing "',filename,'" ....']

[I, map]=imread(filename,'frame','all') ;
[rows, columns, numColorChannels, numImages] = size(I);

Frame = I(:,:,:, 1);
I_k = uint8(255 * ind2rgb(Frame, map));

Ig = rgb2gray(I_k);
[m,n] = size(Ig);
cutedge = [1 1 m n];

criterion = sum(Ig(:,1)) / m;

for i = 1 : m
    if sum(Ig(i,:)) ~= n * criterion
        cutedge(1) = i;
        break
    end
end
for i = 1 : n
    if sum(Ig(:,i)) ~= m * criterion
        cutedge(2) = i;
        break
    end
end
for i = m : -1 : 1
    if sum(Ig(i,:)) ~= n * criterion
        cutedge(3) = i;
        break
    end
end
for i = n : -1 : 1
    if sum(Ig(:,i)) ~= m * criterion
        cutedge(4) = i;
        break
    end
end

Cutedge = cutedge + [-20 -20 20 20];
Cutedge(1) = max(Cutedge(1),1);
Cutedge(2) = max(Cutedge(2),1);
Cutedge(3) = min(Cutedge(3),m);
Cutedge(4) = min(Cutedge(4),n);


for k = 1 : numImages
    Frame = I(:,:,:, k);
    I_cut = Frame(Cutedge(1):Cutedge(3),Cutedge(2):Cutedge(4));
   
    if k == 1
         imwrite(I_cut,map,filename,'gif', 'Loopcount',inf,'DelayTime',1e-6);
    else
         imwrite(I_cut,map,filename,'gif','WriteMode','append','DelayTime',1e-6);
    end
end


