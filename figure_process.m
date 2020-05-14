clc
clear
close all

cut_figure = 1;
cut_gif = 1;

files = dir('E:\Deutschland\Karlsruher Institut f¨¹r Technologie\4. Semester\Compressive Sensing + Filterung\Codes Compressive');

if cut_figure == 1
    for i = 1 : size(files,1)
        NAME = files(i).name;
        l = size(NAME,2);
        if l >= 4
            TYPE = NAME(l-2:l);
            if TYPE == 'png'
                I = cutfigure(NAME);
                imwrite(I,NAME)
            end
        end
    end
end
    
if cut_gif == 1
    for i = 1 : size(files,1)
        NAME = files(i).name;
        l = size(NAME,2);
        if l >= 4
            TYPE = NAME(l-2:l);
            if TYPE == 'gif'
                cutgif(NAME);
            end
        end
    end
end
