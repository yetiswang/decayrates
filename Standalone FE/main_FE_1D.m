function main_FE_1D
% 
% height = [  45 50 55 58 60 62 65 70 73 75 80 83];
% diameter = [ 30 30 30 30 30 30 30 30 30 30 30 30];

height = 75 ;
diameter = 20 ;

metal = {'Au'};
parent = pwd; 
for i = 1 : length(metal)
    for j = 1 : length(height)
        if isequal(height(j),diameter(j))
            FE_GNS_1D(diameter(j),metal{i});
            close all
        else
            FE{j} = FE_GNR_1D(height(j),diameter(j),metal{i},637, 647);
            %FE{j} = FE_GNR_1D(height(j),diameter(j),metal{i});
            close all
        end
        cd(parent)
    end
    cd(parent)
end


end