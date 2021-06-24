function [init4DCoord] = backward_propagation(init4DCoord, fin4DCoord, diffMat, pixcal, resfactor)
pixelDim = size(diffMat);
for i = 1 : pixelDim(1)
    for j = 1 : pixelDim(2)
        if diffMat(i,j) > 0
            index_list = [];
            pixLeftX = pixcal * ((resfactor * j) - (pixelDim(2) / 2));
            pixTopY = pixcal * ((resfactor * -i) + (pixelDim(1) / 2));
            for k = 1:length(fin4DCoord(1,:))
                if fin4DCoord(1,k) <= (pixLeftX + resfactor * pixcal) && fin4DCoord(1,k) >= pixLeftX && fin4DCoord(3,k) <= pixTopY && fin4DCoord(3,k)>= pixTopY - resfactor * pixcal
                    index_list(end+1) = k;
                end
            end
            coordInBox = length(index_list);
            numDiffMat = floor(diffMat(i,j));
            coordInBox-numDiffMat;
%             i
%             j
            if ~isempty(index_list) && coordInBox > numDiffMat
                for m = 1:floor(diffMat(i,j))
                    n = index_list(m);
%                     fin4DCoord(1,n)
%                     fin4DCoord(3,n)
%                     pixLeftX
%                     pixTopY
                    init4DCoord(2, n) = 1e-4*rand-5e-5;
                    init4DCoord(4, n) = 1e-4*rand-5e-5;
                end
            end
        end
    end
end