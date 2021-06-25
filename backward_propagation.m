function [init4DCoord] = backward_propagation(init4DCoord, fin4DCoord, diffMat, pixcal, resfactor, sigma)
pixelDim = size(diffMat);
pos_count = 0;
neg_count = 0;
for i = 1 : pixelDim(1)
    for j = 1 : pixelDim(2)
        if diffMat(i,j) > 0
            index_list = [];
            pixRightX = pixcal * ((resfactor * j) - (pixelDim(2) / 2));
            pixBotY = pixcal * ((resfactor * -i) + (pixelDim(1) / 2));
            for k = 1:length(fin4DCoord(1,:))
                if fin4DCoord(1,k) < pixRightX && fin4DCoord(1,k) >= pixRightX - resfactor * pixcal && fin4DCoord(3,k) > pixBotY && fin4DCoord(3,k)<= pixBotY + resfactor * pixcal
                    index_list(end+1) = k;
                end
            end
            coordInBox = length(index_list);
            numDiffMat = floor(diffMat(i,j));
            if coordInBox < numDiffMat
                neg_count = neg_count + 1;
            elseif coordInBox > numDiffMat
                pos_count = pos_count + 1;
            end
%             i
%             j
            if ~isempty(index_list) && coordInBox > numDiffMat
                for m = 1:floor(diffMat(i,j))
                    n = index_list(m);
%                     fin4DCoord(1,n)
%                     fin4DCoord(3,n)
%                     pixLeftX
%                     pixTopY
                    init4DCoord(2, n) = normrnd(0, sigma);
                    init4DCoord(4, n) = normrnd(0, sigma);
                end
            end
        end
    end
end
neg_count
pos_count