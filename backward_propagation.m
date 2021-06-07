function [init4DCoord] = backward_propagation(init4DCoord, fin4DCoord, diffMat, pixcal, resfactor)
pixelDim = size(diffMat);
for i = 1 : pixelDim(1)
    for j = 1 : pixelDim(2)
        if diffMat(i,j) > 0
            for k = 1:length(fin4DCoord(1,:))
                pixCenterX = (resfactor * pixcal) * (j - (pixelDim(2) / 2) / resfactor);
                pixCenterY = (resfactor * pixcal) * (-i + (pixelDim(1) / 2) / resfactor);
                if fin4DCoord(1,k) < (pixCenterX + (resfactor * pixcal) / 2) && fin4DCoord(1,k) >= (pixCenterX - (resfactor * pixcal) / 2) && fin4DCoord(3,k)< pixCenterY + (resfactor * pixcal) / 2 && fin4DCoord(3,k)>= pixCenterY - (resfactor * pixcal) / 2
                    init4DCoord(2, k) = 1e-4*rand;
                    init4DCoord(4, k) = 1e-4*rand;
                end 
            end
        end
    end
end