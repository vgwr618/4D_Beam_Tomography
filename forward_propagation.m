function [init4DCoord, fin4DCoord, finIm, finImDisp, diffMat] = forward_propagation(init4DCoord, MeasuredIm, R, pixcal, resfactor)
pixelDim = size(MeasuredIm);
fin4DCoord = R * init4DCoord;
finXCoord = fin4DCoord(1,:);
finYCoord = fin4DCoord(3,:);
finIm = zeros(floor(pixelDim(1)/resfactor), floor(pixelDim(2)/resfactor));
for i = 1:length(finXCoord)
    xpos = floor(finXCoord(i) / (resfactor * pixcal) + (pixelDim(2) / 2) / resfactor);
    ypos = floor(-finYCoord(i) / (resfactor * pixcal) + (pixelDim(1) / 2) / resfactor);
    if xpos >= 0 && xpos < 800 / resfactor && ypos >= 0 && ypos < 800 / resfactor
        finIm(ypos,xpos) = finIm(ypos,xpos) + 1;
    end
end
newpart =  sum(finIm,'all');
radius = 1; %pixels
filter = fspecial('disk', radius);
finImDisp = conv2(finIm,filter);
finImDisp = imresize(finImDisp,[800 800]);
finImDisp = uint8(finImDisp);
% sum(finImDisp, 'all')
% sum(finIm, 'all')
finImDisp = (newpart/sum(finImDisp, 'all')) * double(finImDisp);
MeasuredIm = (newpart/sum(MeasuredIm, 'all')) * double(MeasuredIm);
diffMat =  finIm - MeasuredIm;