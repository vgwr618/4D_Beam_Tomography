function [init4DCoord, fin4DCoord, finIm, finImDisp, diffMat, posCoord, negCoord, sigmaRecon, stReconMeasXY] = first_propagation(InitScreenIm, MeasuredIm, R, npar, sigma, pixcal, resfactor)
if sum(InitScreenIm, 'all') == npar
    initXCoord = [];
    initYCoord = [];
    pixelDim = size(InitScreenIm);
    for i = 1 : pixelDim(1)
        for j = 1 : pixelDim(2)
            if InitScreenIm(i,j) > 0
                for k = 1 : InitScreenIm(i,j)
                    initXCoord(end+1) = pixcal * rand + (resfactor * pixcal) * (j - (pixelDim(2) / 2) / resfactor) - pixcal;
                    initYCoord(end+1) = pixcal * rand + (resfactor * pixcal) * (-i + (pixelDim(1) / 2) / resfactor) - pixcal;
                end
            end
        end
    end
    
    initPxCoord = normrnd(0, sigma, [1, length(initXCoord)]);
    initPyCoord = normrnd(0, sigma, [1, length(initXCoord)]);
    init4DCoord = [initXCoord;
                   initPxCoord;
                   initYCoord;
                   initPyCoord];
    fin4DCoord = R * init4DCoord;
    finXCoord = fin4DCoord(1,:);
    finPxCoord = fin4DCoord(2,:);
    finYCoord = fin4DCoord(3,:);
    finPyCoord = fin4DCoord(4,:);
    finIm = zeros(floor(pixelDim(1)/resfactor), floor(pixelDim(2)/resfactor));
    for i = 1:length(finXCoord)
        xpos = floor(finXCoord(i) / (resfactor * pixcal) + (pixelDim(2) / 2) / resfactor);
        ypos = floor(-finYCoord(i) / (resfactor * pixcal) + (pixelDim(1) / 2) / resfactor);
        if xpos > 0 && xpos <= 800 / resfactor && ypos > 0 && ypos <= 800 / resfactor
            finIm(ypos,xpos) = finIm(ypos,xpos) + 1;
        end
    end
    radius = 1; %pixels
    filter = fspecial('disk', radius);
    finImDisp = conv2(finIm,filter);
    finImDisp = imresize(finImDisp,[800 800]);
    finImDisp = uint8(finImDisp);
    sum(finImDisp, 'all')
    sum(finIm, 'all')
    finImDisp = (1/sum(finImDisp, 'all')) * double(finImDisp);
    diffMat =  finImDisp - double(MeasuredIm);
    posX = [];
    posY = [];
    negX = [];
    negY = [];
    for i = 1:pixelDim(1)
        for j = 1:pixelDim(2)
            if diffMat(i,j) < 0
                for k = 1:abs(diffMat(i,j) * npar)
                    negX(end+1) = (resfactor * pixcal) * (j - (pixelDim(2) / 2) / resfactor);
                    negY(end+1) = (resfactor * pixcal) * (-i + (pixelDim(1) / 2) / resfactor);
                end
            elseif diffMat(i,j) > 0
                for k = 1:abs(diffMat(i,j) * npar)
                    posX(end+1) = (resfactor * pixcal) * (j - (pixelDim(2) / 2) / resfactor);
                    posY(end+1) = (resfactor * pixcal) * (-i + (pixelDim(1) / 2) / resfactor);
                end
            end
        end
    end
    posCoord = [posX; 
                posY];
    negCoord = [negX; 
                negY];
            
    measXCoord = [];
    measYCoord = [];
    for i = 1 : pixelDim(1)
        for j = 1 : pixelDim(2)
            if MeasuredIm(i,j) > 0
                for k = 1 : MeasuredIm(i,j)
                    measXCoord(end+1) = (resfactor * pixcal) * (j - (pixelDim(2) / 2) / resfactor);
                    measYCoord(end+1) = (resfactor * pixcal) * (-i + (pixelDim(1) / 2) / resfactor);
                end
            end
        end
    end
    meas4DCoord = zeros(length(measXCoord),4);
    meas4DCoord(:,1) = measXCoord;
    meas4DCoord(:,2) = finPxCoord(1:length(meas4DCoord(:,2)))';
    meas4DCoord(:,3) = measYCoord;
    meas4DCoord(:,4) = finPyCoord(1:length(meas4DCoord(:,4)))';
    sigmaRecon = cov(fin4DCoord');
    stdxRecon = sqrt(sigmaRecon(1,1));
    stdyRecon = sqrt(sigmaRecon(3,3));
    sigmaMeas = cov(meas4DCoord);
    stdxMeas = sqrt(sigmaMeas(1,1));
    stdyMeas = sqrt(sigmaMeas(3,3));
    stReconMeasXY = [stdxRecon, stdyRecon, stdxMeas, stdyMeas];
else
    if sum(InitScreenIm, 'all') ~= npar
        error('Error. \nInitial screen not normalized to match the number of particles.')
    elseif sum(MeasuredIm, 'all') ~= npar
        error('Error. \nMeasured screen not normalized to match the number of particles.')
    end
end


