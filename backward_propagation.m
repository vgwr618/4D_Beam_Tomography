function [init4DCoord] = backward_propagation(init4DCoord, fin4DCoord, diffMat, pixcal, resfactor, sigma)
pixelDim = size(diffMat);
%debugging code
pos_count = 0;
neg_count = 0;

%convert pixel (matrix) coordinates to gpt coordinates
for i = 1 : pixelDim(1)
    for j = 1 : pixelDim(2)
        if diffMat(i,j) > 0
            index_list = [];
            %right and bottom boundary lines for each pixel
            pixRightX = pixcal * ((resfactor * j) - (pixelDim(2) / 2));
            pixBotY = pixcal * ((resfactor * -i) + (pixelDim(1) / 2));
            %check for final simulated particles that are in each pixel box
            for k = 1:length(fin4DCoord(1,:))
                if fin4DCoord(1,k) < pixRightX && fin4DCoord(1,k) >= pixRightX - resfactor * pixcal && fin4DCoord(3,k) > pixBotY && fin4DCoord(3,k)<= pixBotY + resfactor * pixcal
                    index_list(end+1) = k;
                end
            end
            %how many coordinates are in this pixel
            coordInBox = length(index_list);
            %what is the difference between simulated and measured image
            %for this pixel
            numDiffMat = floor(diffMat(i,j));
            
            %debugging code
            if coordInBox < numDiffMat
                neg_count = neg_count + 1;
            elseif coordInBox > numDiffMat
                pos_count = pos_count + 1;
            end
%             i
%             j
            %if there are simulated particles that landed in this box and
            %the number is greater than the difference (which it should be
            %but right now that's not always true, the > statement is just
            %to get the code to run)
            if ~isempty(index_list) && coordInBox > numDiffMat
                %take the difference number and randomly reassign momentum
                %to that amount of particles in this pixel
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