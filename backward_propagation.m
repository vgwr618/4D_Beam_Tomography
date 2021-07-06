function [init4DCoord] = backward_propagation(init4DCoord, fin4DCoord, diffMat, pixcal, resfactor, sigma, R)
pixelDim = size(diffMat);
%debugging code
pos_count = 0;
neg_count = 0;
            
% sign_mat = sign(diffMat);
% negXCoord = zeros(1, sum(sign_mat(:)==-1)));
% negYCoord = zeros(1, sum(sign_mat(:)==-1)));
total_neg = abs(floor(sum(sum(diffMat(diffMat < 0)))));
negXCoord = zeros(1, total_neg);
negYCoord = zeros(1, total_neg);
% negXCoord = [];
% negYCoord = [];
for i = 1 : pixelDim(1)
    for j = 1 : pixelDim(2)
        %right and bottom boundary lines for each pixel
        pixLeftX = pixcal * ((resfactor * j) - (pixelDim(2) / 2));
        pixTopY = pixcal * ((resfactor * -i) + (pixelDim(1) / 2));
        if diffMat(i,j) < 0
            for k = 1:abs(floor(diffMat(i,j)))
                first_zero = find(negXCoord==0, 1, 'first');
                negXCoord(first_zero) = pixLeftX + rand(1,1)* pixcal * resfactor;
                negYCoord(first_zero) = pixTopY - rand(1,1)* pixcal * resfactor;   
            end
        end
    end
end
rand_coord_index = randperm(length(negXCoord));
        
for i = 1 : pixelDim(1)
    for j = 1 : pixelDim(2)
        %right and bottom boundary lines for each pixel
        pixLeftX = pixcal * ((resfactor * j) - (pixelDim(2) / 2));
        pixTopY = pixcal * ((resfactor * -i) + (pixelDim(1) / 2));
        if diffMat(i,j) > 0
            index_list = [];
            %check for final simulated particles that are in each pixel box
            for k = 1:length(fin4DCoord(1,:))
                if fin4DCoord(1,k) >= pixLeftX && fin4DCoord(1,k) < pixLeftX + resfactor * pixcal && fin4DCoord(3,k) <= pixTopY && fin4DCoord(3,k) > pixTopY - resfactor * pixcal
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
            if ~isempty(index_list)
%                 index_list;
                %take the difference number and randomly reassign momentum
                %to that amount of particles in this pixel
                for m = 1:floor(diffMat(i,j))
                    n = index_list(m);
%                     fin4DCoord(1,n)
%                     fin4DCoord(3,n)
%                     pixLeftX
%                     pixTopY
                    init4DCoord(2, n) = (negXCoord(rand_coord_index(1)) - R(1,1) * init4DCoord(1, n)) / R(1,2);
                    init4DCoord(4, n) = (negYCoord(rand_coord_index(1)) - R(3,3) * init4DCoord(3, n)) / R(3,4);
                    rand_coord_index = rand_coord_index([2:end]);
%                     init4DCoord(2, n) = -3 * sigma + rand(1,1)* 6 * sigma;
%                     init4DCoord(4, n) = -3 * sigma + rand(1,1)* 6 * sigma;
                end
            end
        end
    end
end
neg_count
pos_count
length(rand_coord_index)