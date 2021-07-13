guess_sigma = 1e-4;
numImages = 343;
init_gpt_image=imread('C:\Users\vgwr6\Desktop\UCLA\soph\Musumeci_Lab\GPT\GPT_Images\GPT_image_init_-1.5_2.0_-1.5.bmp');
init_gpt_image = uint8(100000/sum(init_gpt_image, 'all')*double(init_gpt_image));
numpart = sum(init_gpt_image, 'all');
initXMomentum = zeros(1, numpart);
initYMomentum = zeros(1, numpart);
[init4DCoord, fin4DCoord, finImFirst, finImDispFirst, diffMat, sigmaRecon, stReconMeasXY, newpart] = first_propagation(init_gpt_image, final_gpt_image, R, numpart, guess_sigma, pixcal, resfactor);
for n = 1:25
    for i = 1.0:0.4:1.0
        for j = 1.0:0.4:1.0
            for k = 1.0:0.4:1.0
                I3 = -i;
                I4 = j;
                I5 = -k;

                RFormat = 'R_%.1f_%.1f_%.1f';
                matrixName = strcat('D:\PBPL\gpt_images_0.4_increment\',sprintf(RFormat, I3, I4, I5),'.csv');
                R = readmatrix(matrixName);

                ImageFormat = 'GPT_image_fin_%.1f_%.1f_%.1f';
                ImageName = strcat('D:\PBPL\gpt_images_0.4_increment\',sprintf(Imageformat, I3, I4, I5),'.csv');
                final_gpt_image=imread('C:\Users\vgwr6\Desktop\UCLA\soph\Musumeci_Lab\GPT\GPT_Images\GPT_image_fin_-1.5_2.0_-0.5.bmp');
                final_gpt_image = numpart/sum(final_gpt_image, 'all')*double(final_gpt_image);

                init4DCoordNew = backward_propagation(init4DCoord, fin4DCoord, diffMat, pixcal, resfactor, guess_sigma, R);
                [init4DCoordNew, fin4DCoord, finIm, finImDisp, diffMat] = forward_propagation(init4DCoordNew, final_gpt_image, R, pixcal, resfactor);

                initXMomentum = initXMomentum + init4DCoordNew(2,:);
                initYMomentum = initYMomentum + init4DCoordNew(4,:);
            end
        end
    end
    init4DCoord(2,:) = initXMomentum / numImages;
    init4DCoord(4,:) = initYMomentum / numImages;
end
            
            