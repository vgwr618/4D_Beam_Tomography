for i = 1:0.2:3.6
    for j = 1:0.2:3.6
        for k = 1:0.2:3.6
            qe=-1.6e-19;
            beamdiv=1e-5; %beam divergence approximately
            Spotsize=2e-4; %for gaussian dist
            dG=1e-6; %energy spread rms control parameter.
            ztime=10e-12; %bunch length in seconds
            G=7;
            Qbeam=2e+6*qe*0;
            numstd=0.25;
            npar=50000;
            I3 = i;
            I4 = -j;
            I5 = k;
            I6 = 0;

            infile = 'Tomography_test2.in';
            gdffile = 'Tomography_test_solution.gdf';
            path_to_gpt='C:\Users\vgwr6\Desktop\UCLA\soph\Musumeci_Lab\General_Particle_Tracer\bin\'; %Must adapt this to your pc
            path_to_in_out_file=' C:\Users\vgwr6\Desktop\UCLA\soph\Musumeci_Lab\GPT\' ; %Must adapt this to your pc
            license=' GPTlicense=1476385047';

            callGPT = strcat(path_to_gpt,'gpt -v -o ',path_to_in_out_file,gdffile,path_to_in_out_file,infile,...
            ' Qtot=',num2str(Qbeam),' npar=',num2str(npar),' Spotsize=',num2str(Spotsize),' beamdiv=',num2str(beamdiv),' G=',num2str(G),' dG=',num2str(dG),...
            ' ztime=',' I3=',num2str(I3),' I4=',num2str(I4),' I5=',num2str(I5),' I6=',num2str(I6),num2str(ztime),license)

            system(callGPT,'-echo') ;
            data=load_gdf(strcat(path_to_in_out_file,gdffile)) ; 
            callGPT2=strcat(path_to_gpt,'gdfa -o',path_to_in_out_file,'result.gdf ',path_to_in_out_file,gdffile,' position ',' stdBy ',' stdBx ',' stdBz',' stdx ',' stdy',' CSgammax',' CSgammay',' CSalphax',' CSalphay', ' CSbetax',' CSbetay',' nemixrms',' nemiyrms',' nemirrms',' nemizrms',' avgz',' stdz',' avgG',' stdG');
            system(callGPT2,'-echo') ;
            stats=load_gdf(strcat(path_to_in_out_file,'result.gdf')); 

            x2ps=data(end).d.x;
            y2ps=data(end).d.y;
            fileloc='D:\PBPL\tomography_gpt_images\';
            resfactor=1;
            saveImage=1;
            pixcal=27e-6; %13.614e-6 old calibration? not really sure where it came from.
            if saveImage==1
                 yag_Image = zeros(floor(800/resfactor),floor(800/resfactor));
                 for l=1:length(x2ps)
                     xpos=floor(x2ps(l)/(resfactor*pixcal)+400/resfactor);
                     ypos=floor(-y2ps(l)/(resfactor*pixcal)+400/resfactor);
                     if xpos>0 && xpos<=800/resfactor && ypos>0 && ypos<=800/resfactor
                         yag_Image(ypos,xpos)=yag_Image(ypos,xpos)+1;
                     end
                 end

                 %renormalize and smooth image
                 maxIm=max(max(yag_Image));
                 yag_Image=yag_Image.*floor(400/maxIm);
                 radius=3; %pixels
                 filter=fspecial('disk',radius);
                 yag_Image=conv2(yag_Image,filter);
                 yag_Image = imresize(yag_Image,[800 800]);
                 yag_Image=uint8(yag_Image);
                 formatSpec = 'GPT_image_fin_%.1f_%.1f_%.1f';
                 imageName = strcat(sprintf(formatSpec, I3, I4, I5),'.bmp');
                 imwrite(yag_Image,strcat(fileloc,imageName));
            end

            %%%WARNING, RUN THIS AFTER YOU HAVE RUN Tomgraphy_testfun%%%%
            screen=1;
            % quads=[1.1,-3.2,2.15];%quad currents
            quads=[-I3,-I4,-I5];
            G=7; %gamma factor
            gammaBeta=sqrt(G^2-1); %gamma times beta
            L=2; %distance from 3rd quad center to final screen.
            position=0; %initial screen position


            %%%triplet_focusing_example is the matrix propagator%%%%
            [ R , R_all, z_all ,pos] = triplet_focusing_example(quads,gammaBeta,L,position);
            formatSpec = 'R_%.1f_%.1f_%.1f';
            matrixName = strcat(sprintf(formatSpec, I3, I4, I5),'.csv');
            writematrix(R, strcat(fileloc, matrixName))
        end
    end
end