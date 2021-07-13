rotX = [];
rotY = [];
for i = 1.0:0.4:1.0
    for j = 1.0:0.4:1.0
        for k = 1.0:0.4:1.0
            I3 = -i;
            I4 = j;
            I5 = -k;
            formatSpec = 'R_%.1f_%.1f_%.1f';
            matrixName = strcat('D:\PBPL\gpt_images_0.4_increment\',sprintf(formatSpec, I3, I4, I5),'.csv');
            R = readmatrix(matrixName);
            S=[1 0 0 0; 
              (R(1,2)*R(2,2)+R(1,1)*R(2,1))/(R(1,1)^2+R(1,2)^2) 1 0 0;
               0 0 1 0;
               0 0 (R(3,4)*R(4,4)+R(3,3)*R(4,3))/(R(3,3)^2+R(3,4)^2) 1];
            E=[sqrt(R(1,1)^2+R(1,2)^2) 0 0 0;
                0 1/sqrt(R(1,1)^2+R(1,2)^2) 0 0;
                0 0 sqrt(R(3,3)^2+R(3,4)^2) 0;
                0 0 0 1/sqrt(R(3,3)^2+R(3,4)^2)];
            Rot=[R(1,1)/sqrt(R(1,1)^2+R(1,2)^2) R(1,2)/sqrt(R(1,1)^2+R(1,2)^2) 0 0;
                -R(1,2)/sqrt(R(1,1)^2+R(1,2)^2) R(1,1)/sqrt(R(1,1)^2+R(1,2)^2) 0 0; ...
                0 0 R(3,3)/sqrt(R(3,3)^2+R(3,4)^2) R(3,4)/sqrt(R(3,3)^2+R(3,4)^2);
                0 0 -R(3,4)/sqrt(R(3,3)^2+R(3,4)^2) R(3,3)/sqrt(R(3,3)^2+R(3,4)^2)];
            e1=E(1,1);
            e3=E(3,3);
            shearx=S(2,1);
            sheary=S(4,3);
            rotationsx=acos(Rot(1,1));
            rotationsy=acos(Rot(3,3));
            format = '\n%.1f_%.1f_%.1f: rotX = %.4f, rotY = %.4f';
            fprintf(format, i,j,k,rotationsx, rotationsy)
            rotX(end+1) = rotationsx;
            rotY(end+1) = rotationsy;
        end
    end
end


%%
figure()
scatter(rotX, rotY, 3.5, 'filled')