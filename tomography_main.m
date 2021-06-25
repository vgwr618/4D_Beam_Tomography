clear all
close all

%initial beam paramters
qe=-1.6e-19;
beamdiv=1e-4; %beam divergence approximately
Spotsize=1e-4; %for gaussian dist
dG=1e-6; %energy spread rms control parameter.
ztime=10e-12; %bunch length in seconds
G=7;
Qbeam=2e+6*qe*0;
numstd=0.25;
npar=50000; %numbber of particles

%Quad settings (currents) arbitrary for now%
I3=-0.5;
I4=2.2;
I5=-3.5;
I6=5*0;

%what are the phase advances for arbitrary settings?

%% %%%%%%%%%%%%%%%%%%%%%%%%%Run GPT%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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

%% Initial Yag Screen
u=zeros(2,length(data(1).d.x));
x2ps=data(1).d.x;
y2ps=data(1).d.y;

fileloc='C:\Users\vgwr6\Desktop\UCLA\soph\Musumeci_Lab\GPT\Test_Images\'
resfactor=1;
saveImage=1;
pixcal=27e-6; % 13.614e-6 old calibration? not really sure where it came from.
if saveImage==1
     yag_Image_init = zeros(floor(800/resfactor),floor(800/resfactor));
     for l=1:length(x2ps)
         xpos=floor(x2ps(l)/(resfactor*pixcal)+400/resfactor);
         ypos=floor(-y2ps(l)/(resfactor*pixcal)+400/resfactor);
         if xpos>=0 && xpos<800/resfactor && ypos>=0 && ypos<800/resfactor
             yag_Image_init(ypos,xpos)=yag_Image_init(ypos,xpos)+1;
         end
     end
     %renormalize and smooth image
     maxIm=max(max(yag_Image_init));
     yag_Image_init=yag_Image_init.*floor(400/maxIm);
     radius=1; %pixels
     filter=fspecial('disk',radius);
     yag_Image_init=conv2(yag_Image_init,filter);
     yag_Image_init=imresize(yag_Image_init,[800 800]);
     yag_Image_init=uint8(yag_Image_init);
     imageName = strcat('GPT_image_init_1.1_-3.2_2.15','.bmp');
     imwrite(yag_Image_init,strcat(fileloc,imageName));
end
init_gpt_image=imread('C:\Users\vgwr6\Desktop\UCLA\soph\Musumeci_Lab\GPT\Test_Images\GPT_image_init_1.1_-3.2_2.15.bmp');
init_gpt_image = uint8(100000/sum(init_gpt_image, 'all')*double(init_gpt_image));
numpart = sum(init_gpt_image, 'all')
figure()
imagesc(init_gpt_image)
colorbar


%% Final Yag Screen
u=zeros(2,length(data(end).d.x));
x2ps=data(end).d.x;
y2ps=data(end).d.y;

fileloc='C:\Users\vgwr6\Desktop\UCLA\soph\Musumeci_Lab\GPT\Test_Images\'
resfactor=1;
saveImage=1;
pixcal=27e-6; %13.614e-6 old calibration? not really sure where it came from.
if saveImage==1
     yag_Image = zeros(floor(800/resfactor),floor(800/resfactor));
     for l=1:length(x2ps)
         xpos=floor(x2ps(l)/(resfactor*pixcal)+400/resfactor);
         ypos=floor(-y2ps(l)/(resfactor*pixcal)+400/resfactor);
         if xpos>=0 && xpos<800/resfactor && ypos>=0 && ypos<800/resfactor
             yag_Image(ypos,xpos)=yag_Image(ypos,xpos)+1;
         end
     end
     
     %renormalize and smooth image
     maxIm=max(max(yag_Image));
     yag_Image=yag_Image.*floor(400/maxIm);
     radius=1; %pixels
     filter=fspecial('disk',radius);
     yag_Image=conv2(yag_Image,filter);
     yag_Image = imresize(yag_Image,[800 800]);
     yag_Image=uint8(yag_Image);
     imageName = strcat('GPT_image_fin2_1.1_-3.2_2.15','.bmp');
     imwrite(yag_Image,strcat(fileloc,imageName));
end

final_gpt_image=imread('C:\Users\vgwr6\Desktop\UCLA\soph\Musumeci_Lab\GPT\Test_Images\GPT_image_fin2_1.1_-3.2_2.15.bmp');
% final_gpt_image = numpart/sum(final_gpt_image, 'all')*double(final_gpt_image);
figure()
imagesc(final_gpt_image)
colorbar

%% R Matrix
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

%% compare GPT to the matrix code
%%% need to pull data from GPT simulation to construct the initial sigma
%%% matrix.
x=data(screen).d.x;
xp=data(screen).d.Bx;
y=data(screen).d.y;
yp=data(screen).d.By;
phsp=zeros(length(x),4);
phsp(:,1)=x;
phsp(:,2)=xp;
phsp(:,3)=y;
phsp(:,4)=yp;
sigma0=cov(phsp);


%Initialize arrays to append C and S matrix elements at each z
Cpointsx=zeros(1,length(R_all));
Spointsx=zeros(1,length(R_all));
Cpointsy=zeros(1,length(R_all));
Spointsy=zeros(1,length(R_all));
%Initialize arrays to append Spotsizes at each z
stdx=zeros(1,length(R_all));
stdxp=zeros(1,length(R_all));
stdy=zeros(1,length(R_all));
stdyp=zeros(1,length(R_all));


for j=1:length(Spointsx)
Cpointsx(j)=R_all(1,1,j);
Spointsx(j)=R_all(1,2,j);
Cpointsy(j)=R_all(3,3,j);
Spointsy(j)=R_all(3,4,j);
S_aux = R_all(:,:,j)*sigma0*R_all(:,:,j)';
stdx(j) = sqrt(S_aux(1,1));
stdxp(j) = sqrt(S_aux(1,2));
stdy(j) = sqrt(S_aux(3,3));
stdyp(j) = sqrt(S_aux(3,4));
end

figure
plot(z_all,stdx,'LineWidth',2)
hold on
plot(z_all,stdy,'LineWidth',2)
hold on
plot(stats.d.position(2:end),stats.d.stdx(2:end),'--','LineWidth',2)
hold on
plot(stats.d.position(2:end),stats.d.stdy(2:end),'--','LineWidth',2)
legend('\sigma_x(z)','\sigma_y(z)','\sigma_x(z) (GPT)','\sigma_y(z) (GPT)','Location','best')
xlabel('z(m)')
set(gca,'FontSize',15)

%% first iteration
guess_sigma = 1e-4;
[init4DCoordOld, fin4DCoord, finIm, finImDisp, diffMat, sigmaRecon, stReconMeasXY, newpart] = first_propagation(init_gpt_image, final_gpt_image, R, numpart, guess_sigma, pixcal, resfactor);
%% check difference matrix's total positive intensity difference and negative intensity difference
sum(sum(diffMat(diffMat > 0)))
sum(sum(diffMat(diffMat < 0)))
sum(diffMat,'all')
%% plot first iteration final image and difference image
figure()
imagesc(finImDisp)
colorbar
figure()
imagesc(diffMat)
colorbar

%% initialize convergence history
ConvHist = [sum(sum(diffMat(diffMat > 0)))]
%% algorithm loop
for i = 1 : 5
    init4DCoord = backward_propagation(init4DCoordOld, fin4DCoord, diffMat, pixcal, resfactor, guess_sigma);
    [init4DCoordOld, fin4DCoord, finIm, finImDisp, diffMat] = forward_propagation(init4DCoord, final_gpt_image, R, pixcal, resfactor);
%     figure()
%     imagesc(finImDisp)
%     colorbar
%     figure()
%     imagesc(diffMat)
%     colorbar
    figure
    scatter(init4DCoord(2,:),init4DCoord(4,:),0.1)
    ConvHist(end+1) = sum(sum(diffMat(diffMat > 0)))
    i
end
%% plot convergence history on the scale of 0 to max possible error
partIndex = linspace(1,length(ConvHist),length(ConvHist));
figure()
scatter(partIndex,ConvHist,5,'filled')
yheight = sum(finImDisp,'all');
ylim([0 yheight])
