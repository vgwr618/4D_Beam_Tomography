clear all
close all

%initial beam paramters
qe=-1.6e-19;
beamdiv=1e-5; %beam divergence approximately
Spotsize=10e-4; %for gaussian dist
dG=1e-6; %energy spread rms control parameter.
ztime=10e-12; %bunch length in seconds
G=7;
Qbeam=2e+6*qe*0;
numstd=0.25;
npar=50000; %numbber of particles

%Quad settings (currents) arbitrary for now%
I3=1.1;
I4=-3.2;
I5=2.15;
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


%% Output data
screen1=801;

%% GPT envelope solution
figure(1)
plot(stats.d.position(2:end),stats.d.stdx(2:end),'LineWidth',2)
hold on
plot(stats.d.position(2:end),stats.d.stdy(2:end),'LineWidth',2)
legend('\sigma_x(z) (GPT)','\sigma_y(z) (GPT)','Location','best')
xlabel('z(m)')
ylabel('Spot sizes (m)')
set(gca,'FontSize',15)

%% initial
% figure(4)
% ax1=subplot(2,3,1)
% scatter(data(1).d.x*1e+6,data(1).d.Bx*1e+3,0.5)
% ax = gca;
% ax.YDir = 'normal'
% xlabel('x (um)')
% ylabel('\beta_x (mrad)')
% set(gca,'FontSize',15)
% 
% 
% ax2=subplot(2,3,2)
% scatter(data(1).d.y*1e+6,data(1).d.By*1e+3,0.5)
% ax = gca;
% ax.YDir = 'normal'
% xlabel('y (um)')
% ylabel('\beta_y (mrad)')
% set(gca,'FontSize',15)
% 
% ax3=subplot(2,3,3)
% scatter(data(1).d.y*1e+6,data(1).d.Bx*1e+3,0.5)
% ax = gca;
% ax.YDir = 'normal'
% xlabel('y (um)')
% ylabel('\beta_x (mrad)')
% set(gca,'FontSize',15)
% 
% ax5=subplot(2,3,4)
% scatter(data(1).d.x*1e+6,data(1).d.By*1e+3,0.5)
% ax = gca;
% ax.YDir = 'normal'
% xlabel('x (um)')
% ylabel('\beta_y (mrad)')
% set(gca,'FontSize',15)
% 
% ax5=subplot(2,3,5)
% scatter(data(1).d.Bx*1e+3,data(1).d.By*1e+3,0.5)
% ax = gca;
% ax.YDir = 'normal'
% xlabel('\beta_x (mrad)')
% ylabel('\beta_y (mrad)')
% set(gca,'FontSize',15)
% axis equal
% 
% ax6=subplot(2,3,6)
% scatter(data(1).d.x*1e+6,data(1).d.y*1e+6,0.5)
% ax = gca;
% ax.YDir = 'normal'
% xlabel('x (um)')
% ylabel('y (um)')
% set(gca,'FontSize',15)
% axis equal

%% xy
figure()
scatter(data(1).d.Bx*1e+3,data(1).d.By*1e+3,0.5)
ax = gca;
ax.YDir = 'normal'
xlabel('x (um)')
ylabel('\beta_x (mrad)')
set(gca,'FontSize',15)

%% %%%%%%%PHASE SPACE @ screen%%%%%%%%%%%%
figure(111)
ax1=subplot(2,3,1)
scatter(data(screen1).d.x*1e+6,data(screen1).d.Bx*1e+3,0.5)
ax = gca;
ax.YDir = 'normal'
xlabel('x (um)')
ylabel('\beta_x (mrad)')
set(gca,'FontSize',15)


ax2=subplot(2,3,2)
scatter(data(screen1).d.y*1e+6,data(screen1).d.By*1e+3,0.5)
ax = gca;
ax.YDir = 'normal'
xlabel('y (um)')
ylabel('\beta_y (mrad)')
set(gca,'FontSize',15)

ax3=subplot(2,3,3)
scatter(data(screen1).d.y*1e+6,data(screen1).d.Bx*1e+3,0.5)
ax = gca;
ax.YDir = 'normal'
xlabel('y (um)')
ylabel('\beta_x (mrad)')
set(gca,'FontSize',15)

ax5=subplot(2,3,4)
scatter(data(screen1).d.x*1e+6,data(screen1).d.By*1e+3,0.5)
ax = gca;
ax.YDir = 'normal'
xlabel('x (um)')
ylabel('\beta_y (mrad)')
set(gca,'FontSize',15)

ax5=subplot(2,3,5)
scatter(data(screen1).d.Bx*1e+3,data(screen1).d.By*1e+3,0.5)
ax = gca;
ax.YDir = 'normal'
xlabel('\beta_x (mrad)')
ylabel('\beta_y (mrad)')
set(gca,'FontSize',15)
axis equal

ax6=subplot(2,3,6)
scatter(data(screen1).d.x*1e+6,data(screen1).d.y*1e+6,0.5)
ax = gca;
ax.YDir = 'normal'
xlabel('x (um)')
ylabel('y (um)')
set(gca,'FontSize',15)
axis equal

%% xy
figure(11111)
scatter(data(end).d.x*1e+6,data(end).d.y*1e+6,0.5)
ax = gca;
ax.YDir = 'normal'
xlabel('x (um)')
ylabel('y (um)')
set(gca,'FontSize',15)
axis equal

%% projections
figure(1111)
ax1=subplot(1,2,1)
h=histogram(data(screen1).d.x*1e+6,1000)
h.EdgeColor = 'none';
hold on
xfits=h.BinEdges;
xfits(end)=[];
yfits=h.Values;
plot(xfits,yfits)

ax1=subplot(1,2,2)
h=histogram(data(screen1).d.x*1e+6,1000)
h.EdgeColor = 'none';
hold on
xfits=h.BinEdges;
xfits(end)=[];
yfits=h.Values;
plot(xfits,yfits)

%%
%need to include code here that calculates the transport matrix and phase
%advance. Must check that transport applied to initial phase space
%distribution yields final phase space in agreement with GPT.

% generate uniform rxy dist.
% n = npar;
% R = Spotsize;
% % Now create the set of points.
% t = 2*pi*rand(n,1);
% r = R*sqrt(rand(n,1));
% x = r.*cos(t);
% y = r.*sin(t);
% %generate GBx GBy
GB = sqrt(G^2-1);
sig = beamdiv*GB;
% mu = 0;
% ll = -3*sig;
% ul = 3*sig;
% pd = makedist('Normal');
% 
% gaussian = makedist('Normal','mu', 0,'sigma',sig);
% t = truncate(gaussian,ll,ul);
% Bx = random(t,n,1);
% By = random(t,n,1);
% 
% figure(2)
% ax1=subplot(2,3,1)
% scatter(x*1e+6,Bx*1e+3,0.5)
% ax = gca;
% ax.YDir = 'normal'
% xlabel('x (um)')
% ylabel('\beta_x (mrad)')
% set(gca,'FontSize',15)
% 
% 
% ax2=subplot(2,3,2)
% scatter(y*1e+6,By*1e+3,0.5)
% ax = gca;
% ax.YDir = 'normal'
% xlabel('y (um)')
% ylabel('\beta_y (mrad)')
% set(gca,'FontSize',15)
% 
% ax3=subplot(2,3,3)
% scatter(y*1e+6,Bx*1e+3,0.5)
% ax = gca;
% ax.YDir = 'normal'
% xlabel('y (um)')
% ylabel('\beta_x (mrad)')
% set(gca,'FontSize',15)
% 
% ax5=subplot(2,3,4)
% scatter(x*1e+6,By*1e+3,0.5)
% ax = gca;
% ax.YDir = 'normal'
% xlabel('x (um)')
% ylabel('\beta_y (mrad)')
% set(gca,'FontSize',15)
% 
% ax5=subplot(2,3,5)
% scatter(Bx*1e+3,By*1e+3,0.5)
% ax = gca;
% ax.YDir = 'normal'
% xlabel('\beta_x (mrad)')
% ylabel('\beta_y (mrad)')
% set(gca,'FontSize',15)
% axis equal
% 
% ax6=subplot(2,3,6)
% scatter(x*1e+6,y*1e+6,0.5)
% ax = gca;
% ax.YDir = 'normal'
% xlabel('x (um)')
% ylabel('y (um)')
% set(gca,'FontSize',15)
% axis equal

x = data(1).d.x;
y = data(1).d.y;
Bx = data(1).d.Bx;
By = data(1).d.By;
%% Transport Matrix
% init = [x' 
%     Bx'
%     y'
%     By'];
% 
% l = 0.0768;
% k03 = sqrt(0.45*I3);
% k04 = sqrt(0.45*I4);
% k05 = sqrt(0.45*I5);
% % k06 = sqrt(0.45*I6);
% 
% quad3 = [cosh(abs(k03)*l) 1/abs(k03)*sinh(abs(k03)*l) 0 0;
%     abs(k03)*sinh(abs(k03)*l) cosh(abs(k03)*l) 0 0;
%     0 0 cos(k03*l) 1/k03*sin(k03*l);
%     0 0 -k03*sin(k03*l) cos(k03*l)]
% quad4 = [cos(k04*l) 1/k04*sin(k04*l) 0 0;
%     -k04*sin(k04*l) cos(k04*l) 0 0;
%     0 0 cosh(abs(k04)*l) 1/abs(k04)*sinh(abs(k04)*l);
%     0 0 abs(k04)*sinh(abs(k04)*l) cosh(abs(k04)*l)]
% quad5 = [cosh(abs(k05)*l) 1/abs(k05)*sinh(abs(k05)*l) 0 0;
%     abs(k05)*sinh(abs(k05)*l) cosh(abs(k05)*l) 0 0;
%     0 0 cos(k05*l) 1/k05*sin(k05*l);
%     0 0 -k05*sin(k05*l) cos(k05*l)]
% %foc6 = [cos(k06*l) 1/k06*sin(k06*l);
% %    -k06*sin(k06*l) cos(k06*l)]
% ld1 = 3.2945-3.191-l/2
% ld2 = 3.3805-3.191-ld1-l-l/2
% ld3 = 3.466-3.191-ld1-l-ld2-l-l/2
% ld4 = (screen1-1)*0.0025-ld1-l-ld2-l-ld3-l
% drift1 = [1 ld1 0 0;
%     0 1 0 0;
%     0 0 1 ld1;
%     0 0 0 1]
% drift2 = [1 ld2 0 0;
%     0 1 0 0;
%     0 0 1 ld2;
%     0 0 0 1]
% drift3 = [1 ld3 0 0;
%     0 1 0 0;
%     0 0 1 ld3;
%     0 0 0 1]
% drift4 = [1 ld4 0 0;
%     0 1 0 0;
%     0 0 1 ld4;
%     0 0 0 1]

%% 
% final = drift4*quad5*drift3*quad4*drift2*quad3*drift1*init;
% flattenr = reshape(final.',1,[]);
% newx = flattenr(1:50000);
% newGBx = flattenr(50001:100000);
% newy = flattenr(100001:150000);
% newGBy = flattenr(150001:200000);
% 
% figure(6)
% ax1=subplot(2,3,1)
% scatter(newx*1e+6,newGBx*1e+3,0.5)
% ax = gca;
% ax.YDir = 'normal'
% xlabel('x (um)')
% ylabel('\beta_x (mrad)')
% set(gca,'FontSize',15)
% 
% ax2=subplot(2,3,2)
% scatter(newy*1e+6,newGBy*1e+3,0.5)
% ax = gca;
% ax.YDir = 'normal'
% xlabel('y (um)')
% ylabel('\beta_y (mrad)')
% set(gca,'FontSize',15)
% 
% ax3=subplot(2,3,3)
% scatter(newy*1e+6,newGBx*1e+3,0.5)
% ax = gca;
% ax.YDir = 'normal'
% xlabel('y (um)')
% ylabel('\beta_x (mrad)')
% set(gca,'FontSize',15)
% 
% ax5=subplot(2,3,4)
% scatter(newx*1e+6,newGBy*1e+3,0.5)
% ax = gca;
% ax.YDir = 'normal'
% xlabel('x (um)')
% ylabel('\beta_y (mrad)')
% set(gca,'FontSize',15)
% 
% ax5=subplot(2,3,5)
% scatter(newGBx*1e+3,newGBy*1e+3,0.5)
% ax = gca;
% ax.YDir = 'normal'
% xlabel('\beta_x (mrad)')
% ylabel('\beta_y (mrad)')
% set(gca,'FontSize',15)
% axis equal
% 
% ax6=subplot(2,3,6)
% scatter(newx*1e+6,newy*1e+6,0.5)
% ax = gca;
% ax.YDir = 'normal'
% xlabel('x (um)')
% ylabel('y (um)')
% set(gca,'FontSize',15)
% axis equal
% 
% figure(7)
% scatter(newx*1e+6,newy*1e+6,0.5)
% ax = gca;
% ax.YDir = 'normal'
% xlabel('x (um)')
% ylabel('y (um)')
% set(gca,'FontSize',15)
% axis equal

%% Simulate Final Yag Screen
u=zeros(2,length(data(end).d.x));
x2ps=data(end).d.x;
y2ps=data(end).d.y;
% ix2ps=data(end).d.x;
% iy2ps=data(end).d.y;
% rot = [cos(pi/2) -sin(pi/2);
%     sin(pi/2) cos(pi/2)];
% rotfin = rot * [ix2ps'; iy2ps'];
% x2ps = rotfin(1,:);
% y2ps = rotfin(2,:);

fileloc='C:\Users\vgwr6\Desktop\UCLA\soph\Musumeci_Lab\GPT\Test_Images\'
resfactor=1;
saveImage=1;
pixcal=13.614e-6; %old calibration? not really sure where it came from.
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
     radius=1; %pixels
     filter=fspecial('disk',radius);
     yag_Image=conv2(yag_Image,filter);
     yag_Image = imresize(yag_Image,[800 800]);
     yag_Image=uint8(yag_Image);
     imageName = strcat('GPT_image_init2_1.1_-3.2_2.15','.bmp');
     imwrite(yag_Image,strcat(fileloc,imageName));
end
final_gpt_image=imread('C:\Users\vgwr6\Desktop\UCLA\soph\Musumeci_Lab\GPT\Test_Images\GPT_image_init2_1.1_-3.2_2.15.bmp');
final_gpt_image = 1/sum(final_gpt_image, 'all')*double(final_gpt_image);
sum(final_gpt_image, 'all')
figure()
imagesc(final_gpt_image)
colorbar
%image is rotated or transposed??

%% initial YAG screen
u=zeros(2,length(data(1).d.x));
x2ps=data(1).d.x;
y2ps=data(1).d.y;
% ix2ps=data(1).d.x;
% iy2ps=data(1).d.y;
% rot = [cos(pi/2) -sin(pi/2);
%     sin(pi/2) cos(pi/2)];
% rotfin = rot * [ix2ps'; iy2ps'];
% x2ps = rotfin(1,:);
% y2ps = rotfin(2,:);

fileloc='C:\Users\vgwr6\Desktop\UCLA\soph\Musumeci_Lab\GPT\Test_Images\'
resfactor=1;
saveImage=1;
pixcal=13.614e-6; % 13.614e-6 old calibration? not really sure where it came from.
if saveImage==1
     yag_Image_init = zeros(floor(800/resfactor),floor(800/resfactor));
     for l=1:length(x2ps)
         xpos=floor(x2ps(l)/(resfactor*pixcal)+400/resfactor);
         ypos=floor(-y2ps(l)/(resfactor*pixcal)+400/resfactor);
         if xpos>0 && xpos<=800/resfactor && ypos>0 && ypos<=800/resfactor
             yag_Image_init(ypos,xpos)=yag_Image_init(ypos,xpos)+1;
         end
     end
%      yag_Image_init(400,401)
     sum(yag_Image_init,'all');
     %renormalize and smooth image
     maxIm=max(max(yag_Image_init));
     yag_Image_init=yag_Image_init.*floor(400/maxIm);
     radius=1; %pixels
     filter=fspecial('disk',radius);
     yag_Image_init=conv2(yag_Image_init,filter);
     yag_Image_init = imresize(yag_Image_init,[800 800]);
     yag_Image_init=uint8(yag_Image_init);
     imageName = strcat('GPT_image_init2_1.1_-3.2_2.15','.bmp');
     imwrite(yag_Image_init,strcat(fileloc,imageName));
end
init_gpt_image=imread('C:\Users\vgwr6\Desktop\UCLA\soph\Musumeci_Lab\GPT\Test_Images\GPT_image_init2_1.1_-3.2_2.15.bmp');
% maxInitIm=max(max(init_gpt_image));
% init_gpt_image=init_gpt_image.*floor(400/maxIm);
% init_gpt_image = floor((npar / sum(init_gpt_image, 'all')) * init_gpt_image);
% init_gpt_image = floor((npar / sum(init_gpt_image, 'all')) * init_gpt_image);
init_gpt_image = uint8(105000/sum(init_gpt_image, 'all')*double(init_gpt_image));
numpart = sum(init_gpt_image, 'all')
figure()
imagesc(init_gpt_image)
colorbar

% figure()
% imagesc(yag_Image_init)
% colorbar

%%
initmatrixName = strcat('initImage_1.1_-3.2_2.15.csv')
% init_mat = floor((1 / sum(init_gpt_image,'all'))* (npar * init_gpt_image));
writematrix(init_gpt_image, strcat(fileloc, initmatrixName))
init_mat = readmatrix('C:\Users\vgwr6\Desktop\UCLA\soph\Musumeci_Lab\GPT\Test_Images\init2Image_1.1_-3.2_2.15.csv');
figure()
imagesc(init_mat)

finmatrixName = strcat('finImage_1.1_-3.2_2.15.csv')
writematrix(final_gpt_image, strcat(fileloc, finmatrixName))
fin_mat = readmatrix('C:\Users\vgwr6\Desktop\UCLA\soph\Musumeci_Lab\GPT\Test_Images\fin2Image_1.1_-3.2_2.15.csv');
figure()
imagesc(fin_mat)

sum(init_mat,'all')
sum(fin_mat,'all')
%%
[init4DCoord, fin4DCoord, finIm, finImDisp, diffMat, posCoord, negCoord, sigmaRecon, stReconMeasXY] = forward_propagation(init_gpt_image, final_gpt_image, R, numpart, 3e-20, pixcal, resfactor);
figure()
imagesc(finImDisp)
colorbar
figure()
imagesc(diffMat)
colorbar

%% One iteration
[init4DCoord, fin4DCoord, finIm, finImDisp, diffMat, posCoord, negCoord, sigmaRecon, stReconMeasXY] = first_propagation(init_gpt_image, final_gpt_image, R, numpart, 3e-20, pixcal, resfactor);
figure()
imagesc(finImDisp)
colorbar
figure()
imagesc(diffMat)
colorbar

figure
histogram(init4DCoord(1,:),100)
hold on
histogram(data(1).d.x,100)

figure
histogram(init4DCoord(3,:),100)
hold on
histogram(data(1).d.y, 100)
%% first iteration
[init4DCoordOld, fin4DCoord, finIm, finImDisp, diffMat, posCoord, negCoord, sigmaRecon, stReconMeasXY] = first_propagation(init_gpt_image, final_gpt_image, R, numpart, 3e-20, pixcal, resfactor);
sum(diffMat, 'all')
%% algorithm loop
for i = 1 : 5
    init4DCoord = backward_propagation(init4DCoordOld, fin4DCoord, diffMat, pixcal, resfactor);
    [init4DCoordOld, fin4DCoord, finIm, finImDisp, diffMat] = forward_propagation(init4DCoord, final_gpt_image, R, pixcal, resfactor);
    figure()
    imagesc(finImDisp)
    colorbar
    figure()
    imagesc(diffMat)
    colorbar
    ConvHist = []
    ConvHist(end+1) = sum(diffMat, 'all')
    i
end
%% ALL BELOW IS SCRATCH

figure
histogram(init4DCoordOld(2,:),100)
figure
histogram(data(1).d.Bx,100)

figure
histogram(init4DCoordOld(4,:),100)
hold on
histogram(data(1).d.By, 100)
%%
diffmatrixName = strcat('diffMat_1.1_-3.2_2.15.csv')
writematrix(diffMat, strcat(fileloc, diffmatrixName))
%%
init_gpt_image1 = init_gpt_image;
sum(init_gpt_image1, 'all')
sum(final_gpt_image, 'all')
% sum(init_gpt_image, 'all')
%%
posX = posCoord(1,:);
posY = posCoord(2,:);
figure()
scatter(posX,posY,0.5)
ax = gca;
ax.YDir = 'normal'
xlabel('\beta_x (mrad)')
ylabel('\beta_y (mrad)')
set(gca,'FontSize',15)
axis equal

negX = negCoord(1,:);
negY = negCoord(2,:);
figure()
scatter(negX,negY,0.5)
ax = gca;
ax.YDir = 'normal'
xlabel('\beta_x (mrad)')
ylabel('\beta_y (mrad)')
set(gca,'FontSize',15)
axis equal

%%
figure
itpPosX = interp1(posX,posY,fin4DCoord(1,:)');
plot(posX,posY,'o',fin4DCoord(1,:),itpPosX ,':.');
title('(Default) Linear Interpolation');
%%
x = [1 2 3 4]
y=[0 4 3 4]
size(x)
size(y)
figure()
scatter(x,y,3.5)