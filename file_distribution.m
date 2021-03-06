%% momentum double gaussian
sd = 1e-4;
mu1 = -2*sd;
x1 = linspace(-4*sd,0,1000);
Px1 = 1/(2*pi*sd)*exp(-(x1-mu1).^2/(2*sd^2));
x2 = linspace(1e-10,4*sd,1000);
mu2 = 2*sd;
Px2 = 1/(2*pi*sd)*exp(-(x2-mu2).^2/(2*sd^2));
x = horzcat(x1,x2);
Px = horzcat(Px1,Px2);
Px = Px/norm(Px);
x = x';
Px = Px';
T = table(x, Px);
writetable(T,'pxInput.txt','Delimiter','\t');
%% spatial double gaussian
sd = 1e-4;
mu1 = -2*sd;
x1 = linspace(-4*sd,0,1000);
Px1 = 1/(2*pi*sd)*exp(-(x1-mu1).^2/(2*sd^2));
x2 = linspace(1e-10,4*sd,1000);
mu2 = 2*sd;
Px2 = 1/(2*pi*sd)*exp(-(x2-mu2).^2/(2*sd^2));
x = horzcat(x1,x2);
Px = horzcat(Px1,Px2);
Px = Px/norm(Px);
x = x';
Px = Px';
T = table(x, Px);
writetable(T,'xInput.txt','Delimiter','\t');