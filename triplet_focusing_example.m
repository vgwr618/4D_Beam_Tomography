function [ R , R_all, z_all ,pos] = triplet_focusing_example(quads,gammaBeta,L,position)
% This function calculates the transport matrix from screen4 to the
% microscope objective screen as a function of quad4 and quad5 currents
% taking into account the real field maps (input by tanh function)
bRho=gammaBeta*0.0017 ;
leff=0.0768;
qgrad4=0.45;
qgrad5=0.45;
qgrad6=0.45;


dz = 0.00001;
R = eye(4);
b = 135;

pos_screen4 = position; % 0.555; % position screen4
ds4 = 0.1035; % screen4 to Q4
d45 = 0.1895; % Center Q4 to center Q5 0.086
d56 = 0.2750; % Center Q5 to center Q6 0.085
total_length = L ;

%quad strengths
k = [...
    quads(1)*qgrad4/bRho ...
    quads(2)*qgrad5/bRho ...
    quads(3)*qgrad6/bRho ...
    ];    

%quad positions
pos = pos_screen4 +[...
    ds4 ...
    d45 ...
    d56 ...
    ];
z_all = pos_screen4 + (0:dz:total_length);
R_all = zeros(4,4,length(z_all)); i=1;



for z = z_all
grad = k(1)/2*(tanh(b/2*(leff/2-(z-pos(1))))+tanh(b/2*(leff/2+(z-pos(1)))));
grad = grad + k(2)/2*(tanh(b/2*(leff/2-(z-pos(2))))+tanh(b/2*(leff/2+(z-pos(2)))));
grad = grad + k(3)/2*(tanh(b/2*(leff/2-(z-pos(3))))+tanh(b/2*(leff/2+(z-pos(3)))));


SCx=0;
SCy=0;

    if abs(grad) < 1e-4
        Rdelta = [1 dz 0 0;0 1 0 0;0 0 1 dz;0 0 0 1];
    else

    sKx=sqrt(grad-SCx);
    sKLx=sKx*dz;
    sKy=sqrt(grad-SCy);
    sKLy=sKy*dz;

    Rdelta= [cos(sKLx) sin(sKLx)/sKx 0 0;-sKx*sin(sKLx) cos(sKLx) 0 0;0 0 cosh(sKLy) sinh(sKLy)/sKy;0 0 sKy*sinh(sKLy) cosh(sKLy)];
    end
    R = Rdelta*R;   
    R_all(:,:,i)=R;
    i=i+1;
end

end