%MAIN
%SEVERALMEDIA: this is the main body of the calculation, all the information about the different media and the NA of the imaging system needs to be inserted in this code. NOTE: At some point I was trying to parallelize this code, therefore it has extra comments that are not needed
save('two170coverslip.mat')
%for w = 1:1:12
%w = 1;
 %  load('results.mat');
%matlabpool open 3
%options = optimset('UseParallel','always')
%insert variables
zcenter = [34,35,36,38,39,40,41,39,40,42,43,44];
global M n h k0 k rp thetap;
M = 5; %insert the number of interfaces
n(1:M) = [1.33, 1.33, 1.55, 1.33, 1.35] ; %insert the refractive indeces from left to right
h(1:M-1) = [200 + 170 + 100 + 170, 200 + 170 + 100, 200 + 100, 200]; %insert the distances from the origin to the interfaces from left to right (units in microns)
lambda = 0.940; %insert the wavelength of light (units in microns)
NA = 0.8; %insert the numerical aperture
f = 1; %insert the focal lenght of the objective
lo = 1; %weight constant
xmin = 0; %insert the initial value on x (perpendicular axis)
xmax = 3; %insert the final value on x
xnumb = 31; %insert number of steps x
zmin = 33 - 15; %insert the initial value on z (propagation axis)
zmax = 33 + 15; %insert the final value on z
znumb = 301;

%calculate constants
k0 = 2*pi/lambda; %calculate the wave number in vaccuum
k(1:M) = n(1:M)*k0 ; %calculate the wave number in every material
alpha = asin(NA/n(1)); %calculate the angle of the aperture
Kbig = pi*n(1)*f*lo/lambda ; %calculate weight constant
%I = ones(znumb, xnunb);

%initiate calculation
for x = xmin: (xmax-xmin)/xnumb :xmax
    for z = zmin: (zmax-zmin)/znumb :zmax
        
        
        rp = sqrt(x^2 + z^2);
        if z >= 0
            thetap = atan(x/z);
        end
        if z < 0
            thetap = pi +atan(x/z);
        end
        phyp = 0; %this needs to be changed, but for now it works
        I0N = quad(@(theta1)serothbesselfunct(theta1),0,alpha,0.000001);
        I1N = quad(@(theta1)firstbesselfunct(theta1),0,alpha,0.000001);
        I2N = quad(@(theta1)secondbesselfunct(theta1),0,alpha,0.000001);

        ENx = -1i*Kbig*( I0N + I2N*cos(2*phyp)) ; %electric fiedl in x
        ENy = -1i*Kbig*I2N*sin(2*phyp) ; %electric fiedl in y
        ENz = -2*Kbig*I1N*cos(phyp) ; %electric fiedl in z
        if abs(x) < 0.001 %I have a problem at the point x=0 z=0 because the functions end up being divided by infinity
            if abs(z) < 0.001
                ENx = 0; %electric fiedl in x
                ENy = 0; %electric fiedl in y
                ENz = 0; %electric fiedl in z
            end
        end
        I(int64(((z-zmin)*znumb/(zmax-zmin)) +1),int64(((x-xmin)*xnumb/(xmax-xmin)) +1)) = ENx*conj(ENx)+ENy*conj(ENy)+ENz*conj(ENz);
    end
end

load('two170coverslip.mat')
number = 1;%w;
number = int32(number);
airpocket(number) = h(2) - h(3);
figure
imagesc(I)
%title(airpocket(number))
[numberx(number),maximx] = max(max(abs(I)));
[numberz(number),maximz] = max(max(abs(I')));
maximumdepthposition(number) = maximz*(zmax-zmin)/znumb +zmin;
[numberx(number),numberz(number),maximx,maximumdepthposition(number),airpocket(number),zmin, zmax]
Ifinal(:,:,number) = I(:,:);
zetamin(number) = zmin;
zetamax(number) = zmax;
equizmin(number) = xmin;
equizmax(number) = xmax;
% Ifinal(:,:,number) = Inorm;
% Psffinal(:,:,number) = psf2norm;
save('two170coverslip.mat','Ifinal','airpocket','maximumdepthposition','zetamin','zetamax','equizmin','equizmax')
%print results
%save('results.mat','results'); 
%options = optimset('UseParallel','never')
%matlabpool close
clear all %havent figured out, but something is being recorded
%end