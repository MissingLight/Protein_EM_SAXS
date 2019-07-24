% The program calculates the Fast Fourier Transform and Discrete Fourier
% Transform for a two-dimensional Gaussian, performs an interpolation to
% and fro polar coordinates and Cartesian coordinates.

% 22/9/19
% Cartesian to Polar to Cartesin conversion method 2
% (i) Interpolate in Cartesian coordinates
% (ii) Convert to polar coordinates (define a circle larger than the square
% domain) with npts*npts points
% (iii) Rediscretise in polar coordinates, and redefine the mesh points
% (iv) Interpolate and convert to Cartesian coordinates
% (v) Rediscretise and redefine Cartesian mesh points
% (vi) Apply inverse FFT to new Cartesian domain and compare with original
% signal.

% 18/9/19
% Improving phase through: 
% (i) Removal of low amplitude components and 
% (ii) a phase shifting from 0 to N-1 to -N/2 to N/2-1

% Time discretisation
clc;clear;close all;
% Starting and end-points of the signal, i.e. t0 <= t < tend
t10 = -10;
t1end = 10;
t20 = -10;
t2end = 10;

% Sampling frequency
npts = 128; %i.e. a npts x npts grid
T1o = (t1end - t10);
T2o = (t2end - t20);
T1s = T1o/npts;
T2s = T2o/npts;
f1s = 1/T1s;
f2s = 1/T2s;
t1array = t10:T1s:(t1end-T1s);
t2array = t20:T2s:(t2end-T2s);

% Frequency discretisation
d1w = (2*pi)/T1o;
d2w = (2*pi)/T2o;
w10 = -npts*(d1w)/2;
w20 = -npts*(d2w)/2;
w1end = (npts)*d1w/2 - d1w;
w2end = (npts)*d2w/2 - d2w;
w1array = w10:d1w:w1end;
w2array = w20:d2w:w2end;

% Nyquist criterion is satisfied since 2*(wend/(2*pi)) < fs; hence no aliasing.

% Field arrays
t = zeros(1,2)';
w = zeros(1,2)';
field = zeros(npts);
analyticf = zeros(npts);
FFTfield = zeros(npts);
DFTfield = zeros(npts);
inverse2 = zeros(npts);

% Input the space field at each point
field = smile(npts,t1array,t2array);

% Analytic Fourier Transform
analyticf = anasmile(npts,t1array,t2array,w1array,w2array);

FFTfield = fft2(field);

% Normalisation of amplitude
FFTfield = (T1o/npts)*(T2o/npts)*FFTfield;

%Shifting operations to center the Gaussian
FFTfield = fftshift(FFTfield);

for j = 0:npts-1
    for k = 0:npts-1
    FFTfield(j+1,k+1) = FFTfield(j+1,k+1)*exp(pi*1i*(j+k));
    end
end

% Obtaining the inverse
for j = 0:npts-1
    for k = 0:npts-1
    inverse2(j+1,k+1) = FFTfield(j+1,k+1)*exp(-pi*1i*(j+k));
    end
end

inverse2 = 1/(T1s*T2s)*ifft2(ifftshift(inverse2));

% Remove small amplitude contributions

tol = 1e-02;

for j = 1:npts
  for i = 1:npts
    if abs(real(FFTfield(i,j))) < tol && abs(imag(FFTfield(i,j))) < tol
        FFTfield(i,j) = 0.00 + 0.00i;
    end
    if abs(real(DFTfield(i,j))) < tol && abs(imag(DFTfield(i,j))) < tol
        DFTfield(i,j) = 0.00 + 0.00i;
    end
    if abs(real(analyticf(i,j))) < tol && abs(imag(analyticf(i,j))) < tol
        analyticf(i,j) = 0.00 + 0.00i;
    end
  end
end

% Plotting signals

% Magnitude plots for analytic vs FFT
figure(1)
pcolor(w1array,w2array,abs(FFTfield));
shading interp
%hold on
figure(2)
surf(w1array,w2array,abs(analyticf));
shading interp

% Phase plots for analytic vs FFT
%figure(2)
%surf(w1array,w2array,atan2(imag(FFTfield),real(FFTfield)));
%hold on
%surf(w1array,w2array,atan2(imag(analyticf),real(analyticf)));

% Error in phase plot (analytic vs FFT)
figure(3)
surf(w1array,w2array,unwrap((atan2(imag(FFTfield),real(FFTfield))) - (atan2(imag(analyticf),real(analyticf)))));

% Verify Inverse FFT against the input signal
figure(4)
surf(t1array,t2array,field);
hold on
surf(t1array,t2array,abs(inverse2));

% Error plot for inverse FFT vs input signal
figure(5)
surf(t1array,t2array,abs(field-inverse2));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Interpolation and coordinate conversion procedure [Method 1]: 
% No re-discretisation

% Globally interpolated Cartesian grid and assoc. Fourier transform
k1 = 1; % an N x N matrix => ((2^k-1)*(N-1) + N) x ((2^k-1)*(N-1) + N)) matrix
FFTfield_Ci = interp2(FFTfield,k1,'cubic');
w1_Ci = linspace(w10,w1end,size(w1array,2) + (size(w1array,2)-1)*(2^k1 - 1));
w2_Ci = linspace(w20,w2end,size(w2array,2) + (size(w2array,2)-1)*(2^k1 - 1));


% Transform to polar coordinates (a sampled version)
[w1_P,w2_P,FFTfield_P] = cart2pol(w1_Ci,w2_Ci,FFTfield_Ci);

% Interpolated Polar grid and assoc. Fourier transform 
%[keeps absolute positions of mesh points fixed: NO re-discretisation]
k2 = 0; % an N x N matrix => ((2^k-1)*(N-1) + N) x ((2^k-1)*(N-1) + N)) matrix
w1_Pi = linspace(w1_P(1),w1_P(size(w1_P,2)),size(w1_P,2) + (size(w1_P,2)-1)*(2^k2 - 1));
w2_Pi = linspace(w2_P(1),w2_P(size(w2_P,2)),size(w2_P,2) + (size(w2_P,2)-1)*(2^k2 - 1));
FFTfield_Pi = interp2(FFTfield_P,k2,'cubic');

% Transform polar to Cartesian
[w1arraynew,w2arraynew,FFTfield_new] = pol2cart(w1_Pi,w2_Pi,FFTfield_Pi);

FFTfield_2 = zeros(size(FFTfield,2));

countx = 0;
for j = 0:(2^(k1))*(2^(k2)):(size(w1arraynew,2)-1)
    county = 0;
    for k = 0:(2^(k1))*(2^(k2)):(size(w2arraynew,2)-1)
        FFTfield_2(countx+1,county+1) = FFTfield_new(j+1,k+1);
        county = county + 1;
    end
    countx = countx + 1;
end

inverse2_new = zeros(size(FFTfield_2,2));

% Obtaining the inverse transform after the refinements
for j = 0:size(FFTfield_2,2)-1
    for k = 0:size(FFTfield_2,2)-1
    inverse2_new(j+1,k+1) = FFTfield_2(j+1,k+1)*exp(-pi*1i*(j+k));
    end
end

t1array_new = linspace(t10,(t1end-T1s),size(FFTfield_2,2));
t2array_new = linspace(t20,(t2end-T2s),size(FFTfield_2,2));

inverse2_new = 1/(T1s*T2s)*ifft2(ifftshift(inverse2_new));

figure(7)
surf(t1array_new,t2array_new,abs(inverse2_new));
hold on
surf(t1array,t2array,abs(inverse2));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Interpolation and coordinate conversion procedure [Method 2]
% With re-discretisation in Polar Coordinates and Cartesian coordinates

% Locally interpolated Cartesian grid and assoc. Fourier transform

% Define the polar coordinate grid
w1_polar = linspace(-pi,pi,npts);
w2_polar = linspace(sqrt(2)*max([abs(w10) abs(w1end)])/(npts),sqrt(2)*max([abs(w10) abs(w1end)]),npts);

FFTfield_pol = zeros(npts);

% Fourier transforms at each point evaluated by interpolation 

% Districute npts*npts over a circle that is larger than the Cartesian
% domain
for theta = 1:npts
    for rho = 1:npts
        dx = w2_polar(rho)*cos(w1_polar(theta));
        dy = w2_polar(rho)*sin(w1_polar(theta));
        if (w10 <= dx) && (dx <= w1end) && (w20 <= dy) && (dy <= w2end)
            FFTfield_pol(theta,rho) = interp2(w1array,w2array,FFTfield,dx,dy,'cubic');
        else
            % If points fall outside the Cartesian domain, set the
            % transform to zero
            FFTfield_pol(theta,rho) = 0.000 + 0.000i;
        end
    end
end

% Singular point at r = 0, isolate it
singular = interp2(w1array,w2array,FFTfield,0.00,0.00,'cubic');

% Plot the FFT in polar space
figure(8)
surf(w2_polar,w1_polar,abs(FFTfield_pol));

FFTfield_new2 = zeros(npts);
inverse2_new2 = zeros(npts);

% Transform polar to Cartesian

% Fourier transforms at each point evaluated by interpolation 
for i = 1:npts
    for j = 1:npts
        dtheta = atan2(w2array(j),w1array(i));
        drr = sqrt((w2array(j))^2 + (w1array(i))^2);
        if dtheta == 0 && drr == 0
            FFTfield_new2(i,j) = singular;
        else
            FFTfield_new2(i,j) = interp2(w2_polar,w1_polar,FFTfield_pol,drr,dtheta,'cubic');
        end
    end
end

% Plot the FFT absolute value in Cartesian space
figure(9)
surf(w1array,w2array,abs(FFTfield_new2));

% Obtaining the inverse transform
for j = 0:npts-1
    for k = 0:npts-1
    inverse2_new2(j+1,k+1) = FFTfield_new2(j+1,k+1)*exp(-pi*1i*(j+k));
    end 
end

inverse2_new2 = 1/(T1s*T2s)*ifft2(ifftshift(inverse2_new2));

% Plot the interpolation artefacts

figure(10)
%surf(t1array,t2array,abs(inverse2_new2));
%hold on
surf(t1array,t2array,abs(inverse2_new2) - field);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Function definitions

% Multi-dimensional Gaussian
function f = gauss(x,dmu,covar)
deter = det(covar);
f = (1/((2*pi)*sqrt(deter)))*exp(-0.50*(x - dmu)'*(inv(covar)*(x - dmu)));
end 

%Continuous analytic Fourier transform of gauss(x)
function f = gauss2(w,dmu,covar) 
C = dot(w,dmu);
[V,D] = eig(covar);
L = [dot(w,V(:,1)) dot(w,V(:,2))];
lambda_1 = D(1,1);
lambda_2 = D(2,2);

f = exp(-1i*C)*exp(-0.5*(L(1)^2*lambda_1 + L(2)^2*lambda_2));
end 

% Smile function (Gaussian) by Jianhao
function f = smile(npts,t1array,t2array)

n = npts;%the number of sampled point for a single gaussian distribution
[X,Y] = meshgrid(t1array, t2array);
f = zeros(n,n);%matrix to store the sum of each part of gaussian distribution

%the outline of the face
sigma1 = [0.05 0;0 0.05]; %covariance of the outline and the mouth of face
n1 = 100; %number of total gaussian distribution for the outline
dtheta1 = 2*pi/n1;
theta1 = 0:dtheta1:2*pi-dtheta1;
r = 4.5;
mu1 = zeros(100,2);
mu1(:,1) = cos(theta1)*r;
mu1(:,2) = sin(theta1)*r;

for j = 1:n1
    ft = 1.5*mvnpdf([X(:) Y(:)],mu1(j,:),sigma1);
    f = f+reshape(ft,size(X));
end

%the mouse
n2 = 21;
dtheta2 = (11*pi/6-7*pi/6)/(n2-1);
theta = 7*pi/6:dtheta2:11*pi/6;
r=3;
mu2 = zeros(n2,2);
mu2(:,1) = cos(theta)*r;
mu2(:,2) = sin(theta)*r;
for j = 1:n2
    ft = 1.5*mvnpdf([X(:) Y(:)],mu2(j,:),sigma1); 
    f = f+reshape(ft,size(X));
end

%the eyes
sigma2 = sigma1*10;%covariance of the eyes
eyeleft = [-1.8,1];
eyeright = [1.8,1];
ft = 30*mvnpdf([X(:) Y(:)],eyeleft,sigma2); %times 30 here to make it seems darker for sigma2 here is larger than sigma1
f = f+reshape(ft,size(X));
ft = 30*mvnpdf([X(:) Y(:)],eyeright,sigma2); 
f = f+reshape(ft,size(X));
end

% Analytical smile function (Gaussian) by Jianhao
function G_2d = anasmile(npts,t1array,t2array,w1array,w2array)

n = npts; % the number of sampled point for a single gaussian distribution
[X,Y] = meshgrid(t1array,t2array);%range of x and y
[WX,WY] = meshgrid(w1array,w2array);
f = zeros(n,n);

%the outline of the face
sigma1=[0.05 0;0 0.05];%covariance of the outline and the mouse of face
n1 = 100;%number of total gaussian distribution for the outline
dtheta1 = 2*pi/n1;
theta1 = 0:dtheta1:2*pi-dtheta1;
r = 4.5;
mu1 = zeros(100,2);
mu1(:,1) = cos(theta1)*r;
mu1(:,2) = sin(theta1)*r;

for j = 1:n1
    ft = 1.5*mvnpdf([X(:) Y(:)],mu1(j,:),sigma1); 
    f = f+reshape(ft,size(X));
end

%the mouse
n2 = 21;
dtheta2 = (11*pi/6-7*pi/6)/(n2-1);
theta = 7*pi/6:dtheta2:11*pi/6;
r = 3;
mu2 = zeros(n2,2);
mu2(:,1) = cos(theta)*r;
mu2(:,2) = sin(theta)*r;
for j = 1:n2
    ft = 1.5*mvnpdf([X(:) Y(:)],mu2(j,:),sigma1); 
    f = f+reshape(ft,size(X));
end

%the eyes
sigma2=sigma1*10;%covariance of the eyes
eyeleft = [-1.8,1];
eyeright = [1.8,1];
ft=30*mvnpdf([X(:) Y(:)],eyeleft,sigma2); %times 30 here to make it seems darker for sigma2 here is larger than sigma1
f=f+reshape(ft,size(X));
ft=30*mvnpdf([X(:) Y(:)],eyeright,sigma2); 
f=f+reshape(ft,size(X));

G_2d=zeros(n,n);

for i=1:n
    for j=1:n
        WXY=[WX(i,j),WY(i,j)];
        for k=1:n1%for outline
            G_2d(i,j)=G_2d(i,j)+1.5*exp(-0.5*WXY*sigma1*WXY'-1i*WXY*mu1(k,:)');
        end
        for kk=1:n2%for mouse
            G_2d(i,j)=G_2d(i,j)+1.5*exp(-0.5*WXY*sigma1*WXY'-1i*WXY*mu2(kk,:)');
        end
        G_2d(i,j)=G_2d(i,j)+30*exp(-0.5*WXY*sigma2*WXY'-1i*WXY*eyeleft');%lefteye
        G_2d(i,j)=G_2d(i,j)+30*exp(-0.5*WXY*sigma2*WXY'-1i*WXY*eyeright');%righteye
    end 
end

end