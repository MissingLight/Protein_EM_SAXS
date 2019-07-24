clc
clear
close all

N=128;         %the number of sampled point for a single gaussian distribution
t1=-5;t2=5;
T=(t2-t1);            % Sampling interval
F0=1/T;               % Minimum frequency interval
Fs=(N-1)*F0;          % Sampling frequency
[X,Y]=meshgrid((-N/2:N/2-1)*T/N,(-N/2:N/2-1)*T/N);
wx=(-N/2:N/2-1)*2*pi*Fs/N;
wy=wx;
[WX,WY]=meshgrid(wx,wy);
%d=(t2-t1)/(n-1);
%[X,Y]=meshgrid(t1:d:t2,t1:d:t2);%range of x and y
f=zeros(N,N);%matrix to store the sum of each part of gaussian distribution

%the outline of the face
sigma1=[0.05 0;0 0.05];%covariance of the outline and the mouse of face
n1=100;%number of total gaussian distribution for the outline
dtheta1=2*pi/n1;
theta1=0:dtheta1:2*pi-dtheta1;
r1=4.5;
mu1=zeros(n1,2);
mu1(:,1)=cos(theta1)*r1;
mu1(:,2)=sin(theta1)*r1;
for j=1:n1
    ft=1.5*mvnpdf([X(:) Y(:)],mu1(j,:),sigma1); 
    f=f+reshape(ft,size(X));
end

%the mouth
n2=21; 
mouth1=7*pi/6; mouth2=11*pi/6;
dtheta2=(mouth2-mouth1)/(n2-1);
theta=mouth1:dtheta2:mouth2;
r2=3;
mu2=zeros(n2,2);
mu2(:,1)=cos(theta)*r2;
mu2(:,2)=sin(theta)*r2;
for j=1:n2
    ft=1.5*mvnpdf([X(:) Y(:)],mu2(j,:),sigma1); 
    f=f+reshape(ft,size(X));
end

%the eyes
sigma2=sigma1*10; %covariance of the eyes
eyeleft=[-1.8,1];
eyeright=[1.8,1];
ft=30*mvnpdf([X(:) Y(:)],eyeleft,sigma2); %times 30 here to make it seems darker for sigma2 here is larger than sigma1
f=f+reshape(ft,size(X));
ft=30*mvnpdf([X(:) Y(:)],eyeright,sigma2); 
f=f+reshape(ft,size(X));

%eye outline
sigma3=sigma1*0.8;
n3=20;
mu3=zeros(n3,2);
mu4=mu3;
dtheta3=2*pi/n3;
theta3=0:dtheta3:2*pi-dtheta3;
r3=1.2;
mu3(:,1)=eyeleft(1)+cos(theta3)*r3;
mu3(:,2)=eyeleft(2)+sin(theta3)*r3;
mu4(:,1)=eyeright(1)+cos(theta3)*r3;
mu4(:,2)=eyeright(2)+sin(theta3)*r3;
for j=1:n3
    ft1=1.2*mvnpdf([X(:) Y(:)],mu3(j,:),sigma3); 
    ft2=1.2*mvnpdf([X(:) Y(:)],mu4(j,:),sigma3);
    f=f+reshape(ft1,size(X))+reshape(ft2,size(X));
end

%% analytical Fourier

G_ana_comp=zeros(N,N);
for i=1:N
    for j=1:N
        WXY=[WX(i,j),WY(i,j)];           
        G_ana_comp(i,j)=sum(1.5*exp(-0.5*WXY*sigma1*WXY'-1i*WXY*(mu1')))... % face outline       
                       +sum(1.5*exp(-0.5*WXY*sigma1*WXY'-1i*WXY*(mu2')))...% mouth
                       +30*exp(-0.5*WXY*sigma2*WXY'-1i*WXY*(eyeleft'))...%lefteye
                       +30*exp(-0.5*WXY*sigma2*WXY'-1i*WXY*(eyeright'))...%righteye
                       +sum(1.2*exp(-0.5*WXY*sigma3*WXY'-1i*WXY*(mu3')))...%for left eye outline
                       +sum(1.2*exp(-0.5*WXY*sigma3*WXY'-1i*WXY*(mu4')));%for right eye outline  
    end 
end

G_ana_abs=abs(G_ana_comp);
figure(1);
surf(WX,WY,G_ana_abs);
axis([-5 5 -5 5 0 300])
shading interp
xlabel('w1');ylabel('w2');zlabel('Magnitude');
title('Analytical FT - Magnitude');

%% FFT

G_fft_comp=fftshift(fft2(f))*T^2/N^2;
for i=1:N 
    for j=1:N
        G_fft_comp_m(i,j)=G_fft_comp(i,j)*exp(-1i*pi*(i-1)-1i*pi*(j-1)); %modify the phase
    end 
end
G_fft_abs=abs(G_fft_comp_m);
figure(2)
surf(WX,WY,G_fft_abs)
shading interp
xlabel('w1');ylabel('w2');zlabel('Magnitude');
title('FFT - Magnitude')

figure(3)
error_mag=G_ana_abs-G_fft_abs;
surf(WX,WY,abs(error_mag))
shading interp
xlabel('w1');ylabel('w2');zlabel('Error');
title('Error between Analytical and FFT - Magnitude')

%% Analytical phase

G_ana_real=real(G_ana_comp);
G_ana_imag=imag(G_ana_comp);
index=find(abs(G_ana_real)<1e-2&abs(G_ana_imag)<1e-2);
% index1=find(abs(G_ana_real)<1e-2);
% index2=find(abs(G_ana_imag)<1e-2);
G_ana_real(index)=0;
G_ana_imag(index)=0;
G_ana_comp_m=G_ana_real+1i*G_ana_imag;

G_ana_phase=atan2(G_ana_real,G_ana_imag);
figure(4);
mesh(WX,WY,G_ana_phase);view(2);
%shading interp
xlabel('w1');ylabel('w2');zlabel('Phase');
title('Analytical FT - Phase');
colorbar

%% FFT phase

G_fft_real=real(G_fft_comp_m);
G_fft_imag=imag(G_fft_comp_m);
% index1=find(abs(G_fft_real)<1e-2);
% index2=find(abs(G_fft_imag)<1e-2);
index=find(abs(G_fft_real)<1e-2&abs(G_fft_imag)<1e-2);
G_fft_real(index)=0;
G_fft_imag(index)=0;
G_fft_comp_m2=G_fft_real+1i*G_fft_imag;

G_fft_phase=atan2(G_fft_real,G_fft_imag);
figure(5);
mesh(WX,WY,G_fft_phase);view(2);
%shading interp
xlabel('w1');ylabel('w2');zlabel('Phase');
title('FFT - Phase');
colorbar

figure(6)
error_phase=G_ana_phase-G_fft_phase;
surf(WX,WY,error_phase)
shading interp
xlabel('w1');ylabel('w2');zlabel('Error');
title('Error between Analytical and FFT - Phase')
colorbar

%% IFFT

figure (7);
surf(X,Y,f);
shading interp
%axis ([-10 10 -10 10 0 10]);
xlabel('x');ylabel('y');zlabel('Smiley Gaussian');
title('Smiley Gaussian Distributions - Analytical');
%set(h,'edgecolor','none');
colorbar


figure(8)
f_ifft=abs(ifftshift(ifft2(G_fft_comp_m)))*N^2/T^2;
surf(X,Y,f_ifft)
shading interp
colorbar
xlabel('x');ylabel('y');zlabel('Smiley Gaussian');
title('Smiley Gaussian Distributions - IFFT');

figure(9)
surf(X,Y,f_ifft-f)
shading interp
colorbar
xlabel('x');ylabel('y');zlabel('Smiley Gaussian');
title('Error between IFFT and True Gaussian');

%%

% [theta_W,rho_W,G_pol]=cart2pol(WX,WY,G_fft_abs);
% figure