close all
clear
clc

mu=0;  sigma=1;
N=128;                %sample number
t1=-5;t2=5;
T=(t2-t1);            % Sampling interval
F0=1/T;               % Minimum frequency interval
Fs=(N-1)*F0;          % Sampling frequency

t=[t1:(t2-t1)/(N-1):t2]'; % time series
g=gauss(t,mu,sigma);  % sampled values of the  Gaussian distribution 
xx=t1:0.01:t2;
yy=gauss(xx,mu,sigma);
figure(1);
plot(xx,yy,'b-');
hold on
stem(t,g)
xlabel('t(s)');ylabel('Amplitude')
title('Gaussian Function and Sampling');
legend('Gaussina Function','Sampling'); 

w=2*pi*Fs*(0:N-1)/N;
wshift=(-N/2:N/2-1)*2*pi*Fs/N;
G=Fouriergauss(wshift,mu,sigma);  % sampled values of the continuous Fourier transform
G_fft=abs(fft(g));
G_fftshift=fftshift(abs(fft(g)))/Fs;
G_DFTgauss=fftshift(abs(DFT(g)))/Fs;
G_myfft=fftshift(abs(myfft(g)))/Fs;

figure(2);
%1
subplot(2,2,1);
plot(w,G_fft,'o');
xlabel('w');ylabel('Amplitude');
title('FFT N=128')
grid on;
%
subplot(2,2,2);
plot(w(1:N/2),G_fft(1:N/2),'o');
xlabel('w');ylabel('Amplitude');
title('FFT Before Nyquist Frequency')
grid on;
%
subplot(2,2,3);
plot(wshift,G,'o');
xlabel('w');ylabel('Amplitude');
title('Sample of the Continuous Fourier Transform')
grid on;
%
subplot(2,2,4);
plot(wshift',G_fftshift,'o');
xlabel('w');ylabel('Amplitude');
title('FFTSHIFT and Mdified Amplitude')
grid on;

figure(3)
index=find((wshift>-5)&(wshift<5));
w_continuous=t1:0.01:t2;
G_continuous=Fouriergauss(w_continuous,mu,sigma);
plot(w_continuous,G_continuous);
hold on 
plot(wshift(index),G(index),'o','markersize',10);
hold on
plot(wshift(index),G_fftshift(index),'*','color','r');
xlabel('w(Angular frequency)');
ylabel('Amplitude');
legend('Fourier of Gaussian Function','Sample of Continuous FT ','FFT')
title('Comparisons of using FFT')
grid on 

figure(4)
index=find((wshift>-5)&(wshift<5));
w_continuous=t1:0.01:t2;
G_continuous=Fouriergauss(w_continuous,mu,sigma);
plot(w_continuous,G_continuous);
hold on 
plot(wshift(index),G(index),'o','markersize',10);
hold on
plot(wshift(index),G_DFTgauss(index),'+','color','r');
xlabel('w(Angular frequency)');
ylabel('Amplitude');
legend('Fourier of Gaussian Function','Sample of Continuous FT','DFT')
title('Comparisons of using DFT')
grid on 

figure(5)
index=find((wshift>-5)&(wshift<5));
w_continuous=t1:0.01:t2;
G_continuous=Fouriergauss(w_continuous,mu,sigma);
plot(w_continuous,G_continuous);
hold on 
plot(wshift(index),G(index),'o','markersize',10);
hold on
plot(wshift(index),G_DFTgauss(index),'.','markersize',10,'color','r');
xlabel('w(Angular frequency)');
ylabel('Amplitude');
legend('Fourier of Gaussian Function','Sample of Continuous FT','My FTT')
title('Comparisons of using My FFT')
grid on 

tic
Gt_fft=fft(g)/Fs;
t_fft=toc
tic
Gt_DFT=DFT(g)/Fs;
t_dft=toc
tic
Gt_myfft=myfft(g)/Fs;
t_myfft=toc

%%

clc
clear
close all

N=128;
mu=[4.5,0];
sigma=[0.05 0; 0 0.05];
t1=-5;t2=5;
T=(t2-t1);            % Sampling interval
F0=1/T;               % Minimum frequency interval
Fs=(N-1)*F0;          % Sampling frequency
wx=(-N/2:N/2-1)*2*pi*Fs/N;
wy=wx;
[WX,WY]=meshgrid(wx,wy);
[X,Y]=meshgrid((-N/2:N/2-1)*T/N,(-N/2:N/2-1)*T/N);
p=mvnpdf([X(:) Y(:)],mu,sigma);  %the joint probability density, the Z-axis.
p=reshape(p,size(X));            
figure(1)
set(gcf,'Position',get(gcf,'Position').*[1 1 1.3 1])
subplot(2,3,[1 2 4 5])
mesh(X,Y,p),axis tight,title('2-D Gaussian Distribution')
xlabel('x1');ylabel('x2');
subplot(2,3,3)
mesh(X,Y,p),view(2),axis tight,title('Projection on XOY Plane')
xlabel('x1');ylabel('x2');zlabel('Amplitude');
subplot(2,3,6)
mesh(X,Y,p),view([0 0]),axis tight,title('Projection on XOZ Plane')
xlabel('x1');ylabel('z');

for i=1:N 
    for j=1:N
        WXY=[WX(i,j),WY(i,j)];
        G_2d(i,j)=exp(-0.5*WXY*sigma*WXY'-1i*WXY*mu');
        G_2d(i,j)=abs(G_2d(i,j)); 
    end 
end

G_2dfft=abs(fft2(p))/Fs^2;
G_2dfftshift=fftshift(G_2dfft);

figure(2)
subplot(2,2,1);
mesh(WX,WY,G_2d);
xlabel('w1');ylabel('w2');zlabel('Amplitude');
title('Sample of Continuous FT')

subplot(2,2,2);
%imagesc(G_2d);
surf(WX,WY,G_2d),view(2),axis tight,title('Top View of Sample of Continuous FT')
shading interp
xlabel('w1');ylabel('w2');title('Top View of Sample of Continuous FT')
colorbar

subplot(2,2,3);
mesh(WX,WY,G_2dfft);
xlabel('w1');ylabel('w2');zlabel('Amplitude');
title('FFT of the 2-D Gaussian Distribution')

subplot(2,2,4);
surf(WX,WY,G_2dfft),view(2),title('Top View')
shading interp
%imagesc(G_2dfft);
xlabel('w1');ylabel('w2');title('Top View of FFT')
colorbar

figure(3)
index=find((wx>-5)&(wx<5));
set(gcf,'Position',get(gcf,'Position').*[1 1 1.3 1])
subplot(2,3,[1 2 4 5])
% surf(WX,WY,G_2dfftshift);
% shading interp
for i=index(1):index(end)
    for j=index(1):index(end)
        plot3(WX(i,j),WY(i,j),G_2dfftshift(i,j),'.','markersize',10);
        hold on
    end
end
 xlabel('w1');ylabel('w2');title('Shift FFT of 2-D Guassian')

subplot(2,3,3)
w_continuous=t1:0.02:t2;
[W_C1,W_C2]=meshgrid(w_continuous,w_continuous);
for i=1:length(w_continuous) 
    for j=1:length(w_continuous)
        W_C=[W_C1(i,j),W_C2(i,j)];
        G_2d_C(i,j)=exp(-0.5*W_C*sigma*W_C'-1i*W_C*mu');
        G_2d_C(i,j)=abs(G_2d_C(i,j)); 
    end 
end
surf(W_C1,W_C2,G_2d_C),view(2);
xlabel('w1');ylabel('w2');title('Top View of Continuous 2-D Fourier')
shading interp
colorbar

subplot(2,3,6)
surf(WX(index,index),WY(index,index),G_2dfftshift(index,index)),view(2);
xlabel('w1');ylabel('w2');title('Top View of 2-D Shift FFT');
shading interp
colorbar