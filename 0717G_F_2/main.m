%% 1
close all
clear
clc

mu=0;  sigma=1;
N=128;                %sample number
t1=-5*sigma;t2=5*sigma;
T=(t2-t1);            % Sampling interval
deltaT=T/(N-1);
F0=1/T;               % Minimum frequency interval
Fs=1/deltaT;          % Sampling frequency

%t=linspace(t1,t2,N)'; % time series
t=((-N/2:N/2-1)*T/N)';
g=gauss(t,mu,sigma);  % sampled values of the  Gaussian distribution 
xx=t1:0.01:t2;
yy=gauss(xx,mu,sigma);
figure(1);
plot(xx,yy,'b-');
hold on
stem(t,g)
xlabel('t(s)');ylabel('Magnitude')
title('Gaussian Function and Sampling');
legend('Gaussina Function','Sampling'); 

w=2*pi*Fs*(0:N-1)/N;
wshift=(-N/2:N/2-1)*2*pi*Fs/N;
G=Fouriergauss(wshift,mu,sigma);  % sampled values of the continuous Fourier transform
G_fftcomp=fft(g);
G_fftcomp_shift=fftshift(fft(g));
G_fft=abs(fft(g));
%G_fftshift=abs(fftshift(fft(ifftshift(g))))/Fs;
G_fftshift=fftshift(abs(fft(g)))*T/N;
%G_DFTgauss=fftshift(abs(DFT(g)))/Fs;
%G_myfft=fftshift(abs(myfft(g)))/Fs;

figure(2);
set(gcf,'Position',get(gcf,'Position').*[1 1 2 1])

%1
% subplot(2,2,1);
% plot(w,G_fft,'o');
% xlabel('w');ylabel('Amplitude');
% title('FFT N=128')
% grid on;
% %
% subplot(2,2,2);
% plot(w(1:N/2),G_fft(1:N/2),'o');
% xlabel('w');ylabel('Amplitude');
% title('FFT Before Nyquist Frequency')
% grid on;
%
subplot(1,4,[1 2]);
plot(wshift,G,'o-');
xlabel('w');ylabel('Magnitude');
title('Sample of the Continuous Fourier Transform')
grid on;
%
subplot(1,4,[3 4]);
plot(wshift',G_fftshift,'o');
xlabel('w');ylabel('Magnitude');
title('FFTSHIFT and Mdified Amplitude')
grid on;

figure(3)
wleft=t1/sigma^2;wright=t2/sigma^2;
index=find((wshift>wleft)&(wshift<wright));
w_continuous=wleft:0.001:wright;
G_continuous=Fouriergauss(w_continuous,mu,sigma);
plot(w_continuous,G_continuous);
hold on 
plot(wshift(index),G(index),'o','markersize',10);
hold on
plot(wshift(index),G_fftshift(index),'*','color','r');
xlabel('w(Angular frequency)');
ylabel('Magnitude');
legend('Fourier of Gaussian Function','Sample of Continuous FT ','FFT')
title('Comparisons of using FFT')
grid on 

figure(4)
G_error=abs(G_fftshift'-G);
plot(wshift,G_error,'o-');
xlabel('w');ylabel('Error Magnitude')
grid on


% figure(5)
% for k=0:N-1
%     sumcos=0;sumsin=0;    
%     for n=0:N-1
%         sumcos=sumcos+cos(2*pi*k/N*n)*g(n+1);
%         sumsin=sumsin-sin(2*pi*k/N*n)*g(n+1);
%     end
%    G_real(k+1)=sumcos;
%    G_imag(k+1)=sumsin;
% end
% G_expansion=fftshift(sqrt(G_real.^2+G_imag.^2))/Fs;
% G_error2=abs(G_expansion-G);
% plot(wshift,G_error2);

 b=sum(1/sqrt(2*pi*sigma^2)*exp(-t.^2/(2*sigma^2)))*T/N;
abs(1-b)

%%
close all
clear
clc

mu=0;  sigma=1;
N=128;                %sample number
t1=-5*sigma;t2=5*sigma;
T=(t2-t1);            % Sampling interval
deltaT=T/(N-1);
F0=1/T;               % Minimum frequency interval
Fs=1/deltaT;          % Sampling frequency

%t=linspace(t1,t2,N)'; % time series
t=((-N/2:N/2-1)*T/N)';
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
G_fftcomp=fft(g);
G_fftcomp_shift=fftshift(fft(g));
G_fft=abs(fft(g));
%G_fftshift=abs(fftshift(fft(ifftshift(g))))/Fs;
G_fftshift=fftshift(abs(fft(g)))*T/N;

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
wleft=t1/sigma^2;wright=t2/sigma^2;
index=find((wshift>wleft)&(wshift<wright));
w_continuous=wleft:0.1:wright;
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
G_real=real(G_fftcomp_shift);
G_imag=imag(G_fftcomp_shift);
G_phase=atan2(abs(G_imag),abs(G_real));
plot(wshift,G_phase,'o');
grid on
axis([-inf inf -4 4])

figure(5)
G_phase2=angle(abs(G_real));
%G_phase2=unwrap(G_phase2);
plot(wshift,G_phase2,'o');
grid on

figure(6)
G_phase_ana=angle(G);
plot(w,G_phase_ana,'o');
grid on


for k=0:N-1
    sumcos=0;sumsin=0;
    for n=0:N-1
        sumcos=sumcos+cos(2*pi*k/N*n)*g(n+1);
        sumsin=sumsin-sin(2*pi*k/N*n)*g(n+1);
    end
    G_realcos(k+1)=sumcos;
    G_imagsin(k+1)=sumsin;
    aaangle(k+1)=atan2(sumsin,sumcos);
end
figure(7)
plot(wshift,(aaangle),'o')
axis([-inf inf -4 4])
grid on


%% 

clc
clear
close all

N=128;
mu=[0.5,-0.5];
sigma=[1 -0.5; -0.5 1];
t1=-5;t2=5;
T=(t2-t1);            % Sampling interval
F0=1/T;               % Minimum frequency interval
Fs=(N-1)*F0;          % Sampling frequency
wx=(-N/2:N/2-1)*2*pi*Fs/N;
wy=wx;
[WX,WY]=meshgrid(wx,wy);
[X,Y]=meshgrid(-5:10/127:5,-5:10/127:5);
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
        G_2d_phase_ana(i,j)=angle(G_2d(i,j));
        G_2d(i,j)=abs(G_2d(i,j)); 
    end 
end

G_2dfft_comp=fftshift(fft2(p))/Fs^2;
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
 surf(WX,WY,G_2dfftshift);
 shading interp
% for i=index(1):index(end)
%     for j=index(1):index(end)
%         plot3(WX(i,j),WY(i,j),G_2dfftshift(i,j),'.','markersize',10);
%         hold on
%     end
% end
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

% figure(4)
% G_2phase=angle(G_2dfft_comp);
% for i=1:N 
%     for j=1:N
%        plot3(WX(i,j),WY(i,j),G_2phase(i,j),'.','markersize',10);
%        hold on; 
%     end 
% end

figure(5)
G_2phase=angle(G_2dfft_comp);
surf(WX,WY,G_2phase),view(2);
shading interp
colorbar

figure(6)
surf(WX,WY,G_2d_phase_ana),view(2);
shading interp
colorbar

