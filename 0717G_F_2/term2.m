%% Gaussian mu is not 0

close all
clear
clc

mu=1;  sigma=1;
N=128;                %sample number
t1=(-5)*sigma;t2=(5)*sigma;
T=(t2-t1);            % Sampling interval
deltaT=T/(N-1);
F0=1/T;               % Minimum frequency interval
Fs=1/deltaT;          % Sampling frequency

%t=linspace(t1,t2,N)'; % time series
t=((-N/2:N/2-1)*T/N)';
% index_t=find(t>5);
% aaa=t(index_t)-10;
% t(index_t)=[];
% t=[aaa;t];

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
wshift=((-N/2:N/2-1)*2*pi*Fs/N)';
G_anacomp=Fouriergauss(wshift,mu,sigma);  % sampled values of the continuous Fourier transform
G_ana=abs(G_anacomp);
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
plot(wshift,G_ana,'o');
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
plot(w_continuous,abs(G_continuous));
hold on 
plot(wshift(index),G_ana(index),'o','markersize',10);
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
index=find(abs(G_real)<5e-2&abs(G_imag)<5e-2);
G_real(index)=0;
G_imag(index)=0;

for k=1:2:length(G_imag)
    G_real(k+1)=-G_real(k+1);
    G_imag(k+1)=-G_imag(k+1);
end
G_phase=atan2(G_imag,G_real);
plot(wshift,G_phase,'o-');
grid on
title('Phase of FFT')

% figure(5)
% G_phase2=angle(abs(G_real));
% %G_phase2=unwrap(G_phase2);
% plot(wshift,G_phase2,'o');
% grid on

figure(6)
G_real2=real(G_anacomp);
G_imag2=imag(G_anacomp);
index=find(abs(G_real2)<1e-3&abs(G_imag2)<1e-3);
G_real2(index)=0;
G_imag2(index)=0;
G_phase_ana=atan2(G_imag2,G_real2);
plot(wshift,G_phase_ana,'o-');
grid on
title('Phase of Analytical Fourier')

figure(7)
plot(wshift,G_phase_ana-G_phase,'.-')
hold on
plot(wshift,G_phase_ana,'.-','LineWidth',2);
hold on
plot(wshift,G_phase,'o');
legend('Error','Analytical','FFT')
xlabel('w');ylabel('phase');


for k=0:N-1
    sumcos=0;sumsin=0;
    for n=0:N-1
        sumcos=sumcos+cos(2*pi*k/N*n)*g(n+1);
        sumsin=sumsin-sin(2*pi*k/N*n)*g(n+1);
    end
    if abs(sumcos)<1e-2&&abs(sumsin)<1e-2
        sumcos=0;
        sumsin=0;
    end
    G_realcos(k+1)=sumcos;
    G_imagsin(k+1)=sumsin;
    aaangle(k+1)=atan2(sumsin,sumcos);
end
figure(8)
aaangle=fftshift(aaangle);
plot(wshift,(aaangle),'o-')
grid on



figure(9)
series=(-N/2:N/2-1)';
G_fftcomp_m=G_fftcomp_shift.*exp(-1i*pi*series);
G_real_m=real(G_fftcomp_m);
G_imag_m=imag(G_fftcomp_m);
index=find(abs(G_real_m)<1e-4&abs(G_imag_m)<1e-4);
G_real_m(index)=0;
G_imag_m(index)=0;

G_phase_m=atan2(G_imag_m,G_real_m);
plot(wshift,G_phase_m,'o-');
grid on
title('Phase of FFT')

figure(10)
plot(wshift,G_phase_ana-G_phase,'.-')
title('Error of Phase between FFT and Analitical')