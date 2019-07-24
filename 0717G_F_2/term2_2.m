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
colorbar
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
        G_real2=real(G_2d(i,j));
        G_imag2=imag(G_2d(i,j));
        if abs(G_real2)<1e-2&&abs(G_imag2)<1e-2
        G_real2=0;
        G_imag2=0;
        end
        G_2d_phase_ana(i,j)=atan2(G_imag2,G_real2);
        G_2d(i,j)=abs(G_2d(i,j)); 
    end 
end

G_2dfft_comp=fftshift(fft2(p))*T^2/N^2;
G_2dfft=abs(fft2(p))*T^2/N^2;
G_2dfftshift=fftshift(G_2dfft);
for i=1:N 
    for j=1:N
        G_2dfft_comp_m(i,j)=G_2dfft_comp(i,j)*exp(-1i*pi*(i-1)-1i*pi*(j-1));
    end 
end

% figure(4)
% G_2phase=angle(G_2dfft_comp);
% for i=1:N 
%     for j=1:N
%        plot3(WX(i,j),WY(i,j),G_2phase(i,j),'.','markersize',10);
%        hold on; 
%     end 
% end

figure(2)
G_real=real(G_2dfft_comp_m);
G_imag=imag(G_2dfft_comp_m);
index=find(abs(G_real)<1e-2&abs(G_imag)<1e-2);
G_real(index)=0;
G_imag(index)=0;

G_2phase=atan2(G_imag,G_real);
surf(WX,WY,G_2phase),view(2);
shading interp
colorbar
title('FFT Phase')

figure(3)
surf(WX,WY,G_2d_phase_ana),view(2);
shading interp
colorbar
title('Analytical Phase')

figure(4)
surf(WX,WY,G_2d_phase_ana-G_2phase),view(2);
shading interp
colorbar
title('Phase Error')

%% IFFT

ggg=abs(ifftshift(ifft2(G_2dfft_comp_m)))*N^2/T^2;
surf(X,Y,abs(ggg))
shading interp
colorbar