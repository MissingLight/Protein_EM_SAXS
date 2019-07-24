clc
clear
close all

N=128;%the number of sampled point for a single gaussian distribution
t1=-5;
t2=5;
T=(t2-t1);            % Sampling interval
F0=1/T;               % Minimum frequency interval
Fs=(N-1)*F0;          % Sampling frequency
[X,Y]=meshgrid((-N/2:N/2-1)*T/N,(-N/2:N/2-1)*T/N);
wx=(-N/2:N/2-1)*2*pi*Fs/N;
wy=wx;
[WX,WY]=meshgrid(wx,wy);
f=zeros(N,N);

%the outline of the face
sigma1=[0.05 0;0 0.05];%covariance of the outline and the mouse of face
n1=100;%number of total gaussian distribution for the outline
dtheta1=2*pi/n1;
theta1=0:dtheta1:2*pi-dtheta1;
r=4.5;
mu1=zeros(100,2);
mu1(:,1)=cos(theta1)*r;
mu1(:,2)=sin(theta1)*r;
for j=1:n1
    ft=1.5*mvnpdf([X(:) Y(:)],mu1(j,:),sigma1); 
    f=f+reshape(ft,size(X));
end

%the mouth
n2=21; 
mouth1=7*pi/6; mouth2=11*pi/6;
dtheta2=(mouth2-mouth1)/(n2-1);
theta=mouth1:dtheta2:mouth2;
r=3;
mu2=zeros(n2,2);
mu2(:,1)=cos(theta)*r;
mu2(:,2)=sin(theta)*r;
for j=1:n2
    ft=1.5*mvnpdf([X(:) Y(:)],mu2(j,:),sigma1); 
    f=f+reshape(ft,size(X));
end

%the eyes
sigma2=sigma1*10;%covariance of the eyes
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

G_ana=zeros(N,N);
for i=1:N
    for j=1:N
        WXY=[WX(i,j),WY(i,j)];
        for k=1:n1%for outline
            G_ana(i,j)=G_ana(i,j)+1.5*exp(-0.5*WXY*sigma1*WXY'-1i*WXY*mu1(k,:)');
        end
        for kk=1:n2%for mouth
            G_ana(i,j)=G_ana(i,j)+1.5*exp(-0.5*WXY*sigma1*WXY'-1i*WXY*mu2(kk,:)');
        end
            G_ana(i,j)=G_ana(i,j)+30*exp(-0.5*WXY*sigma2*WXY'-1i*WXY*eyeleft');%lefteye
            G_ana(i,j)=G_ana(i,j)+30*exp(-0.5*WXY*sigma2*WXY'-1i*WXY*eyeright');%righteye
        for kk=1:n3%for mouth
            G_ana(i,j)=G_ana(i,j)+1.2*exp(-0.5*WXY*sigma3*WXY'-1i*WXY*mu3(kk,:)');
            G_ana(i,j)=G_ana(i,j)+1.2*exp(-0.5*WXY*sigma3*WXY'-1i*WXY*mu4(kk,:)');
        end
    end 
end

G_ana_abs=abs(G_ana);
figure(1);
surf(WX,WY,G_ana_abs);
shading interp
xlabel('w1');ylabel('w2');zlabel('Magnitude');
title('Sample of Continuous FT(Magnitude)');

G_ana_real=real(G_ana);
G_ana_imag=imag(G_ana);
index=find(abs(G_ana_real)<1e-2&abs(G_ana_imag)<1e-2);
G_ana_real(index)=0;
G_ana_imag(index)=0;

G_ana_phase=atan2(G_ana_real,G_ana_imag);
figure(2);
mesh(WX,WY,G_ana_phase);view(2);
%shading interp
xlabel('w1');ylabel('w2');zlabel('Phase');
title('Sample of Continuous FT(Phase)');
colorbar

figure(3)
scatter3(reshape(WX,[N*N,1]),reshape(WY,[N*N,1]),reshape(G_ana_phase,[N*N,1]))