clc
clear
close all

N=128;         %the number of sampled point for a single gaussian distribution
t1=-5;t2=5;
T=(t2-t1);            % Sampling interval
F0=1/T;               % Minimum frequency interval
Fs=(N-1)*F0;          % Sampling frequency
[X,Y]=meshgrid((-N/2:N/2-1)*T/N,(-N/2:N/2-1)*T/N);
%d=(t2-t1)/(n-1);
%[X,Y]=meshgrid(t1:d:t2,t1:d:t2);%range of x and y
f=zeros(N,N);%matrix to store the sum of each part of gaussian distribution

%the outline of the face
sigma1=[0.05 0;0 0.05];%covariance of the outline and the mouse of face
n1=100;%number of total gaussian distribution for the outline
dtheta1=2*pi/n1;
theta1=0:dtheta1:2*pi-dtheta1;
r1=4.5;
mu1=zeros(100,2);
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


%plot
figure (1);
set(gcf,'Position',get(gcf,'Position').*[1 1 1.3 1]);
subplot(2,3,[1 2 4 5]);
h=surf(X,Y,f);
axis ([-10 10 -10 10 0 10]);
title('The figure of the sum of gaussian distributions');
set(h,'edgecolor','none');

subplot(2,3,3);
h2=surf(X,Y,f);view(2);
axis ([-5 5 -5 5]);
title('Projection on plane XOY');
set(h2,'edgecolor','none');
