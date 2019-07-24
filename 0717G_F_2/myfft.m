function y=myfft(x)
N=length(x);
if N==1
    y=x;
else
    u=x(1:2:end);
    v=x(2:2:end);

    yu=myfft(u);
    yv=myfft(v);
    w=exp(-2*pi*1i/N);
    
    y=zeros(N,1);
    y(1:N/2)=yu(1:N/2)+yv(1:N/2).*((w.^(0:N/2-1)).');
    y(N/2+1:N)=yu(1:N/2)-yv(1:N/2).*((w.^(0:N/2-1)).');
end
end