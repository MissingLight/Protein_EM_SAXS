function f=gauss(x,mu,sigma)
%x: n-1 vector
%mu:n-1 vector mean; 
%sigma: n-n matrix variance;
%n=length(x);
%f=1/((2*pi)^(n/2)*sqrt(det(sigma)))*exp(-0.5*(x-mu)'*sigma^(-1)*(x-mu));
n=1;
f=1/((2*pi)^(n/2)*(det(sigma)))*exp(-0.5*(x-mu).^2*sigma^(-2));
end