function f_hat=Fouriergauss(w,mu,sigma)
%w: 
%mu:mean; 
%sigma: variance;
%n=length(w);
f_hat=exp(-0.5*w.^2*sigma-1i*w*mu);
end