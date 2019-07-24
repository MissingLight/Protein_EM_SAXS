function y=DFT(x)
N=length(x);
omega=exp(-2*pi*1i/N);
j=0:N-1;
DFTmatrix=omega.^(j'*j);
y=DFTmatrix*x;
end