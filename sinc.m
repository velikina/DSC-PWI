function y=sinc(x)
% sinc function sinc(x)= sin(pi*x)/(pi*x)

idx=find(x==0);
x(idx)=1;
y=sin(pi*x)./(pi*x);
y(idx)=1;