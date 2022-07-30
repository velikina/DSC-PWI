
function y = dirichlet(x, n)

% Computes Dirichlet kernel of order n
% y = sin(2*pi(4+1/2)*x)/sin(pi*x)

x_int = (x == floor(x));

y = sin(2*pi*(n+1/2)*x)./(sin(pi*x) + x_int) .* (1-x_int) + (2*n+1)*x_int;
