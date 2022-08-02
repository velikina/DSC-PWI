function [kx, ky, phi] = getRadTr(npr, npts, kxmin, kxmax, ph0, phn);

% [kx, ky] = getRadTr(npr, npts, kxmin, kxmax, ph0, phn);
%
%   Obtain 2D radial sampling pattern
%       
%       IN:
%           npr                     - number of projections to generate
%           npts                    - number of points per projection
%           kxmin, kymax  - extent of k-space (resolution)
%           ph0, phn - starting and ending azimuthal (0 and 2*pi for full
%           coverage)
%       
%       OUT:
%           kx,ky - coordinates of the trajectory points (matrices
%           npr-by-npts)
%       
%       by Alexey Samsonov
%           2005

rr = linspace(kxmin, kxmax, npts);
cntr = ceil(npts/2);
rr = rr - rr(cntr);

if numel(ph0)==1
    phi  = linspace(ph0, phn, npr+1);
    phi  = phi(1:end-1);
else
    phi=ph0;
end

cs = cos(phi);
sn = sin(phi);
kx = []; ky = [];
for ii=1:npr
    kx = [kx rr*cs(ii)];
    ky = [ky rr*sn(ii)];
end
kx=reshape(kx,npts,npr);
ky=reshape(ky,npts,npr);
