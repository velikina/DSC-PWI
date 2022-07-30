function [res, wlut, B] = SetupMRIConvData(varargin);

% SETUPMRICONVDATA 
%
%   INPUT:
%       X0 - vector of X coordinates for sampled points
%       Y0 - vector of Y coordinates for sampled points
%       XS - vector of X coordinates for target points
%       YS - vector of Y coordinates for target points
%       DK - k-space interval used to scale window radius
%       WN - name of the convolution window; DEFAULT: KAISER
%           OPTIONS:
%                   'KAISER' - Kaiser-Bessel window
%                   'GAUSSIAN' - Gaussian window
%                   'BSPLINE' - b-spline window as per J. Velikina
%       L  - radius of the window; radius of 1 corresponds to DK; DEFAULT: 1
%       B  - filter-specific parameter
%               For 'KAISER' - parameter B (DEFAULT: calculated inside procedure)
%                   (see KAISER)
%               For 'GAUSSIAN' - parameter SGM (DEFAULT: calculated inside procedure)
%       VALS - vector of function values at (X0,Y0);
%   OUTPUT:
%       RES     - output depends on the state of VALS variable: 
%               if VALS not exist, output is convolution structure to be used in CONVMRI
%               if VALS exist, output is regridded values
%       WLUT  - convolution window LUT
%       B         -  filter-specific parameter
%
%   Author:
%       Alexei Samsonov
%       2000-2004
%       Scientific Computing and Imaging Institute
%       University of Utah
%

if  nargin<5
    error('Wrong number of arguments');
else
    [x0 y0 xs ys dk] = deal(varargin{1:5});
    x0=x0(:); y0=y0(:); xs=xs(:); ys=ys(:);
end

% default parameters
fltname = 'kaiser'; L = 1;  vals = []; os=2;
if nargin>9 & ~isempty(varargin{10}), os = varargin{10}; end
if nargin>5 & ~isempty(varargin{6}), fltname=varargin{6}; end
if nargin>6 & ~isempty(varargin{7}), L = varargin{7}; end
if nargin>7 & ~isempty(varargin{8}), 
    B = varargin{8}; sgm = B; 
else 
    B = kaiser_par(L,os); sgm = gaussian_par(L);
end
if nargin>8 & ~isempty(varargin{9}), vals = varargin{9}; end

% making LUT
wx = linspace(0, L, 10000);
switch lower(fltname)
    case 'kaiser'
        wlut = besseli(0, B*sqrt(1-(wx/L).^2))/(besseli(0,B));
    case 'gaussian'
        wlut = exp(-wx.^2/(2*sgm*sgm));
    case 'sinc'
        wlut = sinc(wx);
    case 'linear'
        wlut = wx*dk;
    case 'bspline'
        dkx=0.5*wx/L;
        wlut=bspline(dkx);
    otherwise
        error('Unknown convolution window type')
end

% binning points on bins of L size
maxvx = max([max(x0) max(xs)]); 
maxvy = max([max(y0) max(ys)]);
minvx = min([min(x0) min(xs)]); 
minvy = min([min(y0) min(ys)]);

% mapping window radius to the supplied coordinates
Lk = L*dk;

% creating array for binning with bin size corresponding to window radius
nbins = floor([(maxvx-minvx)/Lk (maxvy-minvy)/Lk])+1;

% filling in bins
ix0 = floor((x0-minvx)/Lk)+1; 
iy0 = floor((y0-minvy)/Lk)+1;
[i0 ind] = sortrows([ix0(:) iy0(:)]);
[rowv lstel] = unique(i0, 'rows','legacy');
lind = sub2ind(nbins, rowv(:,1), rowv(:,2));

% start of the range
sbins = zeros(nbins);
sbins(lind) = [1; lstel(1:end-1)+1];

% end of the range
ebins = zeros(nbins);
ebins(lind) = lstel;

is = floor((xs-minvx)/Lk+1); 
js = floor((ys-minvy)/Lk+1);

if (isempty(vals))
    [ivals, iwlut, ires] = get_conv_data_c(x0, y0, xs, ys, sbins, ebins, ind, is, js, wlut, Lk);
    res = struct('ivals', ivals, 'wlut', wlut, 'iwlut', iwlut, 'ires', ires);
else
    [res] = get_conv_data_c(x0, y0, xs, ys, sbins, ebins, ind, is, js, wlut, Lk, vals);
end
end

function beta = kaiser_par(L,os);
% Function to obtain obtimal design parameter BETA for Kaiser-Bessel window with radius L
% (from [1])
beta=pi*sqrt(4*L*L/(os*os)*(os-0.5)^2-0.8);
end

function sgm = gaussian_par(L);
% Function to obtain obtimal design parameter SGM for Gaussian window with radius L
% (from [1])
a = -0.00356904761904715;
b = 0.068877380952377;
c = 0.118093452380955;
wsz = 2*L;
sgm = a*wsz^2+b*wsz+c;
end

function vals=bspline(dkx)
K=9;
% a=1/8;
% dkx=dkx;
coef_x = (sinc(dkx)).^4.*dirichlet(dkx, K);
% coef_y = (sinc(dky)).^4.*dirichlet(dky, K);
vals=coef_x;
end