function [vals] = convmridata(vals0,cnvstr,dns,dimsn);

% CONVMRIDATA
%   Performs fast convolution of MRI data to new sampling points
%
%   INPUT:
%       VALS0 - k-space values at the sampled positions
%       CNVSTR - convolution structure produced by SETUPMRICONVDATA, if no VALS variable is supplied
%       dns - sampling density function (OPTIONAL)
%       dimsn - dimensions of Cartesian grid to reshape on (OPTIONAL)
%
%   OUTPUT:
%      VALS - regridded values
%
%   Author:
%       Alexey Samsonov
%       2000-2004


dims=size(vals0);
if (~exist('dns')|isempty(dns))
    dns=1;
end

if length(dims)>1
    if length(dims)==3
        vals0=reshape(vals0,[dims(1)*dims(2) dims(3)]);
        dims=size(vals0);
    end
    vv=zeros(numel(cnvstr.ires),dims(2));
    for ii=1:dims(2)
        if dims(1)<numel(dns)
            dns0=dns([1:dims(1)]+(ii-1)*dims(1));
        else
            dns0=dns;
        end
        vv(:,ii) = conv_c(vals0(:,ii).*dns0, cnvstr.ivals, cnvstr.wlut, cnvstr.iwlut, cnvstr.ires);
    end
else
    vv = conv_c(vals0.*dns, cnvstr.ivals, cnvstr.wlut, cnvstr.iwlut, cnvstr.ires);
end

if (~exist('dimsn')|isempty(dimsn))
    vals=vv;
else
    vals=reshape(vv,dimsn);
end
