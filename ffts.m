function [res]=ffts(data,dims,dr);

%   [res]=ffts(data,dims,dir);
%
%       Multipurpose FFT function
%               data -  ND array of data
%               dims - vector of dimensions to apply FFT along
%               dir    - direction of FFT ("1" forward, "-1" inverse)
%
%       by Alexey Samsonov
%       U of W, Madison
%       2005

res=data;
if ~exist('dims') | isempty(dims)
    dims=[1 2];
end

if ~exist('dr') | isempty(dr)
    dr=1;
end

for ii=1:length(dims)
    if dr==1
        res=fft(res,[],dims(ii));
    elseif dr==-1
        res=ifft(res,[],dims(ii));
    end
end