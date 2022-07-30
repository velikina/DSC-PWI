function y=fftshifts(x,dims);

y=x;
if ~exist('dims') | isempty(dims)
    dims=[1 2];
end
for ii=1:length(dims)
    y=fftshift(y,dims(ii));
end