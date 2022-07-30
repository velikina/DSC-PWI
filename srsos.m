function res=srsos(ims, nd)

%  res=srsos(ims, nd);
%
%

tt=ndims(ims);
if tt<3
    tt=3;
end
if check_var('nd')
    tt=nd;
end

res=squeeze(sqrt(sum(abs(ims).^2,tt)));