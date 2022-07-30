function res=vect(data,pp)
%
% res=vect(data,pp)
%
if ~exist('pp')
    res=data(:);
else
    res=data(pp);
end
end