function res=spreadData(data, obmask);

% Distributes data from each column of 'data' to each slice of 3D array 'res' using 2D mask array 'obmask'
%   Returns array size(obmask)-by-nch, where nch-number of channels in
%   'data' (second dimension)
%
%  AS, 2006

nch=size(data,2);
res=zeros([size(obmask) nch]);
res(repmat(logical(obmask), [1 1 nch]))=data(:);