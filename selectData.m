function res=selectData(data, obmask);

% Select data from each slice of 3D array 'data' using 2D mask array 'obmask'
%   returns array np-by-nch, where np is number of nonzeros in 'obmask',
%   nch-number of slices in 'data'
%
%  AS, 2006

nch=size(data,3);
dims=[size(data,1) size(data,2)];
if ~exist('obmask')
    obmask=ones(dims);
end
res=reshape((data(repmat(logical(obmask), [1 1 nch]))), [], nch);