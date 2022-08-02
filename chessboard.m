function res=chessboard(dims,dr);

%  res=chessboard(dims);
%  dims contains 1 or 2 elements to define image size
% 
%  AS, UofWM,2006

if ~exist('dr','var') | isempty(dr)
    dr=[1 1];
elseif numel(dr)==1
    dr=[~even1(dr) even1(dr)];
end

if numel(dims)==1
    nx=dims;
    ny=dims;
elseif numel(dims)==2
    nx=dims(1);
    ny=dims(2);
end

res=[-1].^(dr(1)*[1:nx]')*[-1].^(dr(2)*[1:ny]);
end
