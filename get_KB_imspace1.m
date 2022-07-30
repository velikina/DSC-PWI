function CC=get_KB_imspace1(krn, dims, W);

%   CC=get_KB_imspace1(krn, dims, W);
%       
%       krn    - presampled kernel vector (N-by-1)
%       dims  - final image dimensions
%       W - width of interpolation kernel
%       of - overgridding factor
%
%   Alexey Samsonov
%   Univesity of Wisconsin,Madison
%   March 2005

tt=zeros(floor(numel(krn)/(W)*dims(1)),1);
tt(1:numel(krn))=krn;
tt=circshift(tt,[ceil(numel(tt)/2)-ceil(numel(krn)/2) 0]);
cc=ifft(tt);
dap=sinc(linspace(-1,1,dims(1)+1))';
dap=1;
% cc=fftshift(abs(cc([1:dims(1)/2+1  end-dims(1)/2+1:end]))).*dap;
sz=ceil(dims(1)/2*sqrt(2));
cc=fftshift(abs(cc([1:sz+1  end-sz+1:end]))).*dap;
res=zeros(dims);
rr1=-dims(1)/2:dims(1)/2;
rr2=-dims(2)/2:dims(2)/2;
[xx yy]=meshgrid(rr2,rr1);
npr=floor(dims(1)*pi/2);
[X Y]=getRadTr(floor(dims(1)*pi/2),2*sz+1,-sz, sz,0,pi);
Z=repmat(cc(:),[1 npr]);
warning off
CC=griddata(X(:),Y(:),Z(:),xx(:),yy(:));
warning on
CC=reshape(CC,dims+1);
CC=imresize(CC,dims,'bilinear');
CC(find(isnan(CC)))=min(CC(find(~isnan(CC))));
