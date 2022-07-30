function res=fwd_mlt(f, smaps,dirs)

% Implementation of Forward multiplication Ef for pMRI reconstruction
%   
%
%   Alexey Samsonov
%   U of W, April 2005
Nc=size(smaps,3);
Ns=size(f,3);

if ~check_var('dirs')
    dirs=[1 2];
end

if numel(smaps)*Ns*8>200e8
    res=zeros([size(smaps) Ns]);
    for ii=1:Ns
        tmp=smaps.*repmat(f(:,:,ii),[1 1 Nc]);
        res(:,:,:,ii)=ffts(tmp, dirs, 1);
    end
    disp('Long');
else
    tmp=repmat(smaps,[1 1 1 Ns]).*repmat(reshape(f,[size(f,1) size(f,2) 1 Ns]),[1 1 Nc 1]);
    res=ffts(tmp, dirs, 1);
end
return;