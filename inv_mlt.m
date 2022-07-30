function res=inv_mlt(f, smaps,dirs)

% Implementation of Inverse multiplication E'f for pMRI reconstruction
%   res=inv_mlt(f, smaps,dirs)
%
%   Alexey Samsonov
%   U of W, April 2005

if ~check_var('dirs')
    tmp=ffts(f,[1 2],-1);
else
    tmp=ffts(f,dirs,-1);
end

% dims=size(f);
Ns=size(f,4);

% tmp0=zeros([dims]);

tmp0=repmat(conj(smaps),[1 1 1 Ns]).*tmp;
% for ii=1:Ns
%     tmp0(:,:,:,ii)=conj(smaps).*tmp(:,:,:,ii);
% end
res=squeeze(sum(tmp0,3));

return;