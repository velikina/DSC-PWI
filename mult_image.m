function res=mult_image(data, cc)

%   res=mult_image(data, cc);
%       Multiplication of stack of images with one image
%           data  - stack of images 
%           cc      - image to multiply
%
%       by Alexey Samsonov

dims=size(data);
if numel(dims)<3
    N1=1;
else
    N1=dims(3:end);
end

dims1=size(cc);
if numel(dims1)<3
    N2=1;
else
    N2=dims1(3:end);
end

res=zeros([dims(1:2) N1 N2]);

if numel(N1)==1
    for ii=1:N2
        res(:,:,:,ii) = data.*repmat(cc(:,:,ii),[1 1 N1]);
    end
else
    for ii=1:N2
        res(:,:,:,:,ii) = data.*repmat(cc(:,:,ii),[1 1 N1]);
    end
end
    
end


