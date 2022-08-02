function [pd, r2, nrsd]=levenbergT2(data, te, ig, Ni, wordy, is_paral,tol, W)

% [pd, r2, nrsd]=levenbergT2(data, te, ig, Ni, wordy)
%   
%   INPUT:
%   
%   
%   OUTPUT:
%   
%   
%   
%   by Alexey Samsonov
%   2011

npts=size(data,1);
tolLmd=1e4;
if ~check_var('tol')
    tol=1e-6;
end
if ~check_var('wordy')
    wordy=0;
end

if isreal(data)
    is_cmpx=0;
else
    is_cmpx=1;
end

if ~check_var('ig')
    if ~is_cmpx
        ig=repmat([1 0], [npts 1]);
    else
        ig=repmat([1 0 1], [npts 1]); 
    end
elseif numel(ig)==1 % -- linear approximation as initial guess
    A(:,2)=te(:);
    A(:,1)=1;
    ig=repmat([1 20], [npts 1]);
    iA=pinv(A);
    tolLmd=1e1;
    
    y=abs(data);
    y(y<=0)=eps;
    tmp=(iA*log(y).').';
    ig=[exp(tmp(:,1)) -tmp(:,2)];
    
    if is_cmpx && 0
        tt=ig(:,1).*exp(1i*angle(sum(data,2)));
        ig(:,1)=real(tt);
        ig(:,3)=imag(tt);
        pd=ig(:,1)+1i*ig(:,3);
    elseif is_cmpx
        ig(:,3)=0;
    else
        pd=ig(:,1);
    end
    
    r2=ig(:,2);
    nrsd=0;
end

if ~check_var('Ni')
    Ni=10;
end

if ~check_var('W')
    W=ones(size(data));
end
    
te=te(:)';

pd=zeros(npts, 1);
r2=pd;
nrsd=pd;
nms=size(data,2);

pool=create_pool(24,1);

parfor ii=1:npts    
    ww=W(ii,:);
    if is_cmpx
        y=[real(data(ii,:)).*ww imag(data(ii,:)).*ww];
    else
        y=[real(data(ii,:)).*ww];
    end
    y(abs(y)<=0)=eps;
    
    b=ig(ii,:);
    b1=b;
    n2=norm(y(:));
    n3=n2;
    i2=1;
    lmd=1e-1;
    lmd=1;
    if is_cmpx
        J=zeros(2*nms,3);
    else
        J=zeros(nms,2);
    end
    n1=0;
    for jj=1:Ni
        E2=exp(-b(2)*te);
        if is_cmpx
            y0=[b(1)*E2.*ww b(3)*E2.*ww];
        else
            y0=[b(1)*E2.*ww];
        end
        rsd=y-y0;
        n0=norm(rsd);
        if jj>1
            if n0<n1
                lmd=0.1*lmd;
                if lmd<1e-10
                    lmd=1e-10;
                end
                b1=b;
                n2=n0;
                i2=i2+1;
            else
                lmd=10*lmd;
                if lmd>1e14
                    lmd=1e14;
                end
                b=b1;
                E2=exp(-b(2)*te);
                if is_cmpx
                    y0=[b(1)*E2.*ww b(3)*E2.*ww];
                else
                    y0=[b(1)*E2.*ww];
                end
                rsd=y-y0;
                n0=norm(rsd);
            end
        end
        % Jacobian
        J(1:nms,1)=E2';
        J(1:nms,2)=(-b(1)*te.*E2)';
        if is_cmpx
            J(1:nms,3)=0;
            J(nms+1:end,1)=0;
            J(nms+1:end,2)=(-b(3)*te.*E2)';
            J(nms+1:end,3)=E2';
        end
        %
        a=J'*diag(repmat(ww,[1 is_cmpx+1]))*J+eye(size(J,2))*lmd;
        beta=rsd*J;
        %
        dlt=pinv(a)*beta';
        b=b+dlt';
        n1=n0;
       
        if n0/n3<tol || lmd>tolLmd
        
            break;
        end
    end
    if is_cmpx
        pd(ii)=b1(1)+b1(3)*1i;
    else
        pd(ii)=b1(1);
    end
    r2(ii)=b1(2);
    nrsd(ii)=n2;
end

end
