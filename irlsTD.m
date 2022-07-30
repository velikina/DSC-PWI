function [res,im1]=irlsTD(kx0t,ky0t,v0t,smaps,im0t,OPT)

% IMPLEMENTATION OF ITERATIVELE REWEIGHTED MINIMIZATION
%
%  [res, im1]=irlsTD(kx0t,ky0t,v0t,smaps,im0t,OPT)
%   INPUT:
%           kx0t -- cell with kx coordinates for each time frame
%           ky0t -- cell with ky coordinates for each time frame
%           vot  -- cell with k-space values (Npoints x Nreceivers) for
%           each time frame
%           smaps -- sensitivity maps
%           im0  -- initial guess (optional)
%           OPT  -- reconstruction options (defined in get_options)
%   OUTPUT:
%           res -- reconstruction resut
%           im1   - image after first reweigted iteration
%
%   Alexey Samsonov, Julia Velikina 2004-2021
%

dims= [size(smaps,1) size(smaps,2)];
Nc  = size(smaps,3);
Nfr = numel(v0t);

if ~check_var('im0t')
    im0t = zeros([dims(1:2) Nfr]);
else
    if isscalar(im0t)
        im0t=im0t*ones([dims(1:2) Nfr]);
    end
end

wordy = OPT.wordy;
L     = OPT.L;
Ni    = OPT.Ni;
Nrw   = OPT.Nrw;
tol   = OPT.tol;
tolRW = OPT.tolRW;
FOV   = OPT.fov;
of    = OPT.of;

% Temporal filter (First derivative)
isTD1_1=OPT.TD1(1);
isTD1_2=OPT.TD1(2);

% Temporal filter (Second derivative)
isTD2_1=OPT.TD2(1);
isTD2_2=OPT.TD2(2);

% Temporal filter (Third derivative)
isTD3_1=OPT.TD3(1);
isTD3_2=OPT.TD3(2);

% Temporal filter (Fourth derivative)
isTD4_1=OPT.TD4(1);
isTD4_2=OPT.TD4(2);

kxt=cell(Nfr+6,1);
kyt=cell(Nfr+6,1);
kvt=cell(Nfr+6,1);
kxt(4:end-3)=kx0t;
kyt(4:end-3)=ky0t;
kvt(4:end-3)=v0t;
kxt(1:3)=kx0t(4:-1:2);
kyt(1:3)=ky0t(4:-1:2);
kvt(1:3)=v0t(4:-1:2);
kxt(end-2:end)=kx0t(end-1:-1:end-3);
kyt(end-2:end)=ky0t(end-1:-1:end-3);
kvt(end-2:end)=v0t(end-1:-1:end-3);
kx0t=kxt;
ky0t=kyt;
v0t=kvt;
Nfr=numel(v0t);



imt=zeros([dims(1:2) Nfr]);
imt(:,:,4:end-3)=im0t;
imt(:,:,1:3)=im0t(:,:,4:-1:2);
imt(:,:,end-2:end)=im0t(:,:,end-1:-1:end-3);
im0t=imt;



traj   =OPT.traj;
switch lower(traj)
    case 'cart'
        isCart=1;
        cc1=ones(dims);
        of=1;
        dimsOG=ceil(dims*of);
        indtf=[];
        for ii=1:Nfr
            ind = sub2ind(dims,kx0t{ii},ky0t{ii});
            indtf = [indtf; (ii-1)*prod(dims)*Nc+vect(repmat(ind(:),[1 Nc])+repmat(prod(dims)*([1:Nc]-1), [numel(ind) 1]))];
        end
        v0=[];
        for ii=1:Nfr
            tmp=vect(v0t{ii});
            rr=[1:numel(tmp)]+numel(v0);
            v0(rr)=tmp;
        end
        dks=1;
        cnt=ceil(dimsOG/2);
        rt=[-cnt(1)+even1(dimsOG(1)):cnt(1)]-1;
        [kx, ky] = meshgrid(rt,rt);
        rk=dks*sqrt(kx.^2+ky.^2);
        kx  = dks*kx(:);
        ky  = dks*ky(:);
        fh0=fftshifts(0.5+atan(200*(1-rk/max(kx(:))))/pi);
    otherwise
        isCart=0;
        dk=1/FOV;
        dks=dk/of;
        dimsOG=ceil(dims*of);
        cnt=ceil(dimsOG/2);
        rt=[-cnt(1)+even1(dimsOG(1)):cnt(1)]-1;
        [kx, ky] = meshgrid(rt,rt);
        rk=dks*sqrt(kx.^2+ky.^2);
        kx=dks*kx(:);
        ky=dks*ky(:);
        fh0=(0.5+atan(200*(1-rk/max(kx(:))))/pi);
        
        % Setting up convolution structures
        for ii=1:Nfr
            [tmpAC, wlut] = SetupMRIConvData(kx0t{ii},ky0t{ii},kx,ky, dk,'kaiser', L,[]);
            [tmpCA, wlut] = SetupMRIConvData(kx,ky,kx0t{ii},ky0t{ii}, dk,'kaiser', L,[]);
            cnvstrAC{ii}=tmpAC;
            cnvstrCA{ii}=tmpCA;
            dns{ii}=ones(size(kx0t{ii}, 1)*size(kx0t{ii},2), 1);
        end
        
        B=pi*sqrt(4*L*L/(4)*(2-0.5)^2-0.8);
        cc1 = get_KB_imspace1([fliplr(wlut(2:end)) wlut], dims, 2*L);
        cc1=cc1/max(cc1(:));
end

df=ceil((dimsOG-dims)/2);
r1=[1:dims(1)]+df(1);
r2=[1:dims(2)]+df(2);

smapsCC0=zeros([dimsOG Nc]);
smapsCC0(r1,r2,:)=mult_image(smaps,imdiv(cc1));
smapsCC=fftshifts(smapsCC0);
smapsOG=zeros([dimsOG Nc]);
smapsOG(r1,r2,:)=smaps;
smapsOG=fftshifts(smapsOG);

omOG=zeros(dimsOG);
omOG(r1,r2)=ones(dims(1:2));
omOG=fftshifts(omOG);

%  --  scale for regularization
scale = double(sqrt(sum(abs(smaps(:)).^2/numel(smaps))));
mm=false([dimsOG Nfr]);
mm1=mm;
mm(r1,r2,:)=1;
mm=fftshifts(mm);
mm1(r1(1)+2:r1(end)-2,r2(1)+2:r2(end)-2,:)=1;

d=zeros([dimsOG Nfr]);
cb=chessboard(dimsOG);

res=im0t;

nTD1_0=0;
nTD2_0=0;
nTD3_0=0;
nTD4_0=0;

for iRW=1:Nrw
    d(:)=0;
    nr0=norm(res(:));
    
    d(r1,r2,:) = res;
    d=fftshifts(d);
    if ~strcmp(traj,'cart')
        d=mult_image(d,cb);
    end
    
    if iRW<3
        sgmTD1=[];
        sgmTD2=[];
        sgmTD3=[];
        sgmTD4=[];
    end
    
    %******************************************************
    % -- TD reweighting (First derivative)
    WTD1=isTD1_1^2*ones(size(d)+[0 0 2]);
    d2=d;
    if isTD1_1 && nr0~=0
        d3=zeros(size(d2)+[0 0 2]);
        d3(:,:,1)=1*(d2(:,:,1)*1+0*d2(:,:,3));
        d3(:,:,end)=1*(d2(:,:,end)*1+0*d2(:,:,end-2));
        d3(:,:,2:end-1)=d2;
        cr=zeros(size(d3));
        cr(:,:,1:end-1)=diff(d3,1,3);
        WR1=abs(cr).^2;
        sgmTD1=std(sqrt(WR1(:)))*0.6;
        if sgmTD1==0
            W=1;
        else
            W=1./sqrt(WR1/sgmTD1.^2+1);
        end
        WTD1=WTD1.*W;
    end
    WTD1=WTD1+isTD1_2^2;
    nTD1=norm(WTD1(:));
    
    % - TD reweighting (second derivative)
    WTD2=isTD2_1^2*ones(size(d)+[0 0 4]);
    if isTD2_1 && nr0~=0
        d3=zeros(size(d)+[0 0 4]); d3(:,:,1)=d(:,:,1)*3-d(:,:,2)*2; d3(:,:,2)=d(:,:,1)*2-d(:,:,2); d3(:,:,end-1)=d(:,:,end)*2-d(:,:,end-1); d3(:,:,end)=d(:,:,end)*3-d(:,:,end-1)*2;d3(:,:,3:end-2)=d;
        cr=zeros(size(d3));
        cr(:,:,2:end-1)=diff(d3,2,3);
        WR1=abs(cr).^2;
        sgmTD2=std(sqrt(WR1(:)))*0.6;
        if sgmTD2==0
            W=1;
        else
            W=1./sqrt(WR1/sgmTD2.^2+1);
        end
        WTD2=WTD2.*W;
    end
    WTD2=WTD2+isTD2_2^2;
    nTD2=norm(WTD2(:));
    
    % - TD reweighting (third derivative)
    WTD3=isTD3_1^2*ones(size(d)+[0 0 6]);
    if isTD3_1 && nr0~=0
        d3=zeros(size(d)+[0 0 6]);
        d3(:,:,1)=10*d(:,:,1)-15*d(:,:,2)+6*d(:,:,3);
        d3(:,:,2)=6*d(:,:,1)-8*d(:,:,2)+3*d(:,:,3);
        d3(:,:,3)=3*d(:,:,1)-3*d(:,:,2)+d(:,:,3);
        d3(:,:,end-2)=d(:,:,end-2)-3*d(:,:,end-1)+3*d(:,:,end);
        d3(:,:,end-1)=3*d(:,:,end-2)-8*d(:,:,end-1)+6*d(:,:,end);
        d3(:,:,end)=6*d(:,:,end-2)-15*d(:,:,end-1)+10*d(:,:,end);
        d3(:,:,4:end-3)=d;
        cr=zeros(size(d3));
        cr(:,:,3:end-1)=diff(d3,3,3);
        WR1=abs(cr).^2;
        sgmTD3=std(sqrt(WR1(:)))*0.6
        if sgmTD3==0
            W=1;
        else
            W=1./sqrt(WR1/sgmTD3.^2+1);
        end
        WTD3=WTD3.*W;
    end
    WTD3=WTD3+isTD3_2^2;
    nTD3=norm(WTD3(:));
    
    % - TD reweighting (forth derivative)
    WTD4=isTD4_1^2*ones(size(d)+[0 0 8]);
    if isTD4_1 && nr0~=0
        d3=zeros(size(d)+[0 0 8]); d3(:,:,5:end-4)=d;
        d3(:,:,4)=4*d(:,:,1)-6*d(:,:,2)+4*d(:,:,3)-1*d(:,:,4);
        d3(:,:,3)=10*d(:,:,1)-20*d(:,:,2)+15*d(:,:,3)-4*d(:,:,4);
        d3(:,:,2)=20*d(:,:,1)-45*d(:,:,2)+36*d(:,:,3)-10*d(:,:,4);
        d3(:,:,1)=35*d(:,:,1)-84*d(:,:,2)+70*d(:,:,3)-20*d(:,:,4);
        d3(:,:,end)=35*d(:,:,end)-84*d(:,:,end-1)+70*d(:,:,end-2)-20*d(:,:,end-3);
        d3(:,:,end-1)=20*d(:,:,end)-45*d(:,:,end-1)+36*d(:,:,end-2)-10*d(:,:,end-3);
        d3(:,:,end-2)=10*d(:,:,end)-20*d(:,:,end-1)+15*d(:,:,end-2)-4*d(:,:,end-3);
        d3(:,:,end-3)=4*d(:,:,end)-6*d(:,:,end-1)+4*d(:,:,end-2)-1*d(:,:,end-3);
        cr=zeros(size(d3));
        cr(:,:,3:end-2)=diff(d3,4,3);
        WR1=abs(cr).^2;
        sgmTD4=std(sqrt(WR1(:)))*0.6;
        if sgmTD4==0
            W=1;
        else
            W=1./sqrt(WR1/sgmTD4.^2+1);
        end
        WTD4=WTD4.*W;
    end
    WTD4=WTD4+isTD4_2^2;
    nTD4=norm(WTD4(:));
    
    x=d;
    if isCart
        s=zeros([dims Nc Nfr]);
        s(indtf)=v0;
        s=mult_image(s,cb);
        u1a=inv_mlt(s, smapsCC).*mm;
    else
        tmpU1=zeros([dimsOG Nc Nfr]);
        for jj=1:Nfr
            tmpU1(:,:,:,jj)=convmridata(v0t{jj},cnvstrAC{jj}, dns{jj},[dimsOG Nc]);
        end
        u1a=inv_mlt(tmpU1, smapsCC).*mm;
    end
    
    u1=u1a;
    
    % -- Begin CG iterations
    if iRW>numel(Ni)
        Ni0=Ni(end);
    else
        Ni0=Ni(iRW);
    end
    
    if iRW>numel(tol)
        tol0=tol(end);
    else
        tol0=tol(iRW);
    end
    
    for ii = 1:Ni0
        % !
        d00=d;
        tmp=fwd_mlt(d00,smapsCC);
        
        if isCart
            s=tmp;
            s(:)=0;
            s(indtf)=tmp(indtf);
            q0=inv_mlt(s , smapsCC).*mm0;
        else
            for jj=1:Nfr
                s=convmridata(selectData(tmp(:,:,:,jj)), cnvstrCA{jj});
                tmp1=convmridata(s,cnvstrAC{jj}, dns{jj},[dimsOG Nc]);
                q0(:,:,jj)=inv_mlt(tmp1, smapsCC).*mm(:,:,jj);
            end
        end
        
        q=q0;
        
        % ***********************
        % -- applying temporal filters
        if nTD1~=0
            d2=d;
            d3=zeros(size(d2)+[0 0 2]);
            d3(:,:,1)=1*(d2(:,:,1)*1+0*d2(:,:,3));
            d3(:,:,2:end-1)=d2;
            d3(:,:,end)=1*(d2(:,:,end)*1+0*d2(:,:,end-2));
            cr=zeros(size(d3));
            cr(:,:,1:end-1)=diff(d3,1,3);
            cr=flipdim(cr.*WTD1,3);
            cr(:,:,1:end-1)=diff(cr,1,3);
            crt=flipdim(cr,3);
            q=q+crt(:,:,2:end-1);
        end
        
        if nTD2~=0
            d3=zeros(size(d)+[0 0 4]);
            d3(:,:,1)=d(:,:,1)*3-d(:,:,2)*2;
            d3(:,:,2)=d(:,:,1)*2-d(:,:,2);
            d3(:,:,3:end-2)=d;
            d3(:,:,end-1)=d(:,:,end)*2-d(:,:,end-1);
            d3(:,:,end)=d(:,:,end)*3-d(:,:,end-1)*2;
            cr=zeros(size(d3));
            cr(:,:,2:end-1)=diff(d3,2,3);
            cr=flipdim(cr.*WTD2,3);
            cr(:,:,2:end-1)=diff(cr,2,3);
            crt=flipdim(cr,3);
            q=q+crt(:,:,3:end-2);
        end
        
        if nTD3~=0
            d3=zeros(size(d)+[0 0 6]);
            d3(:,:,1)=10*d(:,:,1)-15*d(:,:,2)+6*d(:,:,3);
            d3(:,:,2)=6*d(:,:,1)-8*d(:,:,2)+3*d(:,:,3);
            d3(:,:,3)=3*d(:,:,1)-3*d(:,:,2)+d(:,:,3);
            d3(:,:,4:end-3)=d;
            d3(:,:,end-2)=d(:,:,end-2)-3*d(:,:,end-1)+3*d(:,:,end);
            d3(:,:,end-1)=3*d(:,:,end-2)-8*d(:,:,end-1)+6*d(:,:,end);
            d3(:,:,end)=6*d(:,:,end-2)-15*d(:,:,end-1)+10*d(:,:,end);
            cr=zeros(size(d3));
            cr(:,:,3:end-1)=diff(d3,3,3);
            cr=flipdim(cr.*WTD3,3);
            cr(:,:,3:end-1)=diff(cr,3,3);
            crt=flipdim(cr,3);
            q=q+crt(:,:,4:end-3);
        end
        
        if nTD4~=0
            d3=zeros(size(d)+[0 0 8]);
            d3(:,:,1)=35*d(:,:,1)-84*d(:,:,2)+70*d(:,:,3)-20*d(:,:,4);
            d3(:,:,2)=20*d(:,:,1)-45*d(:,:,2)+36*d(:,:,3)-10*d(:,:,4);
            d3(:,:,3)=10*d(:,:,1)-20*d(:,:,2)+15*d(:,:,3)-4*d(:,:,4);
            d3(:,:,4)=4*d(:,:,1)-6*d(:,:,2)+4*d(:,:,3)-1*d(:,:,4);
            d3(:,:,5:end-4)=d;
            d3(:,:,end-3)=4*d(:,:,end)-6*d(:,:,end-1)+4*d(:,:,end-2)-1*d(:,:,end-3);
            d3(:,:,end-2)=10*d(:,:,end)-20*d(:,:,end-1)+15*d(:,:,end-2)-4*d(:,:,end-3);
            d3(:,:,end-1)=20*d(:,:,end)-45*d(:,:,end-1)+36*d(:,:,end-2)-10*d(:,:,end-3);
            d3(:,:,end)=35*d(:,:,end)-84*d(:,:,end-1)+70*d(:,:,end-2)-20*d(:,:,end-3);
            cr=zeros(size(d3));
            cr(:,:,3:end-2)=diff(d3,4,3);
            cr=flipdim(cr.*WTD4,3);
            cr(:,:,3:end-2)=diff(cr,4,3);
            crt=flipdim(cr,3);
            q=q+crt(:,:,5:end-4);
        end
        
        % -- CG calculations
        if ii==1
            u = u1 - q;
            d = u;
            d0 = u(:)'*u(:);
            if iRW==1
                u0=d0;
            end
        else
            ndq=d(:)'*q(:);
            alpha = d0./ndq;
            u = u - q.*alpha;
            d1=u(:)'*u(:);
            x = x + d.*alpha;
            d = u + d.*(d1./d0);
            d0 = d1;
        end
        
        err=abs(d0./u0);
        
        x1=fftshifts(x);
              
        rnrm=norm(u(:))/sqrt(numel(u));
        fprintf(1,'Iteration #%d, Error %g, iRW=%d, rnrm=%g\n',ii,err,iRW,rnrm);
        
        if err<tol0
            break;
        end
    end
    
    if iRW==Nrw
        x1=ffts(ffts(x, [], 1).*repmat(fh0, [1 1 Nfr]), [], -1);
        x1=mult_image(x1,srsos(smapsOG));
    else
        x1=x;
    end
    
    res=fftshifts(x1);
    res=res(r1,r2,:);
    
    
    
    if ~isCart
        res=mult_image(res,cb(r1,r2));
    end
    
    if iRW==1 && Nrw>1
        clear im1
        im1=ffts(ffts(x, [], 1).*repmat(fh0, [1 1 Nfr]), [], -1);
        im1=mult_image(im1,srsos(smapsOG));
        im1=fftshifts(im1);
        im1=im1(r1,r2,:);
        if ~isCart
            im1=mult_image(im1,cb(r1,r2));
        end
    end
    
    
    %  --  Convergence
    if iRW>1
        if nTD1
            if nTD1_0
                errTD1=norm(vect(WTD1-WTD1_0))/nTD1_0;
            else
                nTD1_0=nTD1;
                errTD1=1;
            end
        else
            errTD1=0;
        end
        
        if nTD1
            fprintf(1,'errTD1=%f  \n',errTD1);
            if errTD1<tolRW
                break;
            end
        end
        WTD1_0=WTD1;
        
        if nTD2
            if nTD2_0
                errTD2=norm(vect(WTD2-WTD2_0))/nTD2_0;
            else
                nTD2_0=nTD2;
                errTD2=1;
            end
        else
            errTD2=0;
        end
        if nTD2
            fprintf(1,'errTD2=%f  \n',errTD2);
            if errTD2<tolRW
                break;
            end
        end
        WTD2_0=WTD2;
        
        if nTD3
            if nTD3_0
                errTD3=norm(vect(WTD3-WTD3_0))/nTD3_0;
            else
                nTD3_0=nTD3;
                errTD3=1;
            end
        else
            errTD3=0;
        end
        if nTD3
            fprintf(1,'errTD3=%f  \n',errTD3);
            if errTD3<tolRW
                break;
            end
        end
        WTD3_0=WTD3;
        
        if nTD4
            if nTD4_0
                errTD4=norm(vect(WTD4-WTD4_0))/nTD4_0;
            else
                nTD4_0=nTD4;
                errTD4=1;
            end
        else
            errTD4=0;
        end
        if nTD4
            fprintf(1,'errTD4=%f  \n',errTD4);
            if errTD4<tolRW
                break;
            end
        end
        WTD4_0=WTD4;
        
        
    end
    
end

if iRW<Nrw || ((Nrw==1) && (ii<Ni(1)))
    x1=ffts(ffts(x1, [], 1).*repmat(fh0, [1 1 Nfr]), [], -1);
    x1=mult_image(x1,srsos(smapsOG));
    res=fftshifts(x1);
    if ~isCart
        res = mult_image(res, cb);
    end
    res=res(r1,r2,:);
end
res=res(:,:,4:end-3);

if ~exist('im1','var')
    im1=[];
else
    im1=im1(:,:,4:end-3);
end

end
