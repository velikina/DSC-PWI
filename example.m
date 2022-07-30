%% 
clear
load example_data
[Nx,Ny,Nrcv]=size(smaps);
[Np,Nsp,Nrcv,Nte,Nfr]=size(kdata);

%% Set up time-resolved data
kxt=cell(Nfr,1);
kyt=cell(Nfr,1);
kvt=cell(Nfr,Nte);
Npass=Nfr/Nsp;
sd = 0; % noise standard deviation
% sd = 2; 
for m=1:Nte
    for kk=1:Npass
        for jj=1:Nsp
            ii=(kk-1)*Nsp+jj;
            kxt{ii}=kx(:,jj);
            kyt{ii}=ky(:,jj);           
            kvt{ii,m}=reshape(kdata(:,jj,:,m,ii),[],Nrcv);
            kvt{ii,m}=kvt{ii,m}+sd*crandn(size(kvt{ii,m}));
        end
    end
end

%%
OPT0=get_options;
OPT0.fov=1;
OPT0.Ni =[120 120];
OPT0.Nrw=5;
OPT0.tol=1e-8;

OPT1=OPT0;
OPT2=OPT0;
OPT3=OPT0;

OPT1.TD1(1)=3; % first temporal difference constraint
OPT2.TD2(1)=3; % second temporal difference constraint
OPT3.TD3(1)=3; % third temporal difference constraint
%% Image reconstruction
ims_td1_l1=zeros(Nx,Ny,Nfr,Nte);
ims_td1_l2=zeros(Nx,Ny,Nfr,Nte);
ims_td2_l1=zeros(Nx,Ny,Nfr,Nte);
ims_td2_l2=zeros(Nx,Ny,Nfr,Nte);
ims_td3_l1=zeros(Nx,Ny,Nfr,Nte);
ims_td3_l2=zeros(Nx,Ny,Nfr,Nte);
for m=1:Nte
    [ims_td1_l1(:,:,:,m),ims_td1_l2(:,:,:,m)]=irlsTD(kxt,kyt,kvt(:,m),smaps,0,OPT1);
    [ims_td2_l1(:,:,:,m),ims_td2_l2(:,:,:,m)]=irlsTD(kxt,kyt,kvt(:,m),smaps,0,OPT2);
    [ims_td3_l1(:,:,:,m),ims_td3_l2(:,:,:,m)]=irlsTD(kxt,kyt,kvt(:,m),smaps,0,OPT3);
end

%% R2* fitting
for jj=1:Nfr
    [pp,rr]=levenbergT2(selectData(sqz(abs(ims_td1_l1(:,:,jj,:))),om),TE,[],40,0,1);
    R2_td1_l1(:,:,jj)=spreadData(rr,om);
    [pp,rr]=levenbergT2(selectData(sqz(abs(ims_td1_l2(:,:,jj,:))),om),TE,[],40,0,1);
    R2_td1_l2(:,:,jj)=spreadData(rr,om);
    [pp,rr]=levenbergT2(selectData(sqz(abs(ims_td2_l1(:,:,jj,:))),om),TE,[],40,0,1);
    R2_td2_l1(:,:,jj)=spreadData(rr,om);
    [pp,rr]=levenbergT2(selectData(sqz(abs(ims_td2_l2(:,:,jj,:))),om),TE,[],40,0,1);
    R2_td2_l2(:,:,jj)=spreadData(rr,om);
    [pp,rr]=levenbergT2(selectData(sqz(abs(ims_td3_l1(:,:,jj,:))),om),TE,[],40,0,1);
    R2_td3_l1(:,:,jj)=spreadData(rr,om);
    [pp,rr]=levenbergT2(selectData(sqz(abs(ims_td3_l2(:,:,jj,:))),om),TE,[],40,0,1);
    R2_td3_l2(:,:,jj)=spreadData(rr,om);
end
%% Computing rCBV 
R2_bl=mean(R2_full(:,:,1:Nsp),3);
rCBV=sum(abs(bsxminus(R2_full(:,:,Nsp+1:end),R2_bl)),3);

R2_td1_l1_bl=mean(R2_td1_l1(:,:,1:Nsp),3);
rCBV_td1_l1=sum(abs(bsxminus(R2_td1_l1(:,:,Nsp+1:end),R2_td1_l1_bl)),3);
R2_td1_l2_bl=mean(R2_td1_l2(:,:,1:Nsp),3);
rCBV_td1_l2=sum(abs(bsxminus(R2_td1_l2(:,:,Nsp+1:end),R2_td1_l2_bl)),3);
R2_td2_l1_bl=mean(R2_td2_l1(:,:,1:Nsp),3);
rCBV_td2_l1=sum(abs(bsxminus(R2_td2_l1(:,:,Nsp+1:end),R2_td2_l1_bl)),3);
R2_td2_l2_bl=mean(R2_td2_l2(:,:,1:Nsp),3);
rCBV_td2_l2=sum(abs(bsxminus(R2_td2_l2(:,:,Nsp+1:end),R2_td2_l2_bl)),3);
R2_td3_l1_bl=mean(R2_td3_l1(:,:,1:Nsp),3);
rCBV_td3_l1=sum(abs(bsxminus(R2_td3_l1(:,:,Nsp+1:end),R2_td3_l1_bl)),3);
R2_td3_l2_bl=mean(R2_td3_l2(:,:,1:Nsp),3);
rCBV_td3_l2=sum(abs(bsxminus(R2_td3_l2(:,:,Nsp+1:end),R2_td3_l2_bl)),3);
