%Importing Data
DesignMatrix=table2array(Design);
Parcel=table2array(ROIvec);
%YM=readcsv('Yoriginal.csv');
%opts = detectImportOptions('Y_cent.csv');
%opts.DataLines = 1;
%opts.VariableNamesLine = 1;
T2 = readtable('Y_cent.csv');%,opts);
%Y_original 284 x 91146
Y_original=table2array(T2)+(table2array(ColMeans))';
%Create InverseTrans_ROI_part_Phi
Parc=unique(Parcel);
Phi_MatrixNew=zeros(size(Parc,1),size(Parcel,1));
for i=1:size(Parc,1)
  index_aux=find(Parcel==i);
  Phi_MatrixNew(i,index_aux)=1;
end

%Y_MatrixC = bsxfun(@minus, Y_original, mean(Y_original));
%[row, col] = find(isnan(table2array(T2)));
Y_MatrixC=table2array(T2);
%First and second spatial
% Clusters - Local Spatial Basis
ROI_PC=cell(size(Phi_MatrixNew,1),1);
RV_PC=cell(size(Phi_MatrixNew,1),1);
YC_Matrix_PC=cell(size(Phi_MatrixNew,1),1);

for i=1:size(Phi_MatrixNew,1)
    ind_roi=find(Phi_MatrixNew(i,:)==1);
    X=double(Y_MatrixC(:,ind_roi));
    [U,S,V] =svd(X,'econ');
    VarV=cumsum(diag(S).^2)/sum(diag(S).^2);
    idx=1:1:size(VarV,1);
    p=min(idx(VarV>=0.5));
    PCs=X*V;
    ROI_PC{i}=PCs(:,1:p);
    RV_PC{i}=V(:,1:p)';
    YC_Matrix_PC{i}=PCs(:,1:p)*V(:,1:p)';
end
X0=[];
for i=1:size(Phi_MatrixNew,1)
    X0=[X0,ROI_PC{i}];
end
% Global Spatial Basis
[U0,S0,V0] =svd(X0,'econ');
VarV0=cumsum(diag(S0).^2)/sum(diag(S0).^2);
idx0=1:1:size(VarV0,1);
p0=min(idx0(VarV0>=0.95));
PC_X0=X0*V0(:,1:p0);
%% Wavelet Transform on time domain
Y_star=PC_X0';
wavespecs0.wavelet='db1';
wavespecs0.ndim=1;
wavespecs0.compress=1;
padding=37;
Y0=[Y_star,repmat(0,size(Y_star,1),padding)]; 
%Add padding until numer of basis is even. This step makes sure that Y=XB+E is the same as
%WY=WXB+WE

[Dtrans,wavespecs]=DWT_rows(Y0,wavespecs0);
D=Dtrans';
    
%XCov = bsxfun(@minus, DesignMatrix, mean(DesignMatrix));
%XCov2 = bsxfun(@rdivide, XCov, std(DesignMatrix));
XCov2=DesignMatrix;
Xt=[XCov2;repmat(0,padding,size(XCov2,2))];

[DX,wavespecsX]=DWT_rows(Xt',wavespecs);
model.X=[ones(size(DX,2),1),DX'];
 %% Model set up
[model.n model.p]=size(model.X);
X=model.X;
mvector=[];
    for i=1:wavespecs.J
        aux0=i*ones(wavespecs.Kj(i),1);
        mvector=[mvector;aux0];
    end
W0=GetW_Michelle(model,D);    
[theta_mle0,beta_mle0,Var_beta0]=regression_mle(W0,model.n);

% Estimate scale parameters
tic
[psi_hat,alpha_hat,ConvergenceMatrix]=EstimateScale_LongMemory(D,50,model,wavespecs);
tictoc_variancecomponents=toc;

scale.newpsi=psi_hat(:,50);
scale.alpha=alpha_hat(:,50);
scale.psi=(10^10)*ones(model.p,size(beta_mle0,2));
mvector2=repmat(mvector',size(D,2),1);

% Initial values
Beta_gls=zeros(size(beta_mle0));
for k=1:size(beta_mle0,2)
scale_newpsi2=repmat(scale.newpsi(k,1),size(D,1),1); %fill up the matrix with scale values- rows are equal
scale_alpha2=repmat(scale.alpha(k,1),size(D,1),1);
SigmaInv=(scale_newpsi2.^(-1)).*(2.^(mvector.*scale_alpha2));
Beta_gls(:,k)=inv(X'*diag(SigmaInv)*X)*X'*diag(SigmaInv)*D(:,k);
end
[W_trans,Wv_trans,Vbetans_trans]=transformMichelle(D,scale,model,1,1,mvector);

% MCMC specifications and hyperparameters
MCMCspecs.nu2_psi_idx=1;
MCMCspecs.minVC=10^(-3);
max_nu2=1e-2; 
extra_scale=1e2;
B=6000;
[p,K]=size(beta_mle0);
MCMC_beta=NaN(B,p,K);
MCMC_scale_psi=NaN(B,p,K);
MCMC_scalenewpsi=NaN(B,K);
   
beta=Beta_gls;
MCMCspecs.maxO=10^20;
gamma=ones(size(beta_mle0));

%Grouping spatial basis
wavespecsPC.K=size(Y_star,1);  
V0Criteria=log2(diag(S0(1:p0,1:p0)));
cutoffsV0Crit=quantile(V0Criteria,[0.2,0.4,0.6,0.8,0.9]);
wavespecsPC.J=6;
wavespecsPC.Kj(1)=sum(V0Criteria>=cutoffsV0Crit(5));
wavespecsPC.Kj(2)=sum(V0Criteria<cutoffsV0Crit(5) & V0Criteria>=cutoffsV0Crit(4));
wavespecsPC.Kj(3)=sum(V0Criteria<cutoffsV0Crit(4) & V0Criteria>=cutoffsV0Crit(3));
wavespecsPC.Kj(4)=sum(V0Criteria<cutoffsV0Crit(3) & V0Criteria>=cutoffsV0Crit(2));
wavespecsPC.Kj(5)=sum(V0Criteria<cutoffsV0Crit(2) & V0Criteria>=cutoffsV0Crit(1));
wavespecsPC.Kj(6)=wavespecsPC.K-sum(wavespecsPC.Kj);
% Teste=[V0Criteria>=cutoffsV0Crit(5),V0Criteria<cutoffsV0Crit(5) & V0Criteria>=cutoffsV0Crit(4),V0Criteria<cutoffsV0Crit(4) & V0Criteria>=cutoffsV0Crit(3),V0Criteria<cutoffsV0Crit(3) & V0Criteria>=cutoffsV0Crit(2),V0Criteria<cutoffsV0Crit(2) & V0Criteria>=cutoffsV0Crit(1)]      
if MCMCspecs.nu2_psi_idx==1 % 1 indicating that nu2_psi depend on a and j.  
     meanop=uneqkron(wavespecsPC.Kj)*diag(1./wavespecsPC.Kj);
    nu2.psi=min(2./(Var_beta0*meanop)/extra_scale,max_nu2);    
else                     % 0 indicating that nu2_psi depend on j and k.     
    nu2.psi=min(2./mean(Var_beta0)'/extra_scale,max_nu2);
end

PiMat=ones(size(Vbetans_trans));
psiparam.a0=2;
psiparam.b0=2;

tic

for i=1:B
[W_trans,Wv_trans,Vbetans_trans]=transformMichelle(D,scale,model,1,1,mvector);
[W_trans2]=transformMichelle2(D,scale,model,1,mvector);
scale.newpsi=UpdatePsi(model,psiparam,W_trans2,K,beta);
[beta,gamma,alpha]=UpdateBetaNoOrthog_Michelle(beta,Vbetans_trans,PiMat,Wv_trans,model,wavespecsPC,MCMCspecs,scale.psi); 
MCMC_beta(i,:,:)=beta;
MCMC_scalenewpsi(i,:)=scale.newpsi;
scale.psi=10*Vbetans_trans; %scale.psi is fixed, it's the prior variance of beta;
MCMC_scale_psi(i,:,:)=scale.psi;
end
tictoc_mcmc=toc;

burnin=2000;
thin=5;
Results_MCMC=reshape(MCMC_beta(burnin:thin:end,:,:),801,p*K);
Results_MCMC_psi=reshape(MCMC_scale_psi(burnin:thin:end,:,:),801,p*K);
posterior_mean=shiftdim(mean(MCMC_beta(burnin:thin:end,:,:),1),1);
posterior_mean_scale=mean(MCMC_scalenewpsi(burnin:thin:end,:),1);

%Saving MCMC Results
%save(strcat(path2,'Results_MCMC_WMtask_CHSB'),'Results_MCMC','Results_MCMC_psi','MCMC_scalenewpsi','RV_PC','V0','p0','p','K','tictoc_mcmc','tictoc_variancecomponents')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Compute Correlations between ROIs

% Computing ROI connectivity

NRoi_new=size(Phi_MatrixNew,1);

wavelet=wavespecs.wavelet;
n=wavespecs.T;
boundary=wavespecs.boundary;
nlevels=wavespecs.nlevels;
extend=1;
[W,Kj]=Get_DWT(wavelet,wavespecs.T,boundary,1,nlevels);
Wtrans=W';

%Compute Sigma Y**
CovarianceTime=zeros(size(Y0,2),size(Y0,2),size(Y_star,1));
newpsi_posteriormean=mean(MCMC_scalenewpsi);
for k=1:size(Y_star,1)
scale_newpsi2=repmat(newpsi_posteriormean(1,k),size(D,1),1); %fill up the matrix with scale values- rows are equal
scale_alpha2=repmat(alpha_hat(k,50)',size(D,1),1);
Sigma=round(scale_newpsi2.*(2.^(mvector.*(-scale_alpha2))),7);
CovarianceTime(:,:,k)=Wtrans*diag(Sigma)*Wtrans';
end

ThetaTime=zeros(size(Y_star,1),1);
for k=1:size(Y_star,1)
   %if range(diag(CovarianceTime(:,:,k)))<=10^(10)
       ThetaTime(k,1)=mode(diag(CovarianceTime(:,:,k)));
   %else
   %    ThetaTime(k,1)=NaN;
   %end
end
 
indexpc=cell(size(Y_star,1),1);
npc_cumu=0;
for i=1:NRoi_new
    sa=npc_cumu;
    npc_cumu=npc_cumu+size(RV_PC{i},1);
    indexpc{i}=sa+1:npc_cumu;
end
RV=zeros(NRoi_new,NRoi_new); %Matrix with correlation between regions
for i=2:NRoi_new
    for j=1:i
    index_combined=[indexpc{i},indexpc{j}];
    ni=size(indexpc{i},2);
    nj=size(indexpc{j},2);
    CovThetaPC_comb=V0(index_combined',1:p0)*(diag(ThetaTime)*eye(p0))*(V0(index_combined',1:p0))';
    CovThetaPC_ij=CovThetaPC_comb(ni+1:end,1:ni); 
    CovThetaPC_i=CovThetaPC_comb(1:ni,1:ni);
    CovThetaPC_j=CovThetaPC_comb(1+ni:end,1+ni:end);
    CorThetaPC_comb=inv(sqrtm(CovThetaPC_i))*CovThetaPC_ij'*inv(sqrtm(CovThetaPC_j));
    %CorThetaPC_i=inv(sqrtm(CovThetaPC_i))*CovThetaPC_i'*inv(sqrtm(CovThetaPC_i));
    %CorThetaPC_j=inv(sqrtm(CovThetaPC_j))*CovThetaPC_j'*inv(sqrtm(CovThetaPC_j));
    RV(i,j)=trace(CorThetaPC_comb*CorThetaPC_comb')/sqrt(ni*nj);
    %RV(i,j)=trace(CorThetaPC_comb*CorThetaPC_comb')/sqrt(trace(CorThetaPC_i*CorThetaPC_i')*trace(CorThetaPC_j*CorThetaPC_j'));
    RV(j,i)=RV(i,j);
    end
end

Connectivity=RV;
path3='/Users/michellemiranda/Documents/MATLAB/Mikayla/';
writematrix(Connectivity, strcat(path3,"RV_Coef.txt"));
figure
imagesc(Connectivity)
title('Estimated Background Connectivity')
colormap('jet')
colorbar('FontSize',12)
%saveas(gcf,strcat(path2,'ConnectivityNotOrdered.png'))
