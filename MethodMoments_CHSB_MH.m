psi_hat=zeros(size(D,2),Nite);
alpha_hat=zeros(size(D,2),Nite);
mvector=[];
    for i=1:wavespecs.J
        %aux0=(wavespecs.J-i+1)*ones(wavespecs.Kj(i),1);
        aux0=i*ones(wavespecs.Kj(i),1);
        mvector=[mvector;aux0];
    end
    cumsumwav=cumsum(wavespecs.Kj);
for k=1:size(D,2)
   % E=(eye(size(D,1))-X(:,1)*inv(X(:,1)'*X(:,1))*X(:,1)')*D(:,k);
     E=(eye(size(D,1))-X*inv(X'*X)*X')*D(:,k);
    
    Ej=cell(wavespecs.J,1);
    thetaj=zeros(wavespecs.J,1);
    Ej{1}=E(1:cumsumwav(1),1);
    for j=2:wavespecs.J
    Ej{j}=E(cumsumwav(j-1)+1:cumsumwav(j),1);
    end
    for j=1:wavespecs.J
        thetaj(j,1)=var(Ej{j});
    end
    lb = [-10,0];
    ub = [40,1];
    A = ones(1,2);
    b = 10^2;
    Options=optimset('Display','none');
    coef = lsqlin([ones(wavespecs.J,1) -log(2)*(1:1:wavespecs.J)'],log(thetaj),A,b,[],[],lb,ub,[0;0],Options);
    %coef=regress(log(thetaj),[ones(wavespecs.J,1) -log(2)*(1:1:9)']);
    psi_hat(k,1)=exp(coef(1));
    alpha_hat(k,1)=coef(2);
 end
 


delta=0.01;  
ConvergenceMatrix=zeros(size(D,2),Nite-1);
for ite=2:Nite
for k=1:size(D,2)
mvector2=repmat(mvector',size(D,2),1);
scale_newpsi2=repmat(psi_hat(k,ite-1),size(D,1),1); %fill up the matrix with scale values- rows are equal
scale_alpha2=repmat(alpha_hat(k,ite-1)',size(D,1),1);
SigmaInv=(scale_newpsi2.^(-1)).*(2.^(mvector2(k,:)'.*scale_alpha2));
%Sigma=(scale_newpsi2).*(2.^(-mvector2(k,:)'.*scale_alpha2));
%E_gls=(eye(size(D,1))-X(:,1)*inv(X(:,1)'*diag(SigmaInv)*X(:,1))*X(:,1)'*diag(SigmaInv))*D(:,k);
E_gls=(eye(size(D,1))-X*inv(X'*diag(SigmaInv)*X)*X'*diag(SigmaInv))*D(:,k);

Ej=cell(wavespecs.J,1);
    thetaj=zeros(wavespecs.J,1);
    Ej{1}=E_gls(1:cumsumwav(1),1);
    for j=2:wavespecs.J
    Ej{j}=E_gls(cumsumwav(j-1)+1:cumsumwav(j),1);
    end
    for j=1:wavespecs.J
        thetaj(j,1)=var(Ej{j});
    end
    lb = [-10,0];
    ub = [40,1];
    A = ones(1,2);
    b = 10^2;
    coef = lsqlin([ones(wavespecs.J,1) -log(2)*(1:1:wavespecs.J)'],log(thetaj),A,b,[],[],lb,ub);
    %coef=regress(log(thetaj),[ones(wavespecs.J,1) -log(2)*(1:1:9)']);
    psi_hat(k,ite)=exp(coef(1));
    alpha_hat(k,ite)=coef(2);
    if (psi_hat(k,ite)-psi_hat(k,ite-1))<=100 && (alpha_hat(k,ite)-alpha_hat(k,ite-1))<=delta
    ConvergenceMatrix(k,ite-1)=1;
    end
    %end
end
end


writematrix(psi_hat,"C:\Users\Owner\Documents\Graduate Documents\Research\HCP\Data\psi_hat50.csv")
writematrix(alpha_hat,"C:\Users\Owner\Documents\Graduate Documents\Research\HCP\Data\alpha_hat50.csv")