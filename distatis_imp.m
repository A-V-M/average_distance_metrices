function [compromise compromiseRDM tau_RV eigenVectors] = distatis_imp(D)
%  D{1} = D_csl_vSpace_48;
%  D{2} = Dcategories_RDM_48;
%  

%D = voxRDM;

nmodels = length(D);
if nmodels > 1; n = length(D{1});else n = length(D);end
 %step 1: cross-product matrix
 for nModel=1:nmodels
 m=ones(n,1)./n;
 cent=eye(n,n)-ones(n,1)*m';
 Dsum=zeros(n,n);
 
 cp_mat = -(1/2)*cent*D{nModel}*cent;
 
 [PS,lS]=eigen(cp_mat);
 S{nModel}=cp_mat.*(lS(1).^(-1)); %normalise by first eigen
 end
 
 %step 2: RV coefficient (cosine)
 
 nn = n*n;
 
 cp_mat_all = zeros(nn,2);
 
 for nModel=1:nmodels
     cp_mat_all(:,nModel) = S{nModel}(:);
 end
 
 num_RV = cp_mat_all'*cp_mat_all;
 diag_val = diag(num_RV);
 nd = length(diag_val);
 den_RV=(repmat(diag_val',nd,1).*repmat(diag_val,1,nd)).^(1/2);
 mat_RV=num_RV./den_RV;
 eig_values = eig(mat_RV);
 
 [P,phi]=eigen(mat_RV);
 G=P*diag(phi.^(1/2));
 
 eigenVectors = P;
 tau_RV = phi;
 tau_RV_perc=round(100 *(phi./sum(phi)));
 
% step 3: Compute the alpha weights
weights=P(:,1)/sum(abs(P(:,1)));


% step 4: Compute the compromise
compromise=zeros(n,n);
for nModel=1:nmodels;
  compromise=compromise+weights(nModel)*S{nModel};
end
compromiseRDM = zeros(n,n);
for x=1:n;for y=1:n; compromiseRDM(x,y) = (compromise(x,:)*compromise(x,:)' + compromise(y,:)*compromise(y,:)') - 2*dot(compromise(x,:),compromise(y,:)); end; end
