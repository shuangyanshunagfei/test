
function [W,Anorms]= Angle_PCA_p3( A,dim)
%% angle PCA
%% Wang Q, Gao Q, Gao X, et al.Angle principal component analysis, IJCAI, 2017:2936-2942.
norms=sqrt(sum(A.^2,2)); %add by yi
Anorms=A./norms; %add by yi
A=Anorms'; %add by yi

 [r,c] = size(A);
 vvv = rand(r,dim);

W=orth(vvv);
numberIteration = 1;
numberIterationMax = 60;


J_value = 1;
J0 = 0;
JJ = [];
gama = 1;
Precision = 1e-8;

while  numberIteration < numberIterationMax
    AQ = sqrt(sum((W'*A).*(W'*A),1));
    BQ = sqrt(sum((A-W*W'*A).*(A-W*W'*A),1));

    id = find (BQ<gama);
    BQ(id) = BQ(id) + gama;

    Lamd = (AQ.^3)./(BQ.^3);
    H = A*(diag(1./(AQ.^(-1))+Lamd./(BQ.^(-1))))*A';
    J1 =sum(Lamd);
   
  
  [eigvector,eigvalue]=eig((H+H')/2);
eigvalue = diag(eigvalue);

[~, index] = sort((eigvalue),'descend');

W = eigvector(:,index(1:dim));
JJ = [JJ,J1];

J_value = abs(J1 - J0)/J0;
J0 = J1;

numberIteration = numberIteration +1;
end
numberIteration = numberIteration+0
plot(JJ)
save('YALEBJJPthree.mat', 'JJ')


