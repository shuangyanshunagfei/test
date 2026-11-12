
function [W,Anorms]= Angle_PCA_p6( A,dim)
%% square l3 angle PCA
%% Shuangyan Yi, Sisi Wang, Feiping Nie, et al. 2026.
norms=sqrt(sum(A.^2,2));
Anorms=A./norms;%A是归一化后的数据
A=Anorms';

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
    BQ=sqrt(sum((A-W*W'*A).*(A-W*W'*A),1));
    AQ=sqrt(sum((W'*A).*(W'*A),1));
    id = find (BQ<gama);
    BQ(id)=BQ(id)+gama;

    Lamd = (AQ.^(6))./(BQ.^(6));
    H = A*(diag(1./(AQ.^(-4))+Lamd./(BQ.^(-4))))*A';
    J1 =sum(Lamd); %目标函数的值
  
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
JJJ=JJ;
save('JJJ.mat', 'JJJ');
 % figure;
 %    plot(JJ, 'b-', 'LineWidth', 1.5);
 %    xlabel('迭代次数');
 %    ylabel('目标函数值');
 %    title(['收敛曲线 (dim = ', num2str(dim), ')']);
 %    grid on;
 % 
 %    saveas(gcf, sprintf('JJ_dim_%d.png', dim));
 %    close(gcf);



% figure(2)
% hold on
% plot(JJ2)

