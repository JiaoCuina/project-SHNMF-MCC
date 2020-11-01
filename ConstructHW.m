function [Lap]=ConstructHW(X,n,q)
%% 1
% N = size(X,1);
% squared_X = sum(X.*X,2);
% transi_X = X*X';    %x_j * x_i^T;
% X_i = repmat(squared_X,1,N); % same value in row
% X_j = repmat(squared_X',N,1); % same value in col
% E1 = -(X_i + X_j - 2*transi_X );
% Dis = exp(E1/(2*0.5));
% Dis = Dis - diag(diag(Dis));
% 
 [mm,nn]=size(X);
% 
% for an=1:1:mm
%     for am=1:1:mm
%         Dis(am,an) =  Dis(am,an)*0.01+ Dis(am,an);
%     end
% end
Dis = dist2(X,X);%计算两两样本之间的欧氏距离
%% 2
H = zeros(n,n);
E = zeros(n,q+1);
for j = 1:n
 Dis_j = Dis(j,:)'; % 取第j个特征
 Dis_j_s = sort(Dis_j); % 对Dis_j进行排序
 [~,a,b] = intersect(Dis_j,Dis_j_s);% 计算Dis_j和Dis_j_s的交集,得到a的下标
 E(j,:) = a(1:q+1,:);
 H(a(1:q+1),j) = 1;
end
    %%寻找交集
    eqqqq=1;
    Ev=zeros(n,n);
    for i = 1:1:n %遍历每个E
     eq=1;
      for ii = 1:1:n                         %寻找第i个样本中所有关联点的圈-----step1 遍历每个点
        if H(ii,i)==1                        %寻找第i个样本中所有关联点的圈-----step2 寻找第i个样本中每个等于1的点所对应的，即以每个点ii为中心，画圈
              Eq(:,eq)=H(:,ii);              %寻找第i个样本中所有关联点的圈-----step3 进入到以ii为圆心的圈，找到所有关联点放入Eq中
              eq=eq+1;                       %第i个样本中第ii个圈遍历完继续向下遍历，eq+1
        end
      end                                
      Eq_sum = sum (Eq,2);                   % 寻找超边---step1 对行向量求和
      [~,Eq_n] = size(Eq);
      for iiii=1:1:n                         % 寻找超边---step2 每一行判断是否等于eq-1即每一行加起来为行数
       if Eq_sum(iiii,1) ==  Eq_n            % 寻找超边---step3 如果求和Eq中第iiii个点等于Eq的行数，认为第iiii行所有列都为1即有关联，这点则放入超边Ev中
         Ev(iiii,eqqqq)=1;                   % 寻找超边---step4  创建第q个超边
       end
      end
      if sum(Ev(:,eqqqq),1) > 0              %第i个样本如果形成超边了  则i放入超边中，超边数+1
          Ev(i,eqqqq)=1;                         
          eqqqq=eqqqq+1;
      end
    end
    Ev1=zeros(n,eqqqq-1);
    Ev1=Ev(:,1:eqqqq-1);
%% 3
W = zeros(n,1);
for i = 1:n
    x = X(i,:);
    e = E(i,:);
    w = 0;
    w_m = 0;
    for j = 1:(q+1)
       x_e = X(e(j),:); 
       w_m =w_m+EuDist2(x,x_e);
    end
    w_m = w_m/(q+1);
    for j = 1:(q+1)
       x_e = X(e(j),:); 
        w = w+exp(-EuDist2(x,x_e)/(w_m^2));
%         w = w+exp(-dist2(x,x_e)/ss);
    end
    W(i,1) = w;
end
%% 4 计算顶点度矩阵
De=diag(sum(Ev1));
W = diag(W);
Dv=diag(sum((W*Ev1)'));
% NewW=H*W*De*H';
% NewW=max(NewW,NewW');
% Dv_I=1./diag(Dv);
% Dv_I=diag(Dv_I);
% Lap=sqrt(Dv_I)*NewW*sqrt(Dv_I);
 Lap = Dv-W*Ev1*(De^(-1))*Ev1';
end