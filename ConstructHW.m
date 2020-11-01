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
Dis = dist2(X,X);%������������֮���ŷ�Ͼ���
%% 2
H = zeros(n,n);
E = zeros(n,q+1);
for j = 1:n
 Dis_j = Dis(j,:)'; % ȡ��j������
 Dis_j_s = sort(Dis_j); % ��Dis_j��������
 [~,a,b] = intersect(Dis_j,Dis_j_s);% ����Dis_j��Dis_j_s�Ľ���,�õ�a���±�
 E(j,:) = a(1:q+1,:);
 H(a(1:q+1),j) = 1;
end
    %%Ѱ�ҽ���
    eqqqq=1;
    Ev=zeros(n,n);
    for i = 1:1:n %����ÿ��E
     eq=1;
      for ii = 1:1:n                         %Ѱ�ҵ�i�����������й������Ȧ-----step1 ����ÿ����
        if H(ii,i)==1                        %Ѱ�ҵ�i�����������й������Ȧ-----step2 Ѱ�ҵ�i��������ÿ������1�ĵ�����Ӧ�ģ�����ÿ����iiΪ���ģ���Ȧ
              Eq(:,eq)=H(:,ii);              %Ѱ�ҵ�i�����������й������Ȧ-----step3 ���뵽��iiΪԲ�ĵ�Ȧ���ҵ����й��������Eq��
              eq=eq+1;                       %��i�������е�ii��Ȧ������������±�����eq+1
        end
      end                                
      Eq_sum = sum (Eq,2);                   % Ѱ�ҳ���---step1 �����������
      [~,Eq_n] = size(Eq);
      for iiii=1:1:n                         % Ѱ�ҳ���---step2 ÿһ���ж��Ƿ����eq-1��ÿһ�м�����Ϊ����
       if Eq_sum(iiii,1) ==  Eq_n            % Ѱ�ҳ���---step3 ������Eq�е�iiii�������Eq����������Ϊ��iiii�������ж�Ϊ1���й������������볬��Ev��
         Ev(iiii,eqqqq)=1;                   % Ѱ�ҳ���---step4  ������q������
       end
      end
      if sum(Ev(:,eqqqq),1) > 0              %��i����������γɳ�����  ��i���볬���У�������+1
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
%% 4 ���㶥��Ⱦ���
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