function [U,V,obj] = SHNMF_MCC(X, beta, alpha, maxiter, intau, intav, option ,L)

%%     This code solves the following problem SHNMF-MCC

    
 %L = D - W;   
 U=intau;
 V=intav;
%  Norm = 2;
%  NormV = 0;
% [U,V] = NormalizeUV(U, V, NormV, Norm);
 %%
%  D1 = X - U*V;
%  QQ=sqrt(sum( D1.*D1,1) + 1e-6);
%  QQA=1./QQ;
%   D2 = diag (QQA);
[m,n] = size(X);
for i=1:maxiter
    iter=maxiter;
    cgm_mn = sum(sum((X-U*V').^2,2),1);
    cgm = sqrt(cgm_mn/0.5*m); %求cgm
    e_g = sum((X-U*V').^2 ,2)/ (2*cgm^2); % X-UV'列求和平方，m*1向量，再除以2*cgm^2
    e_i = exp(-e_g);
    H_ii = e_i ;  %求-rho
    H = diag(H_ii);
 
%%                       
     %updating U  
      U=U.*((H*X*V)./(H*U*V'*V));
      U=max(U,1e-50);
%     U=U./(ones(dim(1),1)*sum(U));
     %updating V
     V=V.*((X'*H*U)./(V*U'*H*U+alpha*L*V +beta));
     DX=X-U*V';
      obj_1 = mean(mean(abs(DX)))/mean(mean(X));
      obj(i) =obj_1;
     if mod(i,20)==0 || i==iter
        if option.dis
            disp(['Iterating >>>>>> ', num2str(i),'th']);
        end
        curRes=obj(i) ;
        if option.residual>=curRes || i==option.iter   
        end
     end 
end


function [U, V] = NormalizeUV(U, V, NormV, Norm)% normv=0；
    K = size(U,2);
   
    if Norm == 2
        if NormV
            norms = max(1e-15,sqrt(sum(V.^2,1)))';% 2范数
            V= V*spdiags(norms.^-1,0,K,K);
            U = U*spdiags(norms,0,K,K);
        else
            norms = max(1e-15,sqrt(sum(U.^2,1)))';
            U = U*spdiags(norms.^-1,0,K,K);
           V = V*spdiags(norms,0,K,K);
         %  V = spdiags(norms,0,K,K)*V;
        end
    else
        if NormV
            norms = max(1e-15,sum(abs(V),1))';%1 范数
            V = V*spdiags(norms.^-1,0,K,K);
            U = U*spdiags(norms,0,K,K);
        else
            norms = max(1e-15,sum(abs(U),1))';
            U = U*spdiags(norms.^-1,0,K,K);
            V = V*spdiags(norms,0,K,K);
        end
    end

function obj = Hf(DX,H)
  obj_N = sum((DX.^2),2);
  obj_H = H*obj_N;
  obj = sum(obj_H,1);