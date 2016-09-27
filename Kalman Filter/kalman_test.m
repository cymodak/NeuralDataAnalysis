function [Z_pred,Error] = kalman_test(A,Q,pie,V,C,R,X)

T = size(X,2);
siz = size(A,1);

mu_final = zeros(siz,T);
mu_inter = zeros(siz,T);
sigma_final = zeros(siz,siz,T);
sigma_inter = zeros(siz,siz,T);

mu_inter(:,1) = pie;
sigma_inter(:,:,1) = V;

% for t = 1:1:T
%     K = sigma_inter(:,:,t)*C'*inv(C*sigma_inter(:,:,t)*C'+R);
%     mu_final(:,t) = mu_inter(:,t) + K*(X(:,t)-C*mu_inter(:,t));
%     sigma_final(:,:,t) = sigma_inter(:,:,t) - K*C*sigma_inter(:,:,t);
%     
%     if(t == T)
%         break
%     end
%     
%     mu_inter(:,t+1) = A*mu_final(:,t);
%     sigma_inter(:,:,t+1) = A*sigma_final(:,:,t)*A' + Q;
% end

for t = 1:1:T
    if (t == 1)
        mu_inter(:,t) = A*mu_inter(:,1);
        sigma_inter(:,:,t) = A*sigma_inter(:,:,t)*A' + Q;
    else
        mu_inter(:,t) = A*mu_final(:,t-1);
        sigma_inter(:,:,t) = A*sigma_final(:,:,t-1)*A' + Q;
    end   
    K = sigma_inter(:,:,t)*C'*inv(C*sigma_inter(:,:,t)*C'+R);
    mu_final(:,t) = mu_inter(:,t) + K*(X(:,t)-C*mu_inter(:,t));
    sigma_final(:,:,t) = sigma_inter(:,:,t) - K*C*sigma_inter(:,:,t);    
end

Z_pred = mu_final;
Error = sigma_final;