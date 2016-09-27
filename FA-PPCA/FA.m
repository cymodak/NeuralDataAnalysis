close all
clear
clc

if(~isdeployed)
    cd(fileparts(which(mfilename)));
end

load ps8_data.mat

q = 1;

d = size(Xsim,2);
N = size(Xsim,1);

for i = 1:1:size(Xsim,1)
    cent_Xsim(i,:) = Xsim(i,:) - mean(Xsim,1);
end
[U,S,V] = svd(cent_Xsim);

W = V(:,1);
psi = S(1:d,1:d);

mu = mean(Xsim);

L = [];
threshold = 1;
while(threshold > 0)
    M = inv(W'*inv(psi)*W + eye(q));
    Zsim3 = zeros(q,N);
    E2 = zeros(q,q,N);
    for i = 1:1:N
        Zsim3(:,i) = M*W'*inv(psi)*(Xsim(i,:)-mu)';
        E2(:,:,i) = M + Zsim3(:,i)*Zsim3(:,i)';
    end
    
    T1 = 0;
    T2 = 0;
    T3 = 0;
    T4 = 0;
    T5 = 0;
    for i = 1:1:N
        T2 = T2 + 1/2*trace(E2(:,:,i));
        T3 = T3 + 1/2*(Xsim(i,:)-mu)*inv(psi)*(Xsim(i,:)-mu)';
        T4 = T4 - (Xsim(i,:)-mu)*inv(psi)*W*Zsim3(:,i);
        T5 = T5 + 1/2*trace(W'*inv(psi)*W*E2(:,:,i));
    end
    temp2 = T2 + T3 + T4 + T5 + N/2*log(det(psi));
    
    L = [L -temp2];
    
    T6 = 0;
    T7 = 0;
    for i = 1:1:N
        T6 = T6 + (Xsim(i,:)-mu)'*Zsim3(:,i)';
        T7 = T7 + E2(:,:,i);
    end
    
    W = T6*inv(T7);
    
    T8 = 0;
    T9 = 0;
    for i = 1:1:N
        T8 = T8 + (Xsim(i,:)-mu)'*(Xsim(i,:)-mu);
        T9 = T9 - W*Zsim3(:,i)*(Xsim(i,:)-mu);
    end
    psi = 1/N*diag(diag(T8 + T9));
    
    T1 = 0;
    T2 = 0;
    T3 = 0;
    T4 = 0;
    T5 = 0;
    for i = 1:1:N
        T2 = T2 + 1/2*trace(E2(:,:,i));
        T3 = T3 + 1/2*(Xsim(i,:)-mu)*inv(psi)*(Xsim(i,:)-mu)';
        T4 = T4 - (Xsim(i,:)-mu)*inv(psi)*W*Zsim3(:,i);
        T5 = T5 + 1/2*trace(W'*inv(psi)*W*E2(:,:,i));
    end
    temp2 = T2 + T3 + T4 + T5 + N/2*log(det(psi));
    
    L = [L -temp2];
    threshold = norm(L(end)-L(end-1))/norm(L(end-1));
    
end

for i = 1:1:size(Xsim,1)
    cent_Xsim(i,:) = Xsim(i,:) - mean(Xsim,1);
end

figure()
plot(L);
title('Data Log-Likelihood');
xlabel('Iterations');
ylabel('Log Likelihood');

plots(Xsim,W,Zsim3');
