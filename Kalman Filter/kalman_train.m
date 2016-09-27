function [A,Q,pie,V,C,R] = kalman_train(Z_train,spikecount_train)

numtrials = size(Z_train,1);
numangles = size(Z_train,2);

A1 = 0;
A2 = 0;
totallength = 0;
for i = 1:1:numtrials
    for j = 1:1:numangles
        temp = Z_train{i,j};
        l = size(temp,2);
        totallength = totallength + l; 
        for k=2:1:l
            A1 = A1 + temp(:,k)*temp(:,k-1)';
            A2 = A2 + temp(:,k-1)*temp(:,k-1)';
        end
    end
end
A = (A1/totallength)*inv(A2/totallength);

Q = 0;
for i = 1:1:numtrials
    for j = 1:1:numangles
        temp = Z_train{i,j}; 
        l = size(temp,2);
        for k=2:1:l
            Q1 = temp(:,k) - A*temp(:,k-1);
            Q = Q + Q1*Q1';
        end
    end
end
Q = Q/(totallength-(numtrials*numangles));

C1 = 0;
C2 = 0;
pie = 0;
for i = 1:1:numtrials
    for j = 1:1:numangles
        temp1 = Z_train{i,j};
        temp2 = spikecount_train{i,j};
        pie = pie + temp1(:,1);
        l = size(temp1,2);
        for k=1:1:l
            C1 = C1 + temp2(:,k)*temp1(:,k)'; 
            C2 = C2 + temp1(:,k)*temp1(:,k)';            
        end
    end
end
C = (C1/totallength)*inv(C2/totallength);
pie = pie/totallength;

R = 0;
V = 0;
for i = 1:1:numtrials
    for j = 1:1:numangles
        temp1 = Z_train{i,j};
        temp2 = spikecount_train{i,j};
        V1 = temp1(:,1)-pie;
        V = V + V1*V1';
        l = size(temp1,2);
        for k=1:1:l
            R1 = temp2(:,k)-C*temp1(:,k);
            R = R + R1*R1';
        end
    end
end
R = R/totallength;
V = V/totallength;


