function [pi_k,mu,sigma] = maximization(data,gamma)

N = sum(gamma);
pi_k = N/sum(N);
k = size(gamma,2);
l = size(data,2);
% disp(N);
temp = data*gamma;

for i=1:1:k
    mu(:,i) = temp(:,i)/N(i);
end

for i=1:1:k
    temp2 = 0;
    for j=1:1:l
        temp1 = data(:,j)-mu(:,i);
        temp2 = temp2 + gamma(j,i)*(temp1*temp1');
    end
    sigma(:,:,i) = temp2/N(:,i);   
end

end