function [gamma] = expectation(data,mu,pi_k,sigma)

gamma = [];
l= size(data,2);
siz = size(data,1);
m = size(mu,2);

for i=1:1:l
    for j=1:1:m
%         p(i,j) = pi_k(j)*(2*pi)^(-siz/2)*det(sigma(:,:,j))^(-1/2)*exp(-(1/2)*(data(:,i)-mu(:,j))'*inv(sigma(:,:,j))*(data(:,i)-mu(:,j)));
          p(i,j) = pi_k(j)*mvnpdf(data(:,i),mu(:,j),sigma(:,:,j));
    end
    temp1 = log(p(i,:));
    temp2 = log(sum(p(i,:)));
    gamma(i,:) = exp(temp1-temp2);
end

end