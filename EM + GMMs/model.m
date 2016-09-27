function [objective] = model(data,mu,sigma,gamma,pi_k)

k = size(gamma,2);
n = size(gamma,1);
siz = size(data,1);

objective = 0;
for i=1:1:n
    sum = 0;
    for j=1:1:k
%         sum = sum + pi_k(j)*(2*pi)^(-siz/2)*det(sigma(:,:,j))^(-1/2)*exp(-(1/2)*(data(:,i)-mu(:,j))'*inv(sigma(:,:,j))*(data(:,i)-mu(:,j)));
          sum = sum + pi_k(j)*mvnpdf(data(:,i),mu(:,j),sigma(:,:,j));
    end
    objective = objective + log(sum);
end

end