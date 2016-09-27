close all
clear 
clc

if(~isdeployed)
    cd(fileparts(which(mfilename)));
end

load ps6_data.mat
mu1 = InitParams1.mu;
sigma1 = InitParams1.Sigma;
pi1 = InitParams1.pi;
data = Spikes;

for i=1:1:size(mu1,2)
    sigma(:,:,i) = sigma1;
end
sigma1 = sigma;

count = 0;
q = [];
while(count<100)
    gamma = expectation(data,mu1,pi1,sigma1);
    obj = model(data,mu1,sigma1,gamma,pi1);
    q = [q obj];
    [pi1,mu1,sigma1] = maximization(data,gamma);
    count = count + 1;
end

figure();
plot(q);
title('Data Log Likehilood');
ylabel('Log Likelihood');
xlabel('Iterations');

class = gamma - 0.5;
class = sign(class);
class = (1+class)/2;
plots(mu1,sigma1,class,data);

