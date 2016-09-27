close all
clear
clc

if(~isdeployed)
    cd(fileparts(which(mfilename)));
end

load ps7_data.mat

mu_init = InitParams.mu;
sigma_init = InitParams.Sigma;
data = Spikes;

figure();
plot(data);

for i = 1:1:size(data,2)
    cent_data(:,i) = data(:,i) - mean(data,2);
end

[U,S,V] = svd(cent_data);
figure();
plot(-U(:,1),'red');
hold on
plot(U(:,2),'green');
hold on
plot(U(:,3),'blue');

score = sum(S);

figure();
plot(score(1:31),'o');

D1 = U(:,1);
D2 = -U(:,2);
D = [D1 D2];
data2 = data'*D;
scatter(data2(:,1),data2(:,2));
K = size(mu_init,2);

V = 4;
siz = floor(size(data,2)/V);

set = cell(1,V);
result = zeros(K,V);

for i=1:1:V
    set{1,i} = data2((i-1)*siz+1:i*siz,:);
end

for i = 1:1:K
    array = (1:V)';    
    for m = 1:1:V
        train = [];
        prob = 1/i*ones(1,i);
        mu = mu_init(:,1:i);
        
        for j = 1:1:i
            sigma(:,:,j) = sigma_init;
        end
        
        for j=1:1:V-1
            train = [train ; set{1,array(j)}];
        end
        train = train';
        test = set{1,array(V)};
        test = test';
        
        count = 0;       
        while(count<10)
            gamma = expectation(train,mu,prob,sigma);
            [prob,mu,sigma] = maximization(train,gamma);
            count = count + 1;
        end
        obj = model(test,mu,sigma,gamma,prob);
        result(i,m) = obj;
        array = circshift(array,1);
    end
    if (i == 3)
        muhat = mu';
        dat = D*muhat';
    end
end

figure();
plot(dat);


