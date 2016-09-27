close all
clear
clc

if(~isdeployed)
    cd(fileparts(which(mfilename)));
end

load ps8_data.mat
data = Xplan;
numangles = 8;

for i = 1:1:size(data,1)
    cent_data(i,:) = data(i,:) - mean(data,1);
end

[U,S,V] = svd(cent_data);
[X Y Z] = eig(cent_data*cent_data');
D = U(:,1:3);
% data3 = D*S(1:3,1:3);
data3 = cent_data*V(:,1:3);
K = size(data3,1);
k = K/numangles;
UM = V(:,1:3);

figure()
for i = 1:1:numangles
    colour = de2bi(i-1,3);
    if (colour == [1 1 1])
        colour = [0.75 0.5 0.5];
    end
    for j = 1:1:k        
        scatter3(data3((i-1)*k+j,1),data3((i-1)*k+j,2),data3((i-1)*k+j,3),'MarkerFaceColor',colour,'MarkerEdgeColor',colour);
        hold on
    end
end   

title('Clustered data in 3-D Space');
xlabel('First PC direction');
ylabel('Second PC direction');
zlabel('Third PC direction');

score = sum(S);

figure();
plot(score(1:97),'o');
title('Squarerooted Eigenvalue Spectrum');
xlabel('Eigenvector');
ylabel('Eigenvalue');

figure();
imagesc(UM);
colorbar;
title('Contribution of Neurons to Principal Components');
xlabel('Eigenvector');
ylabel('Neuron');

for i = 1:1:size(Xsim,1)
    cent_Xsim(i,:) = Xsim(i,:) - mean(Xsim,1);
end
[U1,S1,V1] = svd(cent_Xsim);

Dsim = U1(:,1);
Zsim1 = Dsim*S1(1,1);

plots(Xsim,V1(:,1),Zsim1);
