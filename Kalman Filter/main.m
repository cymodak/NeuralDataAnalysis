close all
clear
clc

if(~isdeployed)
    cd(fileparts(which(mfilename)));
end

load ps9_data.mat

dt = 1;

numtrials = size(train_trial,1);
numangles = size(train_trial,2);
bin_size = 20;

spikecount_train = cell(numtrials,numangles);
spikecount_test = cell(numtrials,numangles);
handPos_train = cell(numtrials,numangles);
handPos_test = cell(numtrials,numangles);
handVel_train = cell(numtrials,numangles);
handVel_test = cell(numtrials,numangles);
Z_train = cell(numtrials,numangles);

for i = 1:1:numtrials
    for j = 1:1:numangles
        temp1 = train_trial(i,j).spikes;
        temp2 = train_trial(i,j).handPos;
        temp3 = test_trial(i,j).spikes;
        temp4 = test_trial(i,j).handPos;
        
        l1 = floor(dt*size(temp1,2)/bin_size);
        for k = 1:1:l1
            spikecount_train{i,j} = [spikecount_train{i,j} sum(temp1(:,(k-1)*bin_size/dt+1:k*bin_size/dt),2)];            
            handPos_train{i,j} = [handPos_train{i,j} dt*temp2(:,k*bin_size/dt)];
            if(k ~= 1)
                handVel_train{i,j} = [handVel_train{i,j} (handPos_train{i,j}(:,k)-handPos_train{i,j}(:,k-1))/bin_size];
            end
        end
        Z_train{i,j} = [handPos_train{i,j}(:,1:l1-1) ; handVel_train{i,j}];
        
        l2 = floor(dt*size(temp3,2)/bin_size);
        for k = 1:1:l2
            spikecount_test{i,j} = [spikecount_test{i,j} sum(temp3(:,(k-1)*bin_size/dt+1:k*bin_size/dt),2)];            
            handPos_test{i,j} = [handPos_test{i,j} dt*temp4(:,k*bin_size/dt)];
            if(k ~= 1)
                handVel_test{i,j} = [handVel_test{i,j} (handPos_test{i,j}(:,k)-handPos_test{i,j}(:,k-1))/bin_size];
            end
        end
        Z_test{i,j} = [handPos_test{i,j}(:,1:l2-1) ; handVel_test{i,j}];
    end
end

[A,Q,pie,V,C,R] = kalman_train(Z_train,spikecount_train);

test_sample = spikecount_test{1,1};
z_true = Z_test{1,1};
[z_pred,z_err] = kalman_test(A,Q,pie,V,C,R,test_sample(:,1:size(z_true,2)));
figure()
for i = 1:1:size(z_true,2)
    plots(z_pred(1:2,i),z_err(1:2,1:2,i));
    hold on
end
hold on
plot(z_pred(1,:),z_pred(2,:),'red','LineWidth',2);
hold on
plot(z_true(1,:),z_true(2,:),'black','LineWidth',2);

t1 = [];
t2 = [];
t3 = [];
t4 = [];
for i=1:1:size(z_true,2)
    t1 = [t1 z_pred(1,i)-sqrt(z_err(1,1,i))];
    t2 = [t2 z_pred(1,i)+sqrt(z_err(1,1,i))];
    t3 = [t3 z_pred(2,i)-sqrt(z_err(2,2,i))];
    t4 = [t4 z_pred(2,i) + sqrt(z_err(2,2,i))];
end
figure()
subplot(2,1,1);
plot(z_true(1,:),'black','LineWidth',2);
hold on
plot(z_pred(1,:),'red','LineWidth',2);
hold on
plot(t1,'--r','LineWidth',2);
hold on
plot(t2,'--r','LineWidth',2)
xlabel('Time (1 unit = 20ms)');
ylabel('Horizontal Position');
title('Predicted Position (red) and Actual Position (black)')
subplot(2,1,2);
plot(z_true(2,:),'black','LineWidth',2);
hold on
plot(z_pred(2,:),'red','LineWidth',2);
hold on
plot(t3,'--r','LineWidth',2);
hold on
plot(t4,'--r','LineWidth',2);
xlabel('Time (1 unit = 20ms)');
ylabel('Vertical Position');
title('Predicted Position (red) and Actual Position (black)')

test_sample = spikecount_test{1,4};
z_true = Z_test{1,4};
[z_pred,z_err] = kalman_test(A,Q,pie,V,C,R,test_sample(:,1:size(z_true,2)));
figure()
for i = 1:1:size(z_true,2)
    plots(z_pred(1:2,i),z_err(1:2,1:2,i));
    hold on
end
hold on
plot(z_pred(1,:),z_pred(2,:),'red','LineWidth',2);
hold on
plot(z_true(1,:),z_true(2,:),'black','LineWidth',2);

t1 = [];
t2 = [];
t3 = [];
t4 = [];
for i=1:1:size(z_true,2)
    t1 = [t1 z_pred(1,i)-sqrt(z_err(1,1,i))];
    t2 = [t2 z_pred(1,i)+sqrt(z_err(1,1,i))];
    t3 = [t3 z_pred(2,i)-sqrt(z_err(2,2,i))];
    t4 = [t4 z_pred(2,i) + sqrt(z_err(2,2,i))];
end
figure()
subplot(2,1,1);
plot(z_true(1,:),'black','LineWidth',2);
hold on
plot(z_pred(1,:),'red','LineWidth',2);
hold on
plot(t1,'--r','LineWidth',2);
hold on
plot(t2,'--r','LineWidth',2)
xlabel('Time (1 unit = 20ms)');
ylabel('Horizontal Position');
title('Predicted Position (red) and Actual Position (black)')
subplot(2,1,2);
plot(z_true(2,:),'black','LineWidth',2);
hold on
plot(z_pred(2,:),'red','LineWidth',2);
hold on
plot(t3,'--r','LineWidth',2);
hold on
plot(t4,'--r','LineWidth',2);
xlabel('Time (1 unit = 20ms)');
ylabel('Vertical Position');
title('Predicted Position (red) and Actual Position (black)')

distance = 0;
l = 0;
for i = 1:1:numtrials
    for j = 1:1:numangles
        test_sample = spikecount_test{i,j};
        z_true = Z_test{i,j};
        l = l + size(z_true,2);
        [z_pred,z_err] = kalman_test(A,Q,pie,V,C,R,test_sample(:,1:size(z_true,2)));
        z1 = z_true(1:2,:);
        z2 = z_pred(1:2,:);
        z3 = z1-z2;
        z4 = trace(sqrt(z3'*z3));
        distance = distance + z4;
    end
end
distance = distance/l;
