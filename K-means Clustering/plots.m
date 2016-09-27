function [] = plots(D,clust,cent,obj)
Fs = 30;
l = size(cent,1);
t = 1/Fs*(1:1:size(cent,2));
for j = 1:1:l
    figure();
    plot(t,250*ones(size(D,1)));
    hold on
    [r, c] = find(clust == j);
    for i = 1:1:length(r)
        plot(t, D(:,r(i)));
        hold on
    end
    p = plot(t,cent(j,:),'red');
    set (p, 'LineWidth', 4);
    xlabel('Time (ms)');
    ylabel('Amplitude (uV)');
    title('Clustered Spikes');
end

figure();
t = 1:1:length(obj);
plot(t/2,obj);
hold on
for i = 1:1:length(t)
    if(rem(i,2) == 1)
        scatter(t(i)/2,obj(i),'o','red')
    else
        scatter(t(i)/2,obj(i),'x','green');
    end
end
title('Objective function');
ylabel('Error (Distance Function)');
xlabel('Iterations');