function []= plots(mu,sigma,class,data,clust)

k = size(mu,2);
n = size(class,1);
t = find(class);
for i = 1:1:k
    figure();
    for j = 1:1:n
        if(t(j)> (i-1)*n)
            if(t(j) > i*n)
                break;
            else
                ind = t(j)-(i-1)*n;
                plot(data(:,ind));
                hold on
            end
        end
    end
    p = plot(mu(:,i),'red');
    set (p, 'LineWidth', 3);
    hold on
    p = plot(mu(:,i)-sqrt(diag(sigma(:,:,i))),'--r');
    set (p, 'LineWidth', 3);
    hold on
    p = plot(mu(:,i)+sqrt(diag(sigma(:,:,i))),'--r');
    set (p, 'LineWidth', 3);
    xlabel('Spike Interval (ms)');
    ylabel('Amplitude (uV)');
    title('Clustered Spikes');
end




    