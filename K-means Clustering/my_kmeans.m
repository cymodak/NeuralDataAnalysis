function [clust,cent,total_distance] = my_kmeans(D, K, start)

% Number of points
len = length(D);
% Distance matrix
dist = zeros(len,1);
init = randi(len,1);
% Matrix to store centroid values
cent = start';

% K-means algorithm
clust = zeros(len,1);
% Matrix for end of iterations
differ = ones(size(cent));
total_distance = [];

while (max(max(abs(differ)>0)))
    % Assigning point to cluster
    count = 1;
    while (count <= K)
        for i = 1:1:count
            for j = 1:1:len
                temp = (D(j,:)-cent(i,:))*(D(j,:)-cent(i,:))';
                if (temp < dist(j) || count == 1)
                    dist(j) = temp;
                    clust(j) = count;
                end
            end
        end
        count = count + 1;
    end
    
    % Computing variance
    temp = 0;
    for i = 1:1:len            
            temp = temp + (D(i,:)-cent(clust(i),:))*(D(i,:)-cent(clust(i),:))';
    end
    total_distance = [total_distance temp];
    
    % Computing new means
    for i = 1:1:K
        sum1 = zeros(1,size(D,2));
        num = 0;
        for j = 1:1:len
            if (clust(j) == i)
                sum1 = sum1 + D(j,:);
                num = num + 1;
            end
        end
        differ(i,:) = sum1/num - cent(i,:);
        cent(i,:) = sum1/num;        
    end   
    % Computing variance
    temp = 0;
    for i = 1:1:len            
            temp = temp + (D(i,:)-cent(clust(i),:))*(D(i,:)-cent(clust(i),:))';
    end
    total_distance = [total_distance temp];
end

end
