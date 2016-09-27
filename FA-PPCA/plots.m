function plots(Xsim,W,Zsim)

figure();
scatter(Xsim(:,1),Xsim(:,2),'black','filled');
hold on
mu = mean(Xsim);
scatter(mu(1),mu(2),'MarkerFaceColor','green','MarkerEdgeColor','green');
hold on

x = -5:15;
m = W(2,1)/W(1,1);
y = m*(x-mu(1)) +mu(2);
plot(x,y,'black');
hold on

for i = 1:1:size(Zsim,1)
    temp = W'.*Zsim(i);
    datacomp(i,:) = temp(:) + mu(:);
    plot(datacomp(i,1),datacomp(i,2),'o','MarkerFaceColor','red','MarkerEdgeColor','red');
    hold on
    A = [Xsim(i,1) datacomp(i,1)];
    B = [Xsim(i,2) datacomp(i,2)];
    plot(A(:),B(:),'red');
    hold on
end   
axis equal;
title('Projections in PC Space');
xlabel('Dimension 1');
ylabel('Dimension 2');