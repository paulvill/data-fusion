addpath_datafusion

%%

% parameters value
a = 1;
b = 0.1;
c1 = 1;
d1 = 0.005;
e1 = 65;


% labeled data points index
l1 = 120;
indl1 = 1:l1;

% unlabeled data points index
l2 = 300;
indl2 = l1+1:l1+l2;

% total number of points
n = l1 + l2;

% matrix containing the common modality (x^1, x^2) and the labels y
M = zeros(3,n);

% noise in the various channels
sigma = [0.02,0.02,0];

% time step - for each data point
% the time steps of unlabeled data points are evenly distributed between 
% 1 and 100
t1 = (1:99/(l2-1):100)';
% the time steps of labeled data points are randomly drawn between 1 and
% 100
t = (1:99/(999):100)';
t2 = t(randi(length(t),l1,1));

% ordering the points according to their time steps
tsort = sort([t1;t2]);
% unlabeled data points
indl2_t = find(ismember(tsort,t1));
% labeled data points
[~ , indl1_temp] = sort(t2);

% computing common modality values based on time steps
M(1,:) = [a*([t2;t1]).*(cos(b*([t2;t1]))+sigma(1)*randn(n,1))];
M(2,:) = [a*([t2;t1]).*(sin(b*([t2;t1]))+sigma(2)*randn(n,1))];

% computing label values based on time steps
M(3,indl1) = c1*t2.*(exp(-d1*(t2-e1).^2)+sigma(3)*randn(l1,1));

% visualizing distribution of data points in the common modality (x^1, x^2)
cmap = gray(l1+l2);
cmap = cmap(end:-1:1,:);

figure,
scatter(M(1,indl2),M(2,indl2),50,cmap(indl2_t,:),'filled');
for i = 1:length(indl2_t)-1
    hold on,
    plot([M(1,indl2(i)),M(1,indl2(i+1))],[M(2,indl2(i)),M(2,indl2(i+1))],'-','Color',cmap(indl2_t(i+1),:),'LineWidth',2);
end
xlabel('x_1');
ylabel('x_2');
set(gca,'fontsize', 24);

% visualizing distribution of labeled and unlabeled data points in the 
% (x^1, x^2, y)-space
col2 = [255,98,58]/255;
close

h = figure;
% unlabeled data points
scatter3(M(1,indl2),M(2,indl2),zeros(l2,1),50,cmap(indl2_t,:),'filled');
for i = 1:length(indl2_t)-1
    hold on,
    plot([M(1,indl2(i)),M(1,indl2(i+1))],[M(2,indl2(i)),M(2,indl2(i+1))],'-','Color',cmap(indl2_t(i+1),:),'LineWidth',2);
end
% labeled data points
for j = 1:length(indl1)
    i = indl1_temp(j);
    ytemp = 0:0.5:(M(3,indl1(i)));
    xtemp = [M(1,indl1(i)).*ones(length(ytemp),1),M(2,indl1(i)).*ones(length(ytemp),1)];
    scatter3(M(1,indl1(i)),M(2,indl1(i)),M(3,indl1(i)),50,col2,'filled');
end
xlabel('x^1');
ylabel('x^2');
zlabel('y');
set(gca,'fontsize', 24);
grid off

% figure is saved
d = 'output';
mkdir(d);
saveas(h,fullfile(d, 'x1_x2_labels'),'png');

% computing pairwise differences between pairs of (x1(i), x2(i)) and 
% (x1(j), x2(j)) and corresponding affinity matrix W = (w_i,j)
dist = distances(M(1:2,:));
W = AffinityFromDistance(dist,10);

% solving the semi-supervised problem
[inv_u] = ssl_estimate(W,indl1,indl2);
fu = inv_u*M(3,indl1)';


col3 = [70,105,220]/255;

% visualizing the predicted and original labels in the (x1, x2, t) space
h = figure;
scatter3(M(1,indl2),M(2,indl2),zeros(l2,1),50,'w','filled');
for j = 1:length(indl1)
    i = indl1_temp(j);
    ytemp = 0:0.5:(M(3,indl1(i)));
    xtemp = [M(1,indl1(i)).*ones(length(ytemp),1),M(2,indl1(i)).*ones(length(ytemp),1)];
    scatter3(M(1,indl1(i)),M(2,indl1(i)),M(3,indl1(i)),50,col2,'filled');
    hold on,
end
hold on,
scatter3(M(1,indl2),M(2,indl2),fu,50,cmap(indl2_t,:),'filled');
for i = 1:length(indl2_t)-1
    hold on,
    plot3([M(1,indl2(i)),M(1,indl2(i+1))],[M(2,indl2(i)),M(2,indl2(i+1))],[fu(i),fu(i+1)],'-','Color',cmap(indl2_t(i+1),:),'LineWidth',4);
end

xlabel('x^1');
ylabel('x^2');
zlabel('y');
set(gca,'fontsize', 24);
grid off

% figure is saved
d = 'output';
mkdir(d);
saveas(h,fullfile(d, 'x1_x2_predicted_labels'),'png');

% visualizing the original and predicted labels as a function of time
h = figure;
hold on,
scatter(t1,fu,50,cmap(indl2_t,:),'filled');
for i = 1:length(indl2_t)-1
    hold on,
    plot([t1(i),t1(i+1)],[fu(i),fu(i+1)],'-','Color',cmap(indl2_t(i+1),:),'LineWidth',4);
end
hold on,
for j = 1:length(indl1)
    i = indl1_temp(j);
    
    scatter(t2(i),M(3,indl1(i)),50,col2,'filled');
    hold on,
end
xlabel('t');
ylabel('y');
set(gca,'fontsize', 24);
grid off

% figure is saved
d = 'output';
mkdir(d);
saveas(h,fullfile(d, 'labels_t'),'png');

