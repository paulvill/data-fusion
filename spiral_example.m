addpath_embtime2d
presets
%%

% parameters
a = 1;
b = 0.1;
c1 = 1;
d1 = 0.005;
e1 = 65;


% labeled data points
l1 = 120;
indl1 = 1:l1;

% unlabeled data points
l2 = 300;
indl2 = l1+1:l1+l2;

% number of points
n = l1 + l2;

% matrix containing all the modalities
M = zeros(3,n);

% noise in the various channels
sigma = [0.02,0.02,0];

% time - for each data point
t1 = (1:99/(l2-1):100)';
t = (1:99/(999):100)';
t2 = t(randi(length(t),l1,1));
tsort = sort([t1;t2]);
indl2_t = find(ismember(tsort,t1));
[~ , indl1_temp] = sort(t2);
indl1_t = find(ismember(tsort,sort(t2)));

% common modality
M(1,:) = [a*([t2;t1]).*(cos(b*([t2;t1]))+sigma(1)*randn(n,1))];
M(2,:) = [a*([t2;t1]).*(sin(b*([t2;t1]))+sigma(2)*randn(n,1))];

% labels
M(3,indl1) = c1*t2.*(exp(-d1*(t2-e1).^2)+sigma(3)*randn(l1,1));

% visualizing (x1, x2)
cmap = gray(l1+l2); %winter(l1+l2);
cmap = cmap(end:-1:1,:);

figure,
scatter(M(1,indl2),M(2,indl2),50,cmap(indl2_t,:),'filled');
for i = 1:length(indl2_t)-1,
    hold on,
    plot([M(1,indl2(i)),M(1,indl2(i+1))],[M(2,indl2(i)),M(2,indl2(i+1))],'-','Color',cmap(indl2_t(i+1),:),'LineWidth',2);
end
xlabel('x_1');
ylabel('x_2');
set(gca,'fontsize', 24);

% visualizing (x1, x2, t)
cmap2 = colormap(brewermap(l1+l2,'reds'));
col2 = [255,98,58]/255;
close

h = figure,
scatter3(M(1,indl2),M(2,indl2),zeros(l2,1),50,cmap(indl2_t,:),'filled');
for i = 1:length(indl2_t)-1,
    hold on,
    plot([M(1,indl2(i)),M(1,indl2(i+1))],[M(2,indl2(i)),M(2,indl2(i+1))],'-','Color',cmap(indl2_t(i+1),:),'LineWidth',2);
end
for j = 1:length(indl1),
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
d = date;
mkdir(d);
saveas(h,[d,'/x1_x2_labels'],'png');

% computing pairwise differences between pairs of (x1(i), x2(i)) and (x1(j), x2(j))
dist = distances(M(1:2,:));
W = AffinityFromDistance(dist,10);
D = diag(sum(W, 2));

% solving the semi-supervised problem
[T, Tu, Tl] = transformation(indl1,indl2);
[inv_u] = ssl_estimate4(W, D, T, l1);
fu = inv_u*M(3,indl1)';


col3 = [70,105,220]/255;

% visualizing the predicted and original labels in the (x1, x2, t) space
h = figure,
scatter3(M(1,indl2),M(2,indl2),zeros(l2,1),50,'w','filled');
for j = 1:length(indl1),
    i = indl1_temp(j);
    ytemp = 0:0.5:(M(3,indl1(i)));
    xtemp = [M(1,indl1(i)).*ones(length(ytemp),1),M(2,indl1(i)).*ones(length(ytemp),1)];
    scatter3(M(1,indl1(i)),M(2,indl1(i)),M(3,indl1(i)),50,col2,'filled');
    hold on,
end
hold on,
scatter3(M(1,indl2),M(2,indl2),fu,50,cmap(indl2_t,:),'filled');
for i = 1:length(indl2_t)-1,
    hold on,
    plot3([M(1,indl2(i)),M(1,indl2(i+1))],[M(2,indl2(i)),M(2,indl2(i+1))],[fu(i),fu(i+1)],'-','Color',cmap(indl2_t(i+1),:),'LineWidth',4);
end

xlabel('x^1');
ylabel('x^2');
zlabel('y');
set(gca,'fontsize', 24);
grid off

% figure is saved
d = date;
mkdir(d);
saveas(h,[d,'/x1_x2_predicted_labels'],'png');

% visualizing the original and predicted labels as a function of time
h = figure,
hold on,
scatter(t1,fu,50,cmap(indl2_t,:),'filled');
for i = 1:length(indl2_t)-1,
    hold on,
    plot([t1(i),t1(i+1)],[fu(i),fu(i+1)],'-','Color',cmap(indl2_t(i+1),:),'LineWidth',4);
end
hold on,
for j = 1:length(indl1),
    i = indl1_temp(j);

   scatter(t2(i),M(3,indl1(i)),50,col2,'filled');
hold on,
end
xlabel('t');
ylabel('y');
set(gca,'fontsize', 24);
grid off

% figure is saved
d = date;
mkdir(d);
saveas(h,[d,'/labels_t'],'png');

