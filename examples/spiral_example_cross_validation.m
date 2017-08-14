%%% SPIRAL EXAMPLE CROSS VALIDATION %%%

% This script uses the same setup as `spiral_example` to illustrate the harmonic
% extension method on a toy example dataset. Here, we examine the variability in
% the reconstruction by cross-validating over K = 5 folds. We also vary the
% number of unlabeled data points (that is points containing only x^1 and x^2,
% not the labels y) and to capture its effect on the reconstruction results.

% Load visualization presets.
presets;

%% Parameters

% Parameters for the synthetic trajectory to be generated.
a = 1;
b = 0.1;
c1 = 1;
d1 = 0.005;
e1 = 65;

% Noise in the various channels.
sigma = [0.02, 0.02, 0];

% Number of labeled data points.
l1 = 120;

% Number of unlabeled data points.
l2 = 300;

% Number of nearest neighbors to use when determining scale for affinity matrix.
num_neighbors = 10;

% Number of cross-validation bins.
K = 5;

% Number of repetitions for each configuration of labeled/unlabeled data
% points.
nrep = 100;

%% Generate data

% labeled data points index
indl1 = 1:l1;

% unlabeled data points index
indl2 = l1+1:l1+l2;

% total number of points
n = l1 + l2;

% matrix containing the common modality (x^1, x^2) and the labels y
M = zeros(3,n);

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

%% K-fold cross validation

% number of unlabeled data points
n_unlbds = (0:10:l2);


% vector where the absolute error is stored
abs_error = zeros(length(n_unlbds),nrep);

for u = 1:length(n_unlbds) 
    fprintf([num2str(n_unlbds(u)),' \n']);
    for r = 1:nrep
        % defining the n_unlbds(u) subsamples of the unlabeled data points
        nb_unlbld = n_unlbds(u);
        indu = indl2;
        indl = indl1;
        
        % extracting randomly nb_unlbld unlabeled data points from the
        % original set of unlabeled data points
        mask_unlbld = randperm(length(indu),nb_unlbld);
        
        % the randomly extracted unlabeled data points are reordered and
        % shifted according to the number of labeled data points
        indu_mask = length(indl)+sort(mask_unlbld);
        % the combined set of labeled and unlabeled data points
        mask_tot = [indl,indu_mask];

        % computing the affinity matrix on the reduced number of data
        % points
        dist = distances(M(1:2,:));
        W = AffinityFromDistance(dist(mask_tot,mask_tot), num_neighbors); 
        
        % defining randomly K subsamples of the labeled data points
        nb_labels = length(indl);
        
        ind_l_sub = zeros(nb_labels/K,K);
        mask = 1:nb_labels;
        
        for k = 1:K
            ind_temp = randperm(nb_labels - (k-1)*nb_labels/K,nb_labels/K);
            ind_l_sub(:,k) = mask(ind_temp);
            mask(ind_temp) = [];
        end
        
        % predicting label value on the on the selected labeled points
        predicted_labels = zeros(nb_labels,1);
        for k = 1:K
            indu = (1:nb_unlbld);
            indl = (1:nb_labels);
            
            % the set of unlabeled data points corresponds to the original
            % unlabeled data points and the labeled data points from the
            % current bin
            indu = [(ind_l_sub(:,k)'),indu+nb_labels];
            
            % the set of labeled data points corresponds to the original
            % labeled data points minus the labeled data points from the
            % current bin
            indl(ind_l_sub(:,k)) = [];

            l1 = length(indl);
            l2 = length(indu);
            
            % the semi-supervised learning problem is solved accordingly
            [inv_u] = ssl_estimate( W,indl,indu);
            
            % we obtain the label prediction
            fu = inv_u*M(3,indl)';
            % which is stored in a vector
            predicted_labels(ind_l_sub(:,k)) = fu(1:length(ind_l_sub(:,k)));
    
        end
        % the absolute error is computed by comparing the predicted labels
        % to the original labels
        abs_error(u,r) = (1/length(predicted_labels))*sum(abs(predicted_labels(:) - M(3,1:nb_labels)'))/(max(M(3,1:nb_labels)) - min(M(3,1:nb_labels)));
        fprintf([num2str(abs_error(u,r)),' \n']);
    end
end


m = mean(abs_error,2);
s = std(abs_error,0,2);

h = figure; 
p = errorbar(n_unlbds,m,s);
set(p,'LineWidth',5,'Color',[0.5 0.5 0.5]);
hold on,
plot(n_unlbds,m,'o','Color',[0.5 0.5 0.5],'MarkerSize',10,'MarkerFaceColor',[0.5 0.5 0.5])
xlabel('Number of unlabeled points');
ylabel('Normalized Absolute Error');

d = 'output';
mkdir(d);
save(fullfile(d, 'abs_error_spiral.mat'),'abs_error');
saveas(h,fullfile(d, 'abs_error_spiral'),'png');
