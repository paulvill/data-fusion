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

%%
% K-fold cross validation
% number of bins
K = 5;

% number of unlabeled data points
n_unlbds = (0:10:300);

% number of repetitions for each configuration of labeled/unlabeled data
% points
nrep = 100;

% vector where the absolute error is stored
abs_error = zeros(length(n_unlbds),nrep);

for u = 1:length(n_unlbds), 
    fprintf([num2str(n_unlbds(u)),' \n']);
    
    for r = 1:nrep,
        % defining the n_unlbds(u) subsamples of the unlabeled data points
        nb_unlbld = n_unlbds(u);
        indu = indl2;
        indl = indl1;
        
        mask_unlbld = randperm(length(indu),nb_unlbld);
        
        indu_mask = length(indl)+sort(mask_unlbld);
        mask_tot = [indl,indu_mask];
        nunlbd = length(indu_mask);
        
        % K-Fold Cross-Validation
        dist = distances(M(1:2,:));
        W = AffinityFromDistance(dist(mask_tot,mask_tot),10);
        D = diag(sum(W, 2));
        
        % defining the K subsamples of the labeled data points
        nb_labels = length(indl);
        
        ind_l_sub = zeros(nb_labels/K,K);
        mask = 1:nb_labels;
        
        for k = 1:K,
            ind_temp = randperm(nb_labels - (k-1)*nb_labels/K,nb_labels/K);
            ind_l_sub(:,k) = mask(ind_temp);
            mask(ind_temp) = [];
        end
        
        % predict label value on the on the selected labeled points
        predicted_labels = zeros(nb_labels,1);
        for k = 1:K,
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
            [T, Tu, Tl] = transformation(indl,indu);
            [inv_u] = ssl_estimate4( W, D, T, l1);
            
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

h = figure, 
errorbar(n_unlbds,m,s,'LineWidth',5,'Color',[0.5 0.5 0.5])
hold on,
plot(n_unlbds,m,'o','Color',[0.5 0.5 0.5],'MarkerSize',10,'MarkerFaceColor',[0.5 0.5 0.5])
xlabel('Number of unlabeled points');
ylabel('Normalized Absolute Error');

d = date;
mkdir(d);
save([d,'/abs_error_spiral.mat'],'abs_error');
saveas(h,[d,'/abs_error_spiral'],'png');
