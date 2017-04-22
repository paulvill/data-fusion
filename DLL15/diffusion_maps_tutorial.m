clear all
close all

%% Create Data Set
% First, we create an artifical data set |x| of $n=7$ data points which lie 
% on a one-dimensional curve in two dimensions.
% This data set can be parameterized by t he angular direction $\theta$.

n = 7;

theta = linspace(0, 3*pi/2, n);

x = [0.9*cos(theta); 0.9*sin(theta)]

figure;
scatter(x(1,:),x(2,:),200, 'k', '.');
xlabel('x(1)');
ylabel('x(2)');

%% Compute Distance Matrix
% Then we compute the $n \times n$ distance matrix |dist_matrix| between data points
% where |dist_matrix(i,j)| is the distance between data point $i$ and data
% point $j$.
% Assuming the data points are stored in columns, this can be done using
% the |dist| function.

dist_matrix = dist(x)

%% Calculate Kernel Matrix
% We then compute the $n \times n$ matrix $W$, with 
% $$ W_{ij} = \exp \left(-\frac{ \| x_i - x_j \|^2}{\epsilon^2} \right) $$
% where $\epsilon$, the kernel scale, is a characteristic distance within
% the data set.
% We typically set $\epsilon$ to be some fraction of the median of the
% pairwise distances. Here we take $\epsilon$ as half the median of the
% pairwise distances.

eps = median(dist_matrix(:))/2

W = exp(-dist_matrix.^2/eps.^2)


%% Compute Diffusion Matrix
% We then compute the $n \times n$ diaganol matrix $D$, with 
% $$ D_{ii} = \sum_{j=1}^n W_{ij} $$ and
% the $n \times n$ matrix $$A = D^{-1}W. $$

D = diag(sum(W))

A = inv(D)*W

%% Calculate Eigenvectors and Eigenvalues
% We then compute the top eigenvectors and eigenvalues of $A$, and order
% them by the magnitude of the eigenvalues.

% compute eigenvectors and eigenvalues
[evecs, evals] = eig(A);

% sort eigenvectors and eigenvalues
[~, I] = sort(abs(diag(evals)), 'descend');
evecs = evecs(:,I)
evals = evals(I,I)

%% Order Data by First (nontrivial) Diffusion Maps Eigenvector
% The first eigenvector is a trivial constant eigenvector with eigenvalue
% 1. The second eigenvector provides an ordering of the data. As expected, 
% we see The second eigenvector is one-to-one
% with the $\theta$, which parameterizes the arclength along the curve. 
% Furthermore, if we color the data by the value of this second eigenvector, 
% we can see that this eigenvector
% does order the data along the one-dimensional nonlinear curve. 

figure;
scatter(theta, evecs(:,2), 50, 'k')
xlabel('\theta (parameterizes arclength along curve)')
ylabel('first diffusion maps coordinate')

figure;
scatter(x(1,:),x(2,:),200, evecs(:,2), '.');
xlabel('x(1)');
ylabel('x(2)');
h = colorbar;
ylabel(h, 'first diffusion maps coordinate')
