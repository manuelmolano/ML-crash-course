function kmeans_lesson_main

% sorry for the long code:)

%% 1. Matlab example using the Iris dataset.
% the iris dataset is one of the most know in pattern recognition
% literature, and it's used broadly. It was collected by R.A. Fisher in
% 1936, and consists of 150 observations (samples)of 4 attributes (sepal
% width and length, petal width and length) of Iris flowers (Setosa,
% Versicolour, and Virginica).


%% 1. Basic kmeans and data representation.
%Let's start by observing 2 features and trying to cluster the subset of
%the data (with only these 2 features)

load fisheriris
X = meas(:,3:4);

figure;
plot(X(:,1),X(:,2),'k*','MarkerSize',5);
title 'Fisher''s Iris Data';
xlabel 'Petal Lengths (cm)';
ylabel 'Petal Widths (cm)';

% We would like, by using kmeans clustering, to identify and discriminate
% between the species of flowers given the samples of the mentioned
% features. Our ground truth is 3 types of flowers, and kmeans will try to
% cluster the data into groups given the multivariate data similarity
[idx,C] = kmeans(X,3);

figure;
plot(X(idx==1,1),X(idx==1,2),'r.','MarkerSize',12)
hold on
plot(X(idx==2,1),X(idx==2,2),'b.','MarkerSize',12)
hold on
plot(X(idx==3,1),X(idx==3,2),'g.','MarkerSize',12)
plot(C(:,1),C(:,2),'kx','MarkerSize',15,'LineWidth',3)
legend('Cluster 1','Cluster 2','Cluster 3','Centroids','Location','NW')
title 'Cluster Assignments and Centroids'
hold off

%% 2. kmeans options and K selection.

% Cluster assessment is done by measuring the distance between each sample
% and its cluster centroid, which can be chosen (euclidean, cityblock,
% cosine, correlation, etc.).

% Most of the times you do not know the amount of clusters (ground truth),
% which makes the analysis a bit harder. Kmeans can be run for different
% "k's" and you can assess which solution (amount of clusters) best
% classifies data. There are A LOT of methods to do this. Also, if you
% think your data could be better clustered, replicate the algorithm a
% couple of times and automatically choose the solution with the least sum
% of distances.
Nk = 10;

for k = 1:Nk
    [idx2{k},C2{k},sumd2{k},D2{k}] = kmeans(X,k,'MaxIter',5000,'Distance','cityblock','Display','final','Replicates',10);
end

% here "idx" is the index of the cluster assignation of the sample. "C" is
% the centroid location in units of the data. "sumd" is the sum of
% point-to-centroid distances within a cluster. "D" is the distance matrix
% from each point to each cluster-centroid.

% Assessment of solutions. 2.1. "Elbow method". One of the many Cluster
% Validity Indexes (CVI). Check how each solution has clusters with data
% close to its centroid, as compared with how distanced are cluster
% centroids.
sd = zeros(Nk,1); btw_cdist = zeros(Nk,1); CVI = zeros(Nk,1);

for k = 2:Nk
    sd(k) = sum(sumd2{k});
    btw_cdist(k) = sum(pdist(C2{k}));
    CVI(k) = sd(k)/btw_cdist(k);
end
plot(2:Nk, CVI(2:Nk))

% 2.2 Silhouette cluster evaluation. The silhouette value for the ith
% point, Si, is defined as 
% Si = (bi-ai)/ max(ai,bi)
% where ai is the average distance from the ith point to the other points in the same cluster as i,
% and bi is the minimum average distance from the ith point to points in a
% different cluster, minimized over clusters.

figure
[silh3,h] = silhouette(X,idx2{3},'cityblock');
h = gca;
h.Children.EdgeColor = [.8 .8 1];
xlabel 'Silhouette Value'
ylabel 'Cluster'

% another option is to evaluated different k's.
eva_silh = evalclusters(X,'kmeans','Silhouette','klist',1:10);
plot(eva_silh)
  
% 2.3 Davis-Boulin index. The Davies-Bouldin criterion is based on a ratio of
% within-cluster and between-cluster distances. The Davies-Bouldin index is
% defined as within-to-between cluster distance ratio for the ith and jth
% clusters.
eva_DB = evalclusters(X,'kmeans','DaviesBouldin','klist',1:10);
plot(eva_DB)



% Lets try generating some data and one last cluster evaluation criteria
mu1 = [2 2];
sigma1 = [0.9 -0.0255; -0.0255 0.9];

mu2 = [5 5];
sigma2 = [0.5 0 ; 0 0.3];

mu3 = [-2, -2];
sigma3 = [1 0 ; 0 0.9];

N = 200;

X2 = [mvnrnd(mu1,sigma1,N);...
     mvnrnd(mu2,sigma2,N);...
     mvnrnd(mu3,sigma3,N)];
 
 gscatter(X2(:,1),X2(:,2))
 
% Evaluate the performance of kmeans with the gAP stat
%2.4 GAP statistic. Is said to be the best for high dimensional data and non-gaussian data. TAKES
%A LOT OF TIME!!!
% The gap criterion formalizes the elbow approach by estimating the “elbow”
% location as the number of clusters with the largest gap value. Therefore,
% under the gap criterion, the optimal number of clusters occurs at the
% solution with the largest local or global gap value within a tolerance
% range. GAP is the expected value of the difference between the real CVI
% (sum of distances), and the CVI if the data were generated from uniform
% distribution.

eva_gap = evalclusters(X2,'kmeans','gap','KList',[1:10]);

eva_gap_iris = evalclusters(X,'kmeans','gap','KList',[1:10]);

% THE RANDOM INITIALIZATION OF CENTROIDS IS NOT A PROBLEM ANYMORE (SINCE
% 2007), AS MATLAB BY DEFAULT (UNLESS SPECIFIED OTHERWISE) USES THE
% KMEANS++ ALGORITHM, WHICH INITIALIZES CENTROIDS THAT IN A POSITION WITH
% HIGH PROBABILITY OF BEING DISTANT FROM THE NEXT CLUSTER.

















end

