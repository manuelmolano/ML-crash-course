function kmeans_lesson_main2

% Now lets see what to do when non of the cluster validation methods
% coincide. This will most probably happen with data with high dimensions,
% and low amount of samples.

% Data comes from resting-state fMRI in the mouse brain. We have two
% datasets, and for each dataset, two independent runs of the kmeans
% algorithm (k = 2:20). Here, 'correlation' is the distance measure; 500
% iterations per replicate, 15 replicates (choosing the best).

% Data comes as follows: every subject has a series of fMRI 3D images taken
% every 1.2 seconds (40 subjects in the first dataset DS1 with 500
% timepoints, 41 subjects in DS2 with 300 timepoints). Each image (sample)
% consists of around 80.000 in-brain voxels (3D version of pixels) For each
% run of kmeans for each dataset, I concatenated the timeseries and
% clustered them according to their spatial similarity (correlation
% distance). Samples in each cluster are averaged into so-called
% 'CO-ACTIVATION PATTERNS' (CAPs).


%Lets load the results from 



clear, clc
%Results directories
main =  'D:\Gutierrez\mice_dFC\scripts_dFC_v5\utilities\kmeansClass';
cd(main)
Ndatasets =2;


Ncaps=2:20;

dataset1_path{1} = 'D:\Gutierrez\mice_dFC\scripts_dFC_v5\utilities\kmeansClass\DS1r1';
dataset1_path{2} = 'D:\Gutierrez\mice_dFC\scripts_dFC_v5\utilities\kmeansClass\DS1r2';

Nsubs = 40;
%% 1. Organize the information from each dataset.
 
% load dataset1 results
res1 = load([dataset1_path{1} '\cap_results.mat']); caps_results11 = res1.cap_results; clear res1
res1 = load([dataset1_path{2} '\cap_results.mat']); caps_results12 = res1.cap_results; clear res1

%% 2 compute CVI and plot it as an initial observation of cluster validity
for i = 1:length(Ncaps)
    relative_distortion1(i) = caps_results11.(['CAPS_' num2str(Ncaps(i))]).relative_distortion;
    relative_distortion2(i) = caps_results12.(['CAPS_' num2str(Ncaps(i))]).relative_distortion;
end


figure
plot(Ncaps(1:end),relative_distortion1,'o-b'),
hold on
plot(Ncaps(1:end),relative_distortion2,'*-b')
ylim([0.9*min([relative_distortion1 relative_distortion2]) 1.1*max([relative_distortion1 relative_distortion2])])
xlim([0 (max(Ncaps)+1)])
title('Cluster validity index - DS 1')

legend('r1','r2')

xlabel('Number of clusters')
ylabel('CVI (AU)')
xlim([0 (max(Ncaps)+1)])


%% 3. For each K number of clusters, compute the matrix of distances between CAPs belonging to each dataset.
for k = 1:length(Ncaps)
        
    tmp1 = caps_results11.(['CAPS_' num2str(Ncaps(k))]);
    tmp2 = caps_results12.(['CAPS_' num2str(Ncaps(k))]);
    [x_av1{k}, x_ind1{k}] = sort(tmp1.cap_occurrence,'descend');

    
    % first crop the maps to 
        
    Dmat{k} = zeros(Ncaps(k));

        for c = 1:Ncaps(k)
            
            for i = 1:Ncaps(k)
                for j = 1:Ncaps(k)
                    
                    Dmat{k}(i,j) = pdist([(spm_vec(tmp1.cap_mean_map(x_ind1{k}(i),:))'); (spm_vec(tmp2.cap_mean_map(j,:))')]);
                    
                end
            end
            
            
        end
end

%% 4. Assign CAPs from dataset1 to CAPs of dataset2.

for k = 1:length(Ncaps)
            
    [assignment{k}, cost{k}] = munkres_HA(Dmat{k});
            
end

%% 5. Compute inter-dataset similarities.

ds_sim = [];
for k = 1:length(Ncaps)
    
    tmp1 = caps_results11.(['CAPS_' num2str(Ncaps(k))]);
    tmp2 = caps_results12.(['CAPS_' num2str(Ncaps(k))]);
            
            ds_sim.cap_map_corr{k} = zeros(Ncaps(k),1);
            ds_sim.cap_map_mae{k} = zeros(Ncaps(k),1);
            ds_sim.cap_map_sae{k} = zeros(Ncaps(k),1);
            ds_sim.cap_map_mse{k} = zeros(Ncaps(k),1);
            ds_sim.cap_map_norm{k} = zeros(Ncaps(k),1);
            
            for c = 1:Ncaps(k)
                    
                    ds_sim.cap_map_corr{k}(c) = corr(spm_vec(tmp1.cap_mean_map(x_ind1{k}(c),:)), spm_vec(tmp2.cap_mean_map(assignment{k}(c),:)));
                    ds_sim.cap_map_mae{k}(c) = mae(spm_vec(tmp1.cap_mean_map(x_ind1{k}(c),:)) - spm_vec(tmp2.cap_mean_map(assignment{k}(c),:)));
                    ds_sim.cap_map_sae{k}(c) = sae(spm_vec(tmp1.cap_mean_map(x_ind1{k}(c),:)) - spm_vec(tmp2.cap_mean_map(assignment{k}(c),:)));
                    ds_sim.cap_map_mse{k}(c) = mse(spm_vec(tmp1.cap_mean_map(x_ind1{k}(c),:)) - spm_vec(tmp2.cap_mean_map(assignment{k}(c),:)));
                    ds_sim.cap_map_norm{k}(c) = norm(spm_vec(tmp1.cap_mean_map(x_ind1{k}(c),:)) - spm_vec(tmp2.cap_mean_map(assignment{k}(c),:)));
                    ds_sim.cap_map_cos{k}(c) = dot(spm_vec(tmp1.cap_mean_map(x_ind1{k}(c),:)),spm_vec(tmp2.cap_mean_map(assignment{k}(c),:)))/(norm(spm_vec(tmp1.cap_mean_map(x_ind1{k}(c),:)),2)*norm(spm_vec(tmp2.cap_mean_map(assignment{k}(c),:)),2));
            end
                        
            
end

%  Plot a boxplot for each K, with the values of spatial similarity between dataset CAPs.
            
C1 = []; C2 = [];C3 = [];C4 = [];C5 = []; C6 = [];
grp1 = []; grp2 = []; grp3 = []; grp4 = []; grp5 = []; grp6 = [];




labels = [];
for k = 1:length(Ncaps)
    
    for c = 1:Ncaps(k)
     cap_name{k}{c} = ['CAP' ' ' num2str(c)];
    end
    
    C1   = [C1 ds_sim.cap_map_corr{k}'];
    C2   = [C2 ds_sim.cap_map_mae{k}'];
    C3   = [C3 ds_sim.cap_map_sae{k}'];
    C4   = [C4 ds_sim.cap_map_mse{k}'];
    C5   = [C5 ds_sim.cap_map_norm{k}'];
    C6   = [C6 ds_sim.cap_map_norm{k}'];
    grp1 = [grp1 k*ones(1,Ncaps(k))];
    grp2 = [grp2 k*ones(1,Ncaps(k))];
    grp3 = [grp3 k*ones(1,Ncaps(k))];
    grp4 = [grp4 k*ones(1,Ncaps(k))];
    grp5 = [grp5 k*ones(1,Ncaps(k))];
    grp6 = [grp6 k*ones(1,Ncaps(k))];
    labels{1,k} = Ncaps(k);
end

figure(1)
boxplot(C1,grp1,'Whisker',1,'Labels',labels)
      title('Between Dataset Cluster Correlation')
      xlabel('Number of clusters')
      ylabel('Spatial Correlation')
      ylim([0 1])




end