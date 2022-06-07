%% Simulate patterns
% This script simulates activiy activity patterns for a RSA demonstration.
% The simulated data follows the paradigm from Haxby et al. (2001)
% Here we simulate simple clustered activity patterns following the 
%%
% Set up environment 
clear;clc;close all
rng('default')
casa = '/Data/CTT/code/SimulatedData';
addpath(genpath(casa));
%% Set up simulation parameters
% number of subjects simulated
nSubjects = 1;
% Experimental conditons: 
nCond = 4; 
cond = {'cat', 'dog', 'hammer', 'saw'};
% Number of betas per condition 
nBetas = 10;
nBetasLoc = 1;

% number of voxels/vertices in the ROI
nVoxels = [100 100]; 
% variance (power) of noise in neural patterns 
neural_noise=10; 

% Initialize simulationOptions variable
simulationOptions = simulationOptions_demo;

%% TASK ACTIVITY PATTERNS
% Specification of clusters of representations for generating patterns 
% 1st level: General pattern spread
clusterSpec_task{1} = neural_noise; 
% 2nd level: Stimulus category
clusterSpec_task{2}{1} = neural_noise; 
clusterSpec_task{2}{2}{1} = neural_noise;          % 2.1   Animate
clusterSpec_task{2}{2}{2} = {neural_noise,nBetas}; % 2.1.1 Cat
clusterSpec_task{2}{2}{3} = {neural_noise,nBetas}; % 2.1.2 Dog
clusterSpec_task{2}{3}{1} = neural_noise;          % 2.2   Inanimate
clusterSpec_task{2}{3}{2} = {neural_noise,nBetas}; % 2.2.1 Saw
clusterSpec_task{2}{3}{3} = {neural_noise,nBetas}; % 2.2.2 Hammer
  
%% lOCALIZER ACTIVITY PATTERNS 
% Specification of clusters of representations for generating patterns 
% 1st level: General pattern spread
clusterSpec_loc{1} = neural_noise; 
% 2nd level: Stimulus category
clusterSpec_loc{2}{1} = neural_noise; 
clusterSpec_loc{2}{2}{1} = neural_noise;             % 2.1   Animate
clusterSpec_loc{2}{2}{2} = {neural_noise,nBetasLoc}; % 2.1.1 Cat
clusterSpec_loc{2}{2}{3} = {neural_noise,nBetasLoc}; % 2.1.2 Dog
clusterSpec_loc{2}{3}{1} = neural_noise;             % 2.2   Inanimate
clusterSpec_loc{2}{3}{2} = {neural_noise,nBetasLoc}; % 2.2.1 Saw
clusterSpec_loc{2}{3}{3} = {neural_noise,nBetasLoc}; % 2.2.2 Hammer    


% Add whole-brain info
simulationOptions.clusterSpec = clusterSpec_task;
simulationOptions.brainVol = [84, 84, 50];
simulationOptions.effectCen = [40, 40, 20];
simulationOptions.nConditions = 40;

simulationOptions_loc = simulationOptions;
simulationOptions_loc.clusterSpec = clusterSpec_loc;
simulationOptions_loc.nConditions = 4;
%% Loop over subjects
for subject=1:nSubjects  
    % Create directory tree
    if not(isfolder(['sub' num2str(subject)]))
        mkdir(['sub' num2str(subject) '/GLM_Task']);
        mkdir(['sub' num2str(subject) '/GLM_Loc1']);
        mkdir(['sub' num2str(subject) '/GLM_Loc2']);
    end
    %% Generate MAIN PARADIGM betas
    betas = rsa.sim.simulateClusteredfMRIData_fullBrain(simulationOptions);
    betas = reshape(betas, [simulationOptions.nConditions, simulationOptions.brainVol]);
    % save betas:
    load("empty_hdr.mat");
    cd(['sub' num2str(subject) '/GLM_Task']);
    for b = 1:size(betas,1)
        hdr_tmp = hdr; 
        hdr_tmp.fname = ['beta_000' num2str(b) '.nii'];
        hdr_tmp.private.dat.fname = ['beta_000' num2str(b) '.nii'];
        spm_write_vol(hdr_tmp,squeeze(betas(b,:,:,:)));
        clear hdr_tmp;
    end
    regressor_names = cell(3,size(betas,1));
    cnt = 1;
    for c = 1:nCond
        for b = 1:nBetas    
            regressor_names{1,cnt} = cond{c};
            regressor_names{2,cnt} = 1;
            regressor_names{3,cnt} = ['Sn(' num2str(b) ') ' cond{c} '*bf(1)'];
            cnt = cnt + 1; 
        end
    end
    save('regressor_names.mat','regressor_names'); 
    clear regressor_names;
    cd ..
    cd ..
    %% Generate LOC 1 betas
    betas_loc1 = nan([nCond,simulationOptions_loc.brainVol]);
    for c = 1:nCond
        i = 1+(c-1)*10;
        ii = c*10;
        betas_loc1(c,:,:,:) = mean(betas(i:ii,:,:,:)); clear i ii
    end
    % save betas:
    load("empty_hdr.mat");
    cd(['sub' num2str(subject) '/GLM_Loc1']);
    for b = 1:size(betas_loc1,1)
        hdr_tmp = hdr; 
        hdr_tmp.fname = ['beta_000' num2str(b) '.nii'];
        hdr_tmp.private.dat.fname = ['beta_000' num2str(b) '.nii'];
        spm_write_vol(hdr_tmp,squeeze(betas_loc1(b,:,:,:)));
        clear hdr_tmp;
    end
    regressor_names = cell(3,size(betas_loc1,1));
    cnt = 1;
    for c = 1:nCond
        for b = 1:nBetasLoc    
            regressor_names{1,cnt} = cond{c};
            regressor_names{2,cnt} = 1;
            regressor_names{3,cnt} = ['Sn(' num2str(b) ') ' cond{c} '*bf(1)'];
            cnt = cnt + 1;
        end
    end
    save('regressor_names.mat','regressor_names'); 
    clear regressor_names;
    cd ..
    cd .. 
    %% Generate LOC 2 betas
    betas_loc2 = rsa.sim.simulateClusteredfMRIData_fullBrain(simulationOptions_loc);
    betas_loc2 = reshape(betas_loc2, [nCond, simulationOptions.brainVol]);
    % save betas:
    load("empty_hdr.mat");
    cd(['sub' num2str(subject) '/GLM_Loc2']);
    for b = 1:size(betas_loc2,1)
        hdr_tmp = hdr; 
        hdr_tmp.fname = ['beta_000' num2str(b) '.nii'];
        hdr_tmp.private.dat.fname = ['beta_000' num2str(b) '.nii'];
        spm_write_vol(hdr_tmp,squeeze(betas_loc2(b,:,:,:)));
        clear hdr_tmp;
    end
    regressor_names = cell(3,size(betas_loc2,1));
    cnt = 1;
    for c = 1:nCond
        for b = 1:nBetasLoc    
            regressor_names{1,cnt} = cond{c};
            regressor_names{2,cnt} = 1;
            regressor_names{3,cnt} = ['Sn(' num2str(b) ') ' cond{c} '*bf(1)'];
            cnt = cnt + 1;
        end
    end
    save('regressor_names.mat','regressor_names'); 
    clear regressor_names;
    cd ..
end

