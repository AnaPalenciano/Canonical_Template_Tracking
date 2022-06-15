%% CTT_template:
%  This script is a template to run canonical template tracking analysis 
%  on fMRI data [add biorxiv doi]. 
%  It uses functions from The Decoding Toolbox (TDT) by Martin Hebart & 
%  Kai Goergen (2015). The CTT functions provided (transres_CTT_*.m) can 
%  be integrated in the toolbox work flow. 
%  
%  This template carries out a simple CTT analysis on simulated fMRI data 
%  (as if it was processed with the software SPM).
%  The activity patterns are estimated from the GLM beta values. 
%% 

%% Set TDT defaults
cfg = decoding_defaults;
% Set additional TDT configuration parameters:
cfg.plot_design=0; 
cfg.results.overwrite = 1; % Overwrite previous results
% cfg.verbose = 2; % you want all information to be printed on screen
% cfg.plot_selected_voxels = 0; % 0: no plotting, 1: every step, 2: every second step, 100: every hundredth step...

%% Set the analysis that should be performed ('searchlight', 'roi')
cfg.analysis = 'roi';
% add additional analyses parameters. E.g.:
% cfg.analysis = 'searchlight
% cfg.searchlight.unit = 'mm';
% cfg.searchlight.radius = 12; 
% cfg.searchlight.spherical = 1;

%% Set the filename of your brain mask (searchlight analysis) or ROI mask:
cfg.files.mask = 'SimulatedData/sub1/anat/mask.nii';
%% Set the output directory where data will be saved
cfg.results.dir = 'SimulatedData/sub1/results_CTT/';

%% 1. Main paradigm data:

%  1.1. Set the filepath where your SPM.mat and all related betas are:
beta_dir_task = 'SimulatedData/sub1/GLM_Task/';


% %  1.2. Set the conditions names to extract the activity patterns. E.g.:
label_name{1} = '*cat*';
label_name{2} = '*dog*';
label_name{3} = '*hammer*';
label_name{4} = '*saw*';
%  And the conditions labels 
labels = 1:length(label_name);

% The following function extracts all beta names and corresponding run
% numbers from the SPM.mat
%regressor_names = design_from_spm(beta_dir_task);
load([beta_dir_task 'regressor_names.mat']);

% 1.3. Create a cfg structure with main paradigm data:
% Duplicate cfg structure, and populate it with the main task information 
cfg_task = cfg;
% Extract all information for the cfg.files structure 
cfg_task = decoding_describe_data(...
    cfg_task,label_name,labels,regressor_names,beta_dir_task);
% Add a variable with the source of the data (i.e.: main paradigm)
cfg_task.files.source = ones(size(cfg_task.files.label)); % 1 = main task
clear regressor_names
%% 2. Localizer(s) task data. 

%%  LOCALIZER 1:

%  2.1. Set the file path where your SPM.mat and all related betas are:
beta_dir_loc1 =  'SimulatedData/sub1/GLM_Loc1/';

%  2.2. Set the conditions names to extract the CANONICAL TEMPLATES:
labelname_loc1{1} = '*cat*';
labelname_loc1{2} = '*dog*';
labelname_loc1{3} = '*hammer*';
labelname_loc1{4} = '*saw*';

%  And the conditions labels. 
%  IMPORTANT: The localizer label should BE THE SAME as in the MAIN TASK 
labels_loc1 = 1:length(labelname_loc1); % i.e: "cat" = 1 in both datasets.

% The following function extracts all beta names and corresponding run
% numbers from the SPM.mat
%regressor_names_loc1 = design_from_spm(beta_dir_loc1);
load([beta_dir_loc1 'regressor_names.mat']);

% 2.3. Create a cfg structure with localizer data:
% Duplicate cfg structure, and populate it with the loc 1 information 
cfg_loc1 = cfg;

% Extract all information for the cfg.files structure 
cfg_loc1 = decoding_describe_data(...
    cfg_loc1,labelname_loc1,labels_loc1,regressor_names,beta_dir_loc1);
% Add a variable with the source of the data (i.e.: main paradigm)
cfg_loc1.files.source = -1*ones(size(cfg_loc1.files.label)); % -1 : loc 1
clear regressor_names
%% LOCALIZER 2 (if it applies)

%  2.1. Set the file path where your SPM.mat and all related betas are:
beta_dir_loc2 =  'SimulatedData/sub1/GLM_Loc2/';

%  2.2. Set the conditions names to extract the CANONICAL TEMPLATES:
labelname_loc2{1} = '*cat*';
labelname_loc2{2} = '*dog*';
labelname_loc2{3} = '*hammer*';
labelname_loc2{4} = '*saw*';

%  And the conditions labels. 
%  IMPORTANT: The localizer label should BE THE SAME as in the MAIN TASK 
labels_loc2 = 1:length(labelname_loc2); % i.e: "cat" = 1 in both datasets.

% The following function extracts all beta names and corresponding run
% numbers from the SPM.mat
%regressor_names_loc2 = design_from_spm(beta_dir_loc2);
load([beta_dir_loc2 'regressor_names.mat']);

% 2.3. Create a cfg structure with localizer data:
% Duplicate cfg structure, and populate it with the loc 2 information 
cfg_loc2 = cfg;

% Extract all information for the cfg.files structure 
cfg_loc2 = decoding_describe_data(...
    cfg_loc2,labelname_loc2,labels_loc2,regressor_names,beta_dir_loc2);
% Add a variable with the source of the data (i.e.: localizer 2)
cfg_loc2.files.source = -2*ones(size(cfg_loc2.files.label)); % -2 : loc 2
clear regressor_names
%% 3. Collapse task + localizers in same cfg structure for the CTT analysis

cfg.files.name  =  ...
    [cfg_task.files.name;  cfg_loc1.files.name; cfg_loc2.files.name]; 
cfg.files.label =  ...
    [cfg_task.files.label; cfg_loc1.files.label; cfg_loc2.files.label];
cfg.files.descr =  ...
    [cfg_task.files.descr cfg_loc1.files.descr cfg_loc2.files.descr];
cfg.files.chunk =  ...
    [cfg_task.files.chunk; cfg_loc1.files.chunk; cfg_loc2.files.chunk];
cfg.files.source = ...
    [cfg_task.files.source; cfg_loc1.files.source; cfg_loc2.files.source];
cfg.files.set = [];
cfg.files.xclass = [];

%% 4. Run CTT:

% We use a similarity design, but we modify it and split the main paradigm 
% and localizer(s) data into training and test datasets.
% This way, the TDT will automatically generate a correlation matrix 
% with rows = main task, columns = canonical templates.
cfg.design = make_design_similarity(cfg);
cfg.design.train_eq_test = 0;
cfg.design.train(cfg.files.source<0) = 0;
cfg.design.test(cfg.files.source>0) = 0;

cfg.decoding.software = 'similarity';
cfg.decoding.method = 'classification';

% The similarity metric can be modify here. For further details, please check: pattern_similarity.m
% Here we use the Fisher-z-transformed correlation similarity:
cfg.decoding.train.classification.model_parameters = 'zcorr'; 

%% CTT implementation:
% Option 1: pearson correlation (to compare templates from two conditions: task-relevant vs. task-irrelevant)
% Option 2: multiple regresssion (to compare templates from two localizers)
cfg.results.output = {'CTT_corr_relevant', 'CTT_corr_irrelevant', ...
    'CTT_MultReg_Loc1', 'CTT_MultReg_Loc2'};
cfg.design.unbalanced_data = 'ok';
results = decoding(cfg);

disp('Localizer 1 - Correlation-based distance')
disp(['Correlation between Main Task and Relevant Category template: ' num2str(results.CTT_corr_relevant.output{1})]);
disp(['Correlation between Main Task and Irrelevant Category template: ' num2str(results.CTT_corr_irrelevant.output{1})]);
disp('Localizers 1 and 2 - Multiple regresssion')
disp(['Beta value for localizer 1 templates: ' num2str(results.CTT_MultReg_Loc1.output{1})]);
disp(['Beta value for localizer 2 templates: ' num2str(results.CTT_MultReg_Loc2.output{1})]);





