% function output = transres_CTT_MultReg_Loc2(decoding_out,~,cfg,data)
% 
% This function fits a multiple regression, predicting the main task 
% activity patterns with a series of canonical templates
%
% Instead of using decoding_out.opt (task - localizer(s) correlation 
% matrix), this function directly uses the variable data, containing
% the activity patterns (i.e., before the correlation) 
% 
% To use it: 
% cfg.results.output = {'transres_CTT_MultReg_Loc1','transres_CTT_MultReg_Loc2'}
%
% A.F.Palenciano, 06/06/2022

function output = transres_CTT_MultReg_Loc2(decoding_out,~,cfg,data)
%% Extract conditions' labels from the main paradigm 
labels_task = cfg.design.label(cfg.files.source>0);
%% Extract conditions' labels from the localizer task
%  Here we will use data from two localizers (labelled 
%  as -1, -2 in cfg.files.source)
labels_loc  = cfg.design.label(cfg.files.source<0);
%% Keep track of the identity of the localizer task
source_loc  = cfg.files.source(cfg.files.source<0);

%% Regressors (X): CANONICAL TEMPLATES + constant term 
X = [data(cfg.files.source<0,:);ones(size(data(1,:)))]';
%% Predicted data (Y): MAIN PARADIGM activity patterns
Y = data(cfg.files.source>0,:)';

% Pre-allocate the variable with beta weights: 
beta_mean = nan(size(labels_task));

%% MULTIPLE REGRESSION. 
%  The regression will be performed iteratively. In each
%  iteration, we'll try to predict an activity pattern from 
%  the main task using the canonical templates from both 
%  localizer tasks.
for i = 1:size(Y,2) 
    betas = regress(Y(:,i),X);
    % Store the beta weight from Localizer 2 task-relevant conditions:
    beta_mean(i)  = mean(betas(and(...
        labels_loc==labels_task(i), source_loc == -2)));
    clear betas_tmp;
end
% return the averaged task-relevant beta values for Localizer 2
output = {mean(beta_mean)};
