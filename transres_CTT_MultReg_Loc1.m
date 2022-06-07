% function output = transres_CTT_MultReg_Loc1(decoding_out,~,cfg,data)
% 
% This function computes a multiple regression, predicting the main task 
% activity patterns with a series of canonical templates
% 
% To use it: 
% cfg.results.output = {'transres_CTT_MultReg'}
%
% A.F.Palenciano, 06/06/2022

function output = transres_CTT_MultReg_Loc1(decoding_out,~,cfg,data)

% Extract conditions from main paradigm
labels_task = cfg.design.label(cfg.files.source>0);
% Extract conditions from localizer task  
labels_loc  = cfg.design.label(cfg.files.source<0);
source_loc  = cfg.files.source(cfg.files.source<0);

% Regressors: CANONICAL TEMPLATES + constant term
X = [data(cfg.files.source<0,:);ones(size(data(1,:)))]';
% Predicted data: MAIN PARADIGM
Y = data(cfg.files.source>0,:)';

% Allocate beta variable: 
beta_mean = nan(size(labels_task));

% Predict each paradigm activity pattern with canonical templates
for i = 1:size(Y,2) 
    betas = regress(Y(:,i),X);
    beta_mean(i)  = mean(betas(and(...
        labels_loc==labels_task(i), source_loc == -1)));
    clear betas_tmp;
end
% return the averaged beta values per each canonical template
output = {mean(beta_mean)};