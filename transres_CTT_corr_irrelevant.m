% function output = transres_CTT_corr(decoding_out,~,cfg,data)
% 
% This function extracts the correlations among a series of canonical
% templates (from on localizer task) with the main paradigm activity 
% patterns paradidgm data (decoding_out.opt).
% 
% The correlation matrix (N betas from the main task * M templates from the
% localizer) is already generated in the main decoding.m script. 
% Here we just average the pertinent cells. 
% 
% To use it: 
% cfg.results.output = {'transres_CTT_corr_relevant', 'transres_CTT_corr_irrelevant'}
%
% A.F.Palenciano, 06/06/2022
%%

function output = transres_CTT_corr_irrelevant(decoding_out,~,cfg,data)

%% Extract conditions' labels from the main paradigm 
labels_task = cfg.design.label(cfg.files.source>0);
%% Extract conditions' labels from the localizer task  
%  Here we just use one localizer, labelled as -1 in cfg.files.source
%  By modifying this index, different localizers will be accessed. 
labels_loc  = cfg.design.label(cfg.files.source==-1); 

%% Find correspondece among task and localizer.
%  Here we will compare the correlation of task-relevant and
%  task-irrelevant templates. The rel variable codes the cells 
%  from the correlation matrix (decoding.opt) with matching labels
%  across task and localizer. 
rel = zeros(length(labels_task),length(labels_loc));
for i = 1:length(labels_loc)
    rel(:,i) = labels_task==labels_loc(i);
end

%% Extract correlation values from cells with non-matching conditions 
%  i.e.: extract correlation with task-irrelevant templates, used as baseline.
corr_irr = decoding_out.opt(rel==0);

% return output 
output = {mean(corr_irr)};
