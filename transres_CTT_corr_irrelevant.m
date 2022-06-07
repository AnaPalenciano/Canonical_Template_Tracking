% function output = transres_CTT_corr(decoding_out,~,cfg,data)
% 
% This function extracts the correlations among a series of canonical
% templates (from on localizer task) with the main paradigm activity 
% patterns paradidgm data (decoding_out.opt).
% 
% The correlation matrix (N betas from the main task * M templates from the
% localizer) is already generated. Here we just average the pertinent
% cells. 
% 
% To use it: 
% cfg.results.output = {'transres_CTT_corr_relevant', 'transres_CTT_corr_irrelevant'}
%
% A.F.Palenciano, 06/06/2022
%%

function output = transres_CTT_corr_irrelevant(decoding_out,~,cfg,data)

% Extract conditions from main paradigm
labels_task = cfg.design.label(cfg.files.source>0);
% Extract conditions from localizer task  
labels_loc  = cfg.design.label(cfg.files.source==-1);

% Find correspondece among tasks
rel = zeros(length(labels_task),length(labels_loc));
for i = 1:length(labels_loc)
    rel(:,i) = labels_task==labels_loc(i);
end

% In case they are used for the baseline, extract correlation values from
% cells with non-matching conditions
corr_irr = decoding_out.opt(rel==0);

% return the optional output only
output = {mean(corr_irr)};