% compare Rhiir model to other published fungal models
clear
clc

modelFile = 'model/iRi1572.mat';
load(modelFile);

base_path = 'results/fungal-models';
ec_path = 'data/corrected-EC-numbers.csv';
ec_translation_tab = readtable(ec_path, 'ReadVariableNames', false);

% EC number sets of fungal models
ec_tab = readtable(fullfile(base_path, 'model_comparison_ec.csv'),...
    'ReadVariableNames', true, 'Delimiter', '\t');
n_comp = size(ec_tab,2);
ec_classes = 7;

% correct EC numbers (update)
disp('correcting EC numbers')

% EC number set of R. irregularis
ec_rhiir = correctEC(model.rxnECNumbers, ec_translation_tab);
ec_rhiir = regexp(ec_rhiir, '\|', 'split');
ec_rhiir = unique([ec_rhiir{:}]');

ec_translation_tab = readtable(ec_path, 'ReadVariableNames', false);
for i=1:n_comp
    ec_tmp = correctEC(ec_tab.(i), ec_translation_tab);
    ec_tmp =  regexp(ec_tmp, '\|', 'split');
    ec_tmp = [ec_tmp{:}]';
    
    if numel(ec_tmp) > size(ec_tab,1)
        append = numel(ec_tmp) - size(ec_tab,1);
        ec_tab = [ec_tab;repmat({''}, append, n_comp)];
    end
end; clear ec_tmp append


% loop over the EC number sets of other fungal models and calculate the
% Jaccard distance
jd = ones(n_comp,ec_classes);
n_ec = zeros(n_comp,ec_classes);
n_ec_rhiir = zeros(1, ec_classes);
i_ec = zeros(n_comp,ec_classes);

disp('calculating Jaccard distances of EC number sets')
for i=1:n_comp
    fprintf('\t> %s\n', ec_tab.Properties.VariableNames{i});
    for j=1:ec_classes
        % if column is not empty
        if numel(unique(ec_tab.(i)))>1
            
            ec_tmp = ec_tab.(i);
            ec_tmp = ec_tmp(~cellfun('isempty', ec_tmp));
            
            % only compare the current EC class
            ec_tmp = regexp(ec_tmp, ['^', num2str(j), '.*'], 'match');
            ec_tmp = [ec_tmp{:}]';
            
            ec_model = regexp(ec_rhiir, ['^', num2str(j), '.*'], 'match');
            ec_model = [ec_model{:}]';
            
            jd(i,j) = 1 - numel(intersect(ec_model, ec_tmp)) / ...
                numel(union(ec_model, ec_tmp));
            i_ec(i,j) = numel(intersect(ec_model, ec_tmp));
            n_ec(i,j) = numel(ec_tmp);
            n_ec_rhiir(j) = numel(ec_model);
        end
    end
end; clear ec_tmp ec_model

writetable(array2table(jd, 'RowNames', ec_tab.Properties.VariableNames,...
    'VariableNames', cellstr(strcat('EC',num2str((1:ec_classes)')))), fullfile(base_path, 'jd_ec_class.txt'),...
    'WriteVariableNames', true, 'WriteRowNames', true, 'Delimiter', '\t')