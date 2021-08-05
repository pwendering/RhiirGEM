% calculate growth rates from hyphae dry weights,  dry weight of parent
% spores and growth duration
% DW(t) = DW(0) * exp(mu * t)
% <==> mu = ln( DW(t) / DW(0) ) / t
topDir = '';

% read growth predictions with FBA and eMOMENT with different carbon sources and
% concentration as in Hildebrandt et al. 2006
mu_emoment = readtable([topDir, 'results\carbon-sources\growth.csv'],...
    'ReadVariableNames', false, 'ReadRowNames', true);
carbonSourceNames = mu_emoment.Properties.RowNames;
concentrations = [10 100 1000];
mu_emoment = table2array(mu_emoment);

mu_fba = readtable([topDir, 'results\carbon-sources\growth-fba.csv'],...
    'ReadVariableNames', false, 'ReadRowNames', true);
mu_fba = table2array(mu_fba);

% Total protein from Table 2 in Hildebrandt et al. 2006 (10.1111/j.1574-6968.2005.00027.x)
p_tot = [
    91	95	96;
    88	98	89;
    102	103	106;
    83	86	75;
]/1000;

% Final hyphal dry weight from Table 2 in Hildebrandt et al. 2006 (10.1111/j.1574-6968.2005.00027.x)
DW_t = [
    3.4	3.6	3.3;
    3.2	3.3	3.1;
    3.7	4	4.5;
    3.3	3.5	3
];

% Number of parent spores from Hildebrandt et al. 2006: started experiment from 1000 spores
n_spores = 1000;
% spore dry weight (Sugiura et al. 2020, 10.1073/pnas.2006948117) (Fig S1)
weight_per_spore = .3E-6;
DW_0 = weight_per_spore*n_spores;
% Growth duration from Hildebrandt et al. 2006: 2.5 month (transform to
% hours)
t = 2.5 * 30 * 24;
% calculate growth rate
mu = log(DW_t / DW_0) / t;

% calculate correlations:

[rho,pval] = corr(reshape(mu,numel(mu),1),reshape(mu_emoment,numel(mu_emoment),1),...
    'type', 'Spearman');
fprintf('Spearman correlation calculated growth rate vs eMOMENT growth rate: %.2g (p=%.2g)\n',...
    rho, pval)
[rho,pval] = corr(reshape(mu,numel(mu),1),reshape(mu_emoment,numel(mu_emoment),1));
fprintf('Pearson correlation calculated growth rate vs eMOMENT growth rate: %.2g (p=%.2g)\n\n',...
    rho, pval)

[rho, pval] = corr(reshape(mu,numel(mu),1),reshape(mu_fba,numel(mu_fba),1),...
    'type', 'Spearman');
fprintf('Spearman correlation calculated growth rate vs FBA growth rate: %.2g (p=%.2g)\n',...
    rho, pval)
[rho, pval] = corr(reshape(mu,numel(mu),1),reshape(mu_fba,numel(mu_fba),1));
fprintf('Pearson correlation calculated growth rate vs FBA growth rate: %.2g (p=%.2g)\n\n',...
    rho, pval)

[rho, pval] = corr(reshape(mu_fba,numel(mu_fba),1),reshape(mu_emoment,numel(mu_emoment),1),...
    'type', 'Spearman');
fprintf('Spearman correlation FBA growth rate vs eMOMENT growth rate: %.2g (p=%.2g)\n',...
    rho, pval)

[rho, pval] = corr(reshape(mu_fba,numel(mu_fba),1),reshape(mu_emoment,numel(mu_emoment),1));
fprintf('Pearson correlation FBA growth rate vs eMOMENT growth rate: %.2g (p=%.2g)\n\n',...
    rho, pval)

[rho, pval] = corr(reshape(mu_emoment,numel(mu_emoment),1),reshape(p_tot,numel(p_tot),1),...
    'type', 'Spearman');
fprintf('Spearman correlation eMOMENT growth rate vs measured protein content: %.2g (p=%.2g)\n',...
    rho, pval)
[rho, pval] = corr(reshape(mu_emoment,numel(mu_emoment),1),reshape(p_tot,numel(p_tot),1));
fprintf('Pearson correlation eMOMENT growth rate vs measured protein content: %.2g (p=%.2g)\n\n',...
    rho, pval)

[rho, pval] = corr(reshape(mu_fba,numel(mu_fba),1),reshape(p_tot,numel(p_tot),1),...
    'type', 'Spearman');
fprintf('Spearman correlation FBA growth rate vs measured protein content: %.2g (p=%.2g)\n',...
    rho, pval)
[rho, pval] = corr(reshape(mu_fba,numel(mu_fba),1),reshape(p_tot,numel(p_tot),1));
fprintf('Pearson correlation FBA growth rate vs measured protein content: %.2g (p=%.2g)\n\n',...
    rho, pval)