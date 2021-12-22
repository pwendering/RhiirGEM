% keep biomass at fixed value?
changeCobraSolver('ibm_cplex', 'all',0);
load(fullfile('model','iRi1574.mat'));


% find minimum uptake of palmitate at optimal biomass
minFlux = fluxVariability(model,'rxnNameList', 'r1007_e0');
model.ub(findRxnIDs(model,'r1007_e0')) = minFlux;


% find optimal growth rate
bio_idx = find(model.c);
sol = optimizeCbModel(model);
mu_opt = sol.f;

carbonSources = {
    'm0564[e0]' % D-glucose
    'm0588[e0]' % D-fructose
    'm0047[e0]' % raffinose
    'm1282[e0]' % melibiose
    'm0366[e0]' % D-xylose (Sugiura et al., PNAS Vol. 117, 2020)
    'm0028[e0]' % acetate (Pfeffer et al., Plant Physiol. Vol. 120, 1999)
    'm1331[e0]' % myristate (Sugiura et al., PNAS Vol. 117, 2020)
    'm0129[e0]' % glycine
    'm1335[e0]' % gycerol (Bago et al., Plant Physiol. Vol. 131, 2003)
    'm0239[c0]' % trehalose
};

oldUptake = {
    'r1533_e0' % D-Glucose
    'r1534_e0' % D-Fructose
    'r1006_e0' % myo-Inositol
    'r1002_e0' % Glycine
    'r1633_e0' % Myristate
    };
cSourceNames = regexprep(model.metNames(findMetIDs(model, carbonSources)),'_','-');
% remove all carbon sources
model = removeRxns(model,oldUptake,'metFlag',false);

% add sink reaction for cytosolic ATP 
model = addSinkReactions(model,'m0464[c0]',0,1000);
model.c(:) = 0;
model.c(end) = 1;

mu_perc = [50 0];

atp_production = zeros(numel(carbonSources),numel(mu_perc));
carbon_source_upt = zeros(numel(carbonSources),numel(mu_perc));

for i=1:numel(carbonSources)
    for j=1:numel(mu_perc)
        tmpModel = addExchangeRxn(model,carbonSources(i));
        tmpModel.lb(bio_idx) = mu_perc(j)*mu_opt/100;
        s = optimizeCbModel(tmpModel);
        if s.stat == 1
            atp_production(i,j) = s.f;
            carbon_source_upt(i,j) = abs(s.x(end));
        end
    end
end

colors = [[142,153,204];[236,176,80]]/255;

labels = categorical(cSourceNames);
labels = reordercats(labels,cSourceNames);

subplot(1,2,1)
b=bar(labels,atp_production');
ylabel('ATP production [mmol/gDW/h]','FontSize',14)
set(gca,'LineWidth',1.3,'FontSize',14)
l=legend(strsplit(strtrim(sprintf('%d%% ',100-mu_perc))),'LineWidth',1,...
    'FontSize',10,'NumColumns',2,'Location','northwest','Box','off');
title(l,'Allowed growth reduction')
b(1).FaceColor = colors(1,:);
b(2).FaceColor = colors(2,:);
text(-0.2,1.01,'A','units','normalized','fontweight','bold','fontsize',14)

subplot(1,2,2)
b=bar(labels,atp_production' ./ carbon_source_upt');
ylabel('ATP yield','FontSize',14)
set(gca,'LineWidth',1.3, 'FontSize',14)
l=legend(strsplit(strtrim(sprintf('%d%% ',100-mu_perc))),'LineWidth',1,...
    'FontSize',10,'NumColumns',2,'Location','northwest','Box','off');
title(l,'Allowed growth reduction')
b(1).FaceColor = colors(1,:);
b(2).FaceColor = colors(2,:);
text(-0.2,1.01,'B','units','normalized','fontweight','bold','fontsize',14)

print(fullfile('results','figures','atp-production.png'),'-painters','-dpng')