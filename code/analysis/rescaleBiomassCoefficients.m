function model = rescaleBiomassCoefficients(model,metID,newMF)
%% rescaleBiomassCoefficients(model,metID,newMF)
% Function to adjust the biomass coefficients if the mass fraction of one
% component is changed. The remaining mass fraction (1-newMassFraction) is
% distributed among all other biomass components according their proportions 
% on the remaining mass fraction (i.e. the relative differences between 
% coefficients are conserved).
% The function requires that all biomass components have a molecular
% formula assigned.
% Input:
%       struct model:               COBRA formatted metabolic model
%       char metID:                 compound ID associated to the new mass
%                                   fraction
%       double newMF:               new mass fraction of the compound above
%                                   must be a value between 0 and 1 [g gDW-1]
% Output:
%       struct model:               model with rescaled biomass coefficients
%   

% find biomass reaction and associated metabolites:
bioRxnIdx = model.c==1;
bioCompIdx = model.S(:,bioRxnIdx)<0;

% exclude the GAM from the calculation as it is mass-balanced by itself
gamCpdIdx = ismember(model.metNames,{'atp','h2o'});
bioCompIdx(gamCpdIdx) = 0;
biomassComponents = model.mets(bioCompIdx);

% initialize array of new coefficients
newCoeff = zeros(numel(biomassComponents),1);

% check if given compound ID is part of the biomass reaction
if ~ismember(metID,biomassComponents)
    error('The given metabolite is not contained in the biomass reaction')
end

% check if the mass fraction is between 0 and 1
if ~newMF>0&&~newMF<1
    error('The new mass fraction must be in the interval (0,1)')
end

% get molecular weights for each biomass component [g mmol-1]
MW = getMW(model.metFormulas(bioCompIdx));
MW = MW/1000;

% get old mass fractions for every biomass component [g gDW-1]
coeff = -model.S(bioCompIdx,bioRxnIdx);
oldMF = coeff.*MW;
bioMassSum = sum(oldMF);

% calculate the new coefficient for the given component
cpdBioIdx = ismember(biomassComponents,metID);
newCoeff(cpdBioIdx) = newMF/MW(cpdBioIdx);

% calculate the new coefficients for the remaining components:
% remaining mass fraction
massRemain = bioMassSum-newMF;

% scaled fractions to total mass
fractions = oldMF(~cpdBioIdx)/sum(oldMF(~cpdBioIdx));

% increase or reduction of mass fractions will be proportional to the
% fractions of each component on total mass
newRemainMF = fractions.*massRemain;

% calculate the stoichiometrix coefficients by dividing by the respective 
% molecular masses
newCoeff(~cpdBioIdx) = newRemainMF./MW(~cpdBioIdx);

% assign new coefficients
model.S(bioCompIdx,bioRxnIdx) = -newCoeff;

end