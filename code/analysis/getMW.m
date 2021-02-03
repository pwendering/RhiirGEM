function MW = getMW(formulae)
% Function to calculate the moelcular weight for the metabolite in a
% metabolic model, can only parse simple formulae without bracket
% Input:
%     cellstr formulae:     molecular formulae
% Output:
%     double MW:            array that contains the molecular weights

% dictionary of atom weights
atomWeights=containers.Map;
elements={'C','H','O','N','S','P','K','Ca','Cl',...
    'Co','Cu','Fe','I','K','Mg','Mn','Mo','Na','Zn'};
weights=[12.01 1.01 16 14.01 32.06 30.97 39.10 40.08 35.45 58.93 63.55 55.84 ...
    126.90 39.10 24.30 54.94 95.94 22.99 65.39];

% initialize molecular weights array
MW=zeros(size(formulae));

for i=1:numel(elements)
    atomWeights(elements{i})=weights(i);
end

for i=1:numel(formulae)
    
    % split elements with associated numbers
    formulaSplit=regexp(formulae{i},'[A-Z][a-z]?(\d+(\.\d+)?)?','match');
    
    
    % extract number of atoms per element
    atomNums=regexp(formulaSplit,'\d+(\.\d+)?','match');
    atomNums=cellfun(@str2double,atomNums,'un',0);
    atomNums(cellfun(@isempty,atomNums))={1};
    atomNums=cell2mat(atomNums);
    
    % find molecular weights for each element
    atomSymbols=cellfun(@(x)regexp(x,'[A-Z][a-z]?','match'),formulaSplit);
    
    MW(i)=0;
    for j=1:numel(atomSymbols)
        MW(i)=MW(i)+atomWeights(atomSymbols{j})*atomNums(j);
    end
    
end


end

