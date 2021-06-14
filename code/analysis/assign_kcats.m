function kcats = assign_kcats(model, kcatFile, taxTerm, outFileName)
%% kcats = assign_kcats(model, kcatFile, taxTerm, outFileName)
% Either reads kcats from file or assigns kcats based on matching EC numbers,
% substrate matches and phylogenetic distance.
% Input:
%   struct model:           COBRA metabolic model structure containing only
%                           irreversible reactions
%   char kcatFile:          path to file containing turnover numbers
%                           columns: E.C.-Number;substrate names;taxonomic
%                           lineage;turnover value;'*'
%                           (format as used in GECKO (Sanchez et al. 2017))
%                           the file has no header line!
%   char taxTerm:           taxonomic term that should be queried for in
%                           the organisms' lineages
%   char outFileName:       (optional) path to the file where kcats values should be
%                           read from or saved after assignment (if empty,
%                           the function will only return the kcats array)
% Output:
%   double kcats:           turnover values for every reaction, if
%                           unassigned: kcat==0

req_fields = {'S', 'rxns', 'metNames', 'rxnECNumbers'};

if nargin < 2
    error('USAGE: assign_kcats(model, filename, kcat_file)')
elseif any(~isfield(model, req_fields))
    error(['The model input must have the fields: ', strjoin(req_fields, ', ')])
elseif exist('outFileName', 'file')
    fid = fopen(outFileName, 'r');
    kcats = fscanf(fid, '%f');
    fclose(fid); clear fid
    
    if numel(model.rxns) ~= numel(kcats)
        disp('A file with reaction kcats exists but dimensions do not match.')
    else
        return
    end
end

if ~exist(kcatFile, 'file')
    error('Cannot read reaction-specific kcats and and mapping file does not exist.')
end

if nargin < 4 || ~ischar(taxTerm)
    taxTerm = '';
end
    
% find the substrates from each reaction
substrates = cell(numel(model.rxns), 1);
for i=1:numel(model.rxns)
    tmpSubstrateIdx = model.S(:, i) < 0;
    tmpSubstrates=model.metNames(tmpSubstrateIdx)';
    % convert substrate names to a universal format
    tmpSubstrates = regexprep(tmpSubstrates, '^\(', '');
    tmpSubstrates = regexprep(tmpSubstrates, '[\ \-\+\(\)\:\;\,\[\]]+', '_');
    tmpSubstrates = regexprep(tmpSubstrates, '__', '_');
    tmpSubstrates = regexprep(tmpSubstrates, '_$', '');
    substrates{i} = lower(regexprep(tmpSubstrates, '^_', ''));
end

kcats = findKcatFromEC(model,substrates,kcatFile,taxTerm);

if ~isempty(outFileName)
    fid = fopen(outFileName, 'w');
    fprintf(fid, '%.10f\n', kcats);
    fclose(fid); clear fid
end


end