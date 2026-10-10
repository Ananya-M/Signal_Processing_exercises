function [all_idx, coi_idx, coi_labels] = get_chan_numbers_labels(chanLabels, desiredLabels)
% chanLabels     : cell array of all channel names from EEG
% desiredLabels  : cell array of channels you want (e.g. {'C3','Cz',...})
%
% Outputs:
% all_idx    : indices of all channels (1:N)
% coi_idx    : indices of desired channels within chanLabels
% coi_labels : the matching channel names found

    all_idx = 1:numel(chanLabels);

    coi_idx = [];
    coi_labels = {};

    for i = 1:numel(desiredLabels)
        idx = find(strcmpi(chanLabels, desiredLabels{i}));
        if ~isempty(idx)
            coi_idx(end+1) = idx; %#ok<AGROW>
            coi_labels{end+1} = chanLabels{idx}; %#ok<AGROW>
        end
    end
end