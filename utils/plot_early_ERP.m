function plot_early_ERP(EEG, dataField)
% plot_early_ERP - Plots averaged ERP for selected channels from EEG struct
%
% Usage:
%   plot_early_ERP(EEG, 'data') 
%   plot_early_ERP(EEG, 'data_without_badtrials') 
%   plot_early_ERP(EEG, 'zscoreddata') 
%
% Inputs:
%   EEG       - EEG struct containing the data
%   dataField - (string) which field of EEG to use ('data', 
%                'data_without_badtrials', or 'zscoreddata')

    % Validate input field
    if ~isfield(EEG, dataField)
        error('Field "%s" not found in EEG struct.', dataField);
    end

    % Extract the EEG data
    EEG_data = EEG.(dataField);

    % Get channel labels
    chanLabels = {EEG.chanlocs.labels};

    % Time vector
    EEG_time = EEG.times;

    % Define channels of interest
    D = {'C3','C2','Cz','CP5','CP2','CP1','CP6','C4','Fz','F4'};

    % Get indices of channels of interest
    [~, coi, coi_chan] = get_chan_numbers_labels(chanLabels, D);

    % Preallocate ERP matrix
    ERP_1 = nan(length(coi), size(EEG_data,2));

    % Compute average ERP per channel
    for i = 1:length(coi)
        ERP_1(i,:) = squeeze(nanmean(EEG_data(coi(i), :, :), 3));
    end

    % Plot results
    figure;
    for i = 1:length(coi_chan)
        plot(EEG_time, ERP_1(i,:), 'LineWidth', 2, ...
            'DisplayName', ['Mean: ', coi_chan{i}]); 
        hold on;
    end

    % Add vertical line at time zero
    xline(0, '--b', 'LineWidth', 2, 'DisplayName','Stimulus');

    % Set labels and title
    xlabel('Time (ms)');
    ylabel('Amplitude (\muV)');
    title(sprintf('ERP (%s)', dataField));

    % Adjust axes
    xlim([-300 300]);
    legend('show', 'Location', 'bestoutside');
    grid on;
end