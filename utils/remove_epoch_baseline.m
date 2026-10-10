function [EEG,baseline_mean]=remove_epoch_baseline(EEG,baseline_time_start,baseline_time_end)
n_trials   = EEG.trials;
n_channels = EEG.nbchan;

% Precompute baseline indices
idx_1 = interp1(EEG.times, 1:length(EEG.times), baseline_time_start, 'nearest');
idx_2 = interp1(EEG.times, 1:length(EEG.times), baseline_time_end, 'nearest');

for i = 1:n_channels
    for j = 1:n_trials
        % Compute baseline mean
        baseline_data = EEG.data(i, idx_1:idx_2, j);
        baseline_mean = mean(baseline_data, 'omitnan');
        % Subtract from all time points of this trial
        EEG.data(i,:,j) = EEG.data(i,:,j) - baseline_mean;
    end
end

