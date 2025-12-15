
% TITLE: H_Simulation_MIM_Hypothesis.m
% DATE: 11/27/2025
% AUTHOR: Tyler M. Wiles, PhD
% EMAIL: twiles@ufl.edu

% DESCRIPTION:
% Estimate the Hurst exponent and Sample Entropy from stride intervals
% using the mind in motion dataset. Outliers are removed based on a
% standard deviation. These results are for a mini example of using data
% that sometimes needs to be cut due to missing heel contacts, and then
% comparing between high and low functioning adults.

clear; close all; clc;

addpath(genpath('FUNCTIONS'));
addpath(genpath('DATA'));
addpath(genpath('MATLAB'));
addpath(genpath('R:\Ferris-Lab\share\MindInMotion\eeglab'));
eeglab % Load EEGLAB
my_directory = 'R:\Ferris-Lab\share\MindInMotion\Data';
output_directory = 'C:\Users\twiles\Desktop\Github-Repositories\MindInMotion\DATA\H_Simulation';

my_folders = dir(fullfile(my_directory, 'H*'));
my_folders = [my_folders; dir(fullfile(my_directory, 'NH*'))];
my_folders = my_folders(~contains({my_folders.name}, '_')); % Remove *_FU, _Airram, etc
files = cell(length(my_folders),1);
parfor i = 1:length(my_folders)
    disp(['Pulling: ', my_folders(i).name])
    this_path = fullfile(my_folders(i).folder, my_folders(i).name, ...
        'EEG', 'Trials');
    files_all = dir(fullfile(this_path, '*.set'));
    files_s = files_all(startsWith({files_all.name}, 'SP'));
    files{i} = files_s;
end
my_files = vertcat(files{:});
my_files(strcmp({my_files.name}, 'TM_low_2.set') & strcmp({my_files.folder}, 'R:\Ferris-Lab\share\MindInMotion\Data\NH3027\EEG\Trials')) = []; % No data for this trial?

min_contacts = 51;

parfor i = 1:length(my_files)

    % Getting filenames that will be useful for the rest of the code
    base_filename = my_files(i).name;
    full_filename = fullfile(my_files(i).folder, '\', base_filename);

    disp(['Analyzing: ', base_filename]); % Display file being used

    % Extracting filenames with conditions
    filename_split = strsplit(base_filename, '_');
    temp_id = strsplit(full_filename, '\'); id{i,:} = temp_id{6};
    temp_id_condition = strsplit(base_filename, '_'); treadmill{i,:} = temp_id_condition{1};
    trial{i,:} = temp_id_condition{3}(1);
    if strcmp(temp_id_condition{1}, 'SP')
        speed{i,:} = temp_id_condition{2};
        terrain{i,:} = 'flat';
    else
        speed{i,:} = 'NaN';
        terrain{i,:} = temp_id_condition{2};
    end

    dat_temp = pop_loadset(full_filename); % Load in the EEGLab structure
    dat = struct2table(dat_temp.event); % Just get time series of gait events

    % Cut timeseries to start and end of trial
    dat_start = find(strcmp(dat.type, 'TrialStart'));
    dat_end = find(strcmp(dat.type, 'TrialEnd'));
    dat = dat(dat_start:dat_end, :);

    % Find left and right heel contacts
    contacts_right = find(strcmp(dat.type, 'RHS'));
    contacts_latency = table2array(dat(contacts_right, 1));

    % If minimum number of contacts are not met then add NaN.
    if length(contacts_right) < min_contacts
        hursts(i,:) = NaN;
        entropies(i,:) = NaN;
        intervals_mean(i,:) = NaN;
        intervals_sd(i,:) = NaN;
        n_intervals = NaN;
    else
        % contacts_latency = contacts_latency(1:min_contacts); % Cut to first # of minimum strides
        contacts_latency = contacts_latency(max(1, numel(contacts_latency) - min_contacts + 1):end); % Cut to last # of minimum strides

        % Calculate step or stride intervals
        intervals = diff(contacts_latency)/dat_temp.srate;

        % Remove intervals that are greater/lower than three standard
        % deviations from the mean. Record number of dropped intervals
        intervals_length_original_temp = numel(intervals);
        intervals_length_original(i,:) = numel(intervals);
        intervals_mean_temp = mean(intervals);
        intervals_mean(i,:) = intervals_mean_temp;
        intervals_sd_temp = std(intervals);
        intervals_sd(i,:) = intervals_sd_temp;
        intervals_thresh_upper = intervals_mean_temp + 3 * intervals_sd_temp;
        intervals_thresh_lower = intervals_mean_temp - 3 * intervals_sd_temp;
        intervals = intervals(intervals >= intervals_thresh_lower & intervals <= intervals_thresh_upper);
        intervals_dropped(i,:) = intervals_length_original_temp - numel(intervals);
        pct_timeseries(i,:) = ((intervals_length_original_temp - numel(intervals)) / intervals_length_original_temp) * 100;

        % Hurst Exponent
        hurst = median(bayesH(intervals, 200));
        hursts(i,:) = hurst;

        % Entropy
        entropy = Samp_En(intervals, 2 , 0.25, std(intervals));
        entropies(i,:) = entropy;

        % Plot/save stride intervals to check
        f = figure('Visible', 'off');
        plot(intervals);
        % ylim([1, 2]);
        save_name = sprintf('%s_%s_%s_%s.png', temp_id{6}, temp_id_condition{1}, temp_id_condition{2}, temp_id_condition{3});
        save_path = fullfile(output_directory, '\FIGURES\', save_name);
        saveas(f, save_path);
        close(f);

    end

end

% Export
mim_results = table(id, treadmill, speed, terrain, trial, ...
    hursts, entropies, ...
    intervals_length_original, intervals_dropped, pct_timeseries, intervals_mean, intervals_sd, ...
    'VariableNames', {'id', 'treadmill', 'speed', 'terrain', 'trial', ...
    'hurst', 'entropy', ...
    'n.intervals.original', 'n.intervals.dropped', 'pct.timeseries.dropped', 'intervals.mean', 'intervals.sd'});
writetable(mim_results, fullfile(output_directory, 'H_Simulation_MIM_Hypothesis.csv'));
