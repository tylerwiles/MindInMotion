
% TITLE: H_Simulation_MIM_Hypothesis.m
% DATE: 11/27/2025
% AUTHOR: Tyler M. Wiles, PhD
% EMAIL: twiles@ufl.edu

% DESCRIPTION:
% Estimate the Hurst exponent and Sample Entropy from stride intervals
% using the mind in motion dataset. Individual files of each trial's stride
% intervals have already been run using H_Simulation_MIM.m. Outliers are
% removed based on a standard deviation. These results are for a mini
% example of using data that sometimes needs to be cut due to missing heel
% contacts, and then comparing between high and low functioning adults.
% Hurst exponent uses median of 200 samples and Entropy uses 35% of the
% standard deviation.

clear; close all; clc;

addpath(genpath('FUNCTIONS'));
addpath(genpath('DATA'));
addpath(genpath('MATLAB'));
my_directory = 'R:\Ferris-Lab\twiles\H_Simulation\MIM';
output_directory = 'C:\Users\twiles\Desktop\Github-Repositories\MindInMotion\DATA\H_Simulation';
my_files = dir(fullfile(my_directory, '*.csv'))

min_stride_intervals = 50;

parfor i = 1:length(my_files)

    % Getting filenames that will be useful for the rest of the code
    base_filename = my_files(i).name;
    full_filename = fullfile(my_files(i).folder, '\', base_filename);

    disp(['Analyzing: ', base_filename]); % Display file being used

    % Extracting filenames with conditions
    filename_split = strsplit(base_filename, '_');
    temp_id = strsplit(full_filename, '\');
    temp_id_condition = strsplit(base_filename, '_');
    id{i,:} = temp_id_condition{1};
    treadmill{i,:} = temp_id_condition{2};
    trial{i,:} = temp_id_condition{5}(1);
    if strcmp(temp_id_condition{2}, 'SP')
        speed{i,:} = temp_id_condition{3};
        terrain{i,:} = 'flat';
    else
        speed{i,:} = 'NaN';
        terrain{i,:} = temp_id_condition{5};
    end

    dat = readmatrix(full_filename); % Read in data

    % If minimum number of contacts are not met then add NaN.
    if length(dat) < min_stride_intervals
        hursts(i,:) = NaN;
        entropies(i,:) = NaN;
        intervals_mean(i,:) = NaN;
        intervals_sd(i,:) = NaN;
        n_intervals = NaN;
    else
        
        dat = dat(max(1, numel(dat) - min_stride_intervals + 1):end); % Cut to last # of minimum stride intervals

        % Remove intervals that are greater/lower than three standard
        % deviations from the mean. Record number of dropped intervals
        intervals_length_original_temp = numel(dat);
        intervals_length_original(i,:) = numel(dat);

        intervals_mean_temp = mean(dat);
        intervals_sd_temp = std(dat);

        intervals_thresh_upper = intervals_mean_temp + 3 * intervals_sd_temp;
        intervals_thresh_lower = intervals_mean_temp - 3 * intervals_sd_temp;
        dat = dat(dat >= intervals_thresh_lower & dat <= intervals_thresh_upper); % Removed here

        intervals_mean(i,:) = mean(dat); % Get mean and SD after removing outliers
        intervals_sd(i,:) = std(dat);

        intervals_dropped(i,:) = intervals_length_original_temp - numel(dat);
        pct_timeseries(i,:) = ((intervals_length_original_temp - numel(dat)) / intervals_length_original_temp) * 100;

        % Hurst Exponent
        hurst = median(bayesH(dat, 200));
        hursts(i,:) = hurst;

        % Entropy
        entropy = Samp_En(dat, 2 , 0.35, std(dat));
        entropies(i,:) = entropy;

        % Plot/save stride intervals to check
        f = figure('Visible', 'off');
        plot(dat);
        % ylim([1, 2]);
        save_name = sprintf('%s_%s_%s_%s_%s.png', temp_id_condition{1}, temp_id_condition{2}, temp_id_condition{3}, temp_id_condition{4}, temp_id_condition{5}(1));
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
