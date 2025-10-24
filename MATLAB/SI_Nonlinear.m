
% TITLE: SI_Nonlinear.m
% DATE: 10/22/2025
% AUTHOR: Tyler M. Wiles, PhD
% EMAIL: twiles@ufl.edu

% DESCRIPTION:
% Calculate Hurst Exponent and Sample Entropy from stride intervals.
% R:\Ferris-Lab\share\MindInMotion\IMU_Validation\Code
% R:\Ferris-Lab\share\MindInMotion\Scripts

% NOTES/QUESTIONS:
% What does *_FU, _Airram, etc mean?
% Where do I obtain TM belt speed?
% Obtain metrics for both legs and then compare between legs in R.
% Everyone is cut to the same number of strides for both conditions. Do one
% for each condition?
% Correlate ratings of perceived stability with hurst and entorpy?

clear; close all; clc;

addpath(genpath('FUNCTIONS'));
addpath(genpath('DATA'));
addpath(genpath('MATLAB'));
my_directory = 'R:\Ferris-Lab\share\MindInMotion\Data';
output_directory = 'C:\Users\tyler_3r9w1ip\Desktop\Github_Repositories\MindinMotion\DATA';

%% Get all subject folders then files per participant
my_folders = dir(fullfile(my_directory, 'H*'));
my_folders = [my_folders; dir(fullfile(my_directory, 'NH*'))];
my_files = [];
files = cell(length(my_folders),1);
parfor i = 1:length(my_folders)
    disp(['Pulling: ', my_folders(i).name])
    this_path = fullfile(my_folders(i).folder, my_folders(i).name, ...
        'Loadsol', 'Processed_Trials', 'Stride_by_Stride_Gait_Structs', '*R.mat');
    files{i} = dir(this_path);
end
my_files = vertcat(files{:});

%% Find minimum number of strides

parfor i = 1:length(my_files)

    % Getting filenames that will be useful for the rest of the code
    base_filename = my_files(i).name;
    full_filename = fullfile(my_files(i).folder, '\', base_filename);

    disp(['Analyzing: ', base_filename]); % Display file being used

    % Extracting filenames with conditions
    filename_split = strsplit(base_filename, '_');
    temp_id = strsplit(full_filename, '\'); id{i,:} = temp_id{6};
    temp_id_condition = strsplit(base_filename, '_'); treadmill{i,:} = temp_id_condition{1};
    trial{i,:} = temp_id_condition{3};
    if strcmp(temp_id_condition{1}, 'SP')
        speed{i,:} = temp_id_condition{2};
        terrain{i,:} = 'NaN';
    else
        speed{i,:} = 'NaN';
        terrain{i,:} = temp_id_condition{2};
    end

    dat = load(full_filename);

    n_stride_intervals(i,:) = length([dat.RarrayOfGaitStruct.GaitCycleDur]); % Extract pre-calculated stride intervals

end

min_strides = table(id, treadmill, speed, terrain, trial, n_stride_intervals, ...
    'VariableNames', {'id', 'treadmill', 'speed', 'terrain', 'trial', 'n_stride_intervals'});
writetable(min_strides, fullfile(output_directory, 'min_strides.csv'));

%% Hurst and Entropy

min_strides = 64;

parfor i = 1:length(my_files)

    % Getting filenames that will be useful for the rest of the code
    base_filename = my_files(i).name;
    full_filename = fullfile(my_files(i).folder, '\', base_filename);

    disp(['Analyzing: ', base_filename]); % Display file being used

    % Extracting filenames with conditions
    filename_split = strsplit(base_filename, '_');
    temp_id = strsplit(full_filename, '\'); id{i,:} = temp_id{6};
    temp_id_condition = strsplit(base_filename, '_'); treadmill{i,:} = temp_id_condition{1};
    trial{i,:} = temp_id_condition{3};
    if strcmp(temp_id_condition{1}, 'SP')
        speed{i,:} = temp_id_condition{2};
        terrain{i,:} = 'NaN';
    else
        speed{i,:} = 'NaN';
        terrain{i,:} = temp_id_condition{2};
    end

    dat = load(full_filename); % Load data

    stride_intervals = [dat.RarrayOfGaitStruct.GaitCycleDur]'; % Extract pre-calculated stride intervals

    % If minimum number of stride_intervals are not met then add NaN to
    % nonlinear results.
    if length(stride_intervals) < min_strides
        hursts(i,:) = NaN;
        entropies(i,:) = NaN;
    else

        stride_intervals = stride_intervals(1:min_strides); % Cut to first # of minimum strides

        % Hurst Exponent
        hurst = median(bayesH(stride_intervals, 200));
        hursts(i,:) = hurst;

        % Entropy
        entropy = Samp_En(stride_intervals, 2 , 0.25, std(stride_intervals));
        entropies(i,:) = entropy;

        % Plot/save stride intervals to check
        f = figure('Visible', 'off');
        plot(stride_intervals);
        % ylim([1, 2]);
        save_name = sprintf('%s_%s_%s_%s.png', temp_id{6}, temp_id_condition{1}, temp_id_condition{2}, temp_id_condition{3});
        save_path = fullfile(output_directory, '\FIGURES\', save_name);
        saveas(f, save_path);
        close(f);

    end

end

% Export
mim_results = table(id, treadmill, speed, terrain, trial, hursts, entropies, ...
    'VariableNames', {'id', 'treadmill', 'speed', 'terrain', 'trial', 'hurst', 'entropy'});
writetable(mim_results, fullfile(output_directory, 'MIM_Nonlinear_Results.csv'));
