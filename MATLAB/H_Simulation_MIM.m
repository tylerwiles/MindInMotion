
% TITLE: H_Simulation_MIM.m
% DATE: 11/27/2025
% AUTHOR: Tyler M. Wiles, PhD
% EMAIL: twiles@ufl.edu

% DESCRIPTION:

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

drops = 0.02:0.02:0.10;
min_vec = [51 101 151 201];

n_files = length(my_files);
n_drops = numel(drops);
n_min = numel(min_vec);

% one H per file x min_contacts setting
h_original = NaN(n_files, n_min);
id = cell(n_files, n_min);
treadmill = cell(n_files, n_min);
speed = cell(n_files, n_min);
terrain = cell(n_files, n_min);
trial = cell(n_files, n_min);

% per file x min_contacts x %drop
drops_used = NaN(n_files, n_min, n_drops);
strides = NaN(n_files, n_min, n_drops);
strides_cut = NaN(n_files, n_min, n_drops);
strides_new_length = NaN(n_files, n_min, n_drops);
h_random = NaN(n_files, n_min, n_drops);
h_contiguous = NaN(n_files, n_min, n_drops);

parfor i = 1:length(my_files)

    % Getting filenames that will be useful for the rest of the code
    base_filename = my_files(i).name;
    full_filename = fullfile(my_files(i).folder, '\', base_filename);

    % Extracting filenames with conditions
    filename_split = strsplit(base_filename, '_');
    temp_id = strsplit(full_filename, '\');
    temp_id_condition = strsplit(base_filename, '_');

    dat_temp = pop_loadset(full_filename); % Load in the EEGLab structure
    dat = struct2table(dat_temp.event); % Just get time series of gait events

    % Cut timeseries to start and end of trial
    dat_start = find(strcmp(dat.type, 'TrialStart'));
    dat_end = find(strcmp(dat.type, 'TrialEnd'));
    dat = dat(dat_start:dat_end, :);

    % Find left and right heel contacts
    contacts_right = find(strcmp(dat.type, 'RHS'));
    contacts_latency = table2array(dat(contacts_right, 1));

    % Loop over minimum contacts thresholds
    for k = 1:n_min

        min_contacts = min_vec(k);

        % only run for this threshold if we have enough contacts
        if numel(contacts_latency) >= min_contacts

            % Cut to last min_contacts
            contacts_latency_cut = contacts_latency(max(1, numel(contacts_latency) - min_contacts + 1) : end);

            % Calculate step or stride intervals
            intervals = diff(contacts_latency_cut) / dat_temp.srate;

            % H for original series for this (file, min_contacts) pair
            h_original(i,k) = median(bayesH(intervals, 200));

            % local row buffers for this file & threshold
            drops_used_row = NaN(1, n_drops);
            strides_row = NaN(1, n_drops);
            strides_cut_row = NaN(1, n_drops);
            strides_new_length_row = NaN(1, n_drops);
            h_random_row = NaN(1, n_drops);
            h_contiguous_row = NaN(1, n_drops);
            id{i} = temp_id{6};
            treadmill{i} = temp_id_condition{1};
            speed{i} = temp_id_condition{2};
            terrain{i} = 'flat';
            trial{i} = temp_id_condition{3}(1);

            n_total = length(intervals); % total strides for this file

            % Cut strides randomly based on %
            for j = 1:numel(drops)

                drops_used_row(j) = drops(j); % Current % of strides to be droppped

                n_keep = n_total - floor(n_total * drops(j)); % How many strides to keep

                strides_row(j) = n_total; % Number of strides
                strides_cut_row(j) = n_total - n_keep; % Number of strides cut

                intervals_random = intervals(sort(randperm(n_total, n_keep))); % Cut the dataframe

                strides_new_length_row(j) = length(intervals_random);  % Length of new time series

                h_random_row(j) = median(bayesH(intervals_random, 200)); % Get H

                % Do the same for contiguous strides
                n_drop = n_total - n_keep; % Number of strides to cut

                start_idx = randi(n_total - n_drop + 1); % choose random start for contiguous block to cut

                % Create index of those to keep, then assign those that will be removed a zero
                keep_mask = true(n_total, 1);
                keep_mask(start_idx:start_idx + n_drop - 1) = false;

                intervals_contiguous = intervals(keep_mask); % Cut based on index above

                h_contiguous_row(j) = median(bayesH(intervals_contiguous, 200)); % Get H

            end

            % assign rows once per file & threshold (parfor-compatible slicing)
            strides(i,k,:) = strides_row;
            drops_used(i,k,:) = drops_used_row;
            strides_cut(i,k,:) = strides_cut_row;
            strides_new_length(i,k,:) = strides_new_length_row;
            h_random(i,k,:) = h_random_row;
            h_contiguous(i,k,:) = h_contiguous_row;

        end

    end

end

% Combine outputs from parfor into one long table
n_files = length(my_files);
n_drops = numel(drops);
n_min = numel(min_vec);

% Flatten per-(file, min, drop) matrices into column vectors
strides_col = strides(:);
drops_col = drops_used(:);
strides_cut_col = strides_cut(:);
strides_new_length_col = strides_new_length(:);
h_random_col = h_random(:);
h_contiguous_col = h_contiguous(:);

% Build corresponding min_contacts vector (same shape as strides, etc.)
min_mat = repmat(reshape(min_vec, 1, n_min, 1), [n_files, 1, n_drops]);
min_col = min_mat(:);

% Repeat original H for each drop level
h_original_rep = repmat(h_original, 1, 1, n_drops);
h_original_col = h_original_rep(:);

% Do the same for the id condition variables
id(:, 2:end) = repmat(id(:,1), 1, size(id,2)-1);
id_rep = repmat(id, 1, 1, n_drops);
id_col = id_rep(:);

treadmill(:, 2:end) = repmat(treadmill(:,1), 1, size(treadmill,2)-1);
treadmill_rep = repmat(treadmill, 1, 1, n_drops);
treadmill_col = treadmill_rep(:);

speed(:, 2:end) = repmat(speed(:,1), 1, size(speed,2)-1);
speed_rep = repmat(speed, 1, 1, n_drops);
speed_col = speed_rep(:);

terrain(:, 2:end) = repmat(terrain(:,1), 1, size(terrain,2)-1);
terrain_rep = repmat(terrain, 1, 1, n_drops);
terrain_col = terrain_rep(:);

trial(:, 2:end) = repmat(trial(:,1), 1, size(trial,2)-1);
trial_rep = repmat(trial, 1, 1, n_drops);
trial_col = trial_rep(:);

% Build table in desired column order
results = table(id_col, treadmill_col, speed_col, terrain_col, trial_col, ...
    strides_col, drops_col, strides_cut_col, strides_new_length_col, ...
    h_original_col, h_random_col, h_contiguous_col, ...
    'VariableNames', {'id', 'treadmill', 'speed', 'terrain', 'trial', ...
    'strides', 'drops', 'strides.cut', 'strides.new.length', ...
    'h.original', 'h.random', 'h.contiguous'});

% % Remove skipped trials / thresholds
% valid_idx = ~isnan(h_random_col);
% results = results(valid_idx,:);

writetable(results, fullfile(output_directory, 'H_Simulation_MIM_Results.csv'));
