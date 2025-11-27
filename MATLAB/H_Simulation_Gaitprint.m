
% TITLE: H_Simulation_Gaitprint.m
% DATE: November 26, 2025
% AUTHOR: Tyler M. Wiles, PhD
% EMAIL: twiles@ufl.edu

% DESCRIPTION:
% Get the Hurst exponent from right stride intervals of all heatlhy adults
% that completed four minute walking trials with bad heel contacts.
% Randomly and contiguously drops 2-10% of the full trial's heel contacts
% before getting H.

clear; close all; clc;

addpath(genpath('FUNCTIONS'));
addpath(genpath('DATA'));
addpath(genpath('MATLAB'));
my_directory = 'E:\Research\Gaitprint\CLEAN DATA';
output_directory = 'C:\Users\tyler_3r9w1ip\Desktop\Github_Repositories\MindInMotion\DATA\H_Simulation';
my_files = dir(fullfile(my_directory, 'S*'));
my_files = my_files(~contains({my_files.name}, 'G05'));

ref = readtable("C:\Users\tyler_3r9w1ip\Desktop\Github_Repositories\MindInMotion\DATA\H_Simulation\Gaitprint_trial_characteristics.csv");
ref(:,19) = strcat(ref{:,1}, {'_'}, ref{:,2}, {'_'}, ref{:,3}, {'_'}, ref{:,4}, {'_'}, ref{:,5});

drops = 0.02:0.02:0.10;

n_files = length(my_files);
n_drops = numel(drops);

id_condition = cell(n_files,1);
h_original = NaN(n_files,1);

drops_used = NaN(n_files,n_drops);
strides = NaN(n_files,n_drops);
strides_cut = NaN(n_files,n_drops);
strides_new_length = NaN(n_files,n_drops);
h_random = NaN(n_files,n_drops);
h_contiguous = NaN(n_files,n_drops);

parfor i = 1:length(my_files)

    % Getting file names that will be useful for the rest of the code
    base_filename = my_files(i).name;
    full_filename = fullfile(my_directory, base_filename);

    % Extracting filenames with conditions
    filename_split = strsplit(base_filename,'.');

    % Saved figure names based on id conditions
    png_filename = filename_split{1};

    disp('Analyzing'); % Display file being used
    disp(base_filename);

    % Skip trials with bad heel contacts or not 4 minute trials
    ref_index = find(strcmp(ref{:,19}, png_filename));
    if ref{ref_index, "contacts"} == 0 || ref{ref_index, "full_trial"} == 0
        continue
    else

        id_condition{i,1} = png_filename;

        dat = readmatrix(full_filename); % Read in data

        [right_contact_locs, ~, ~, ~] = Gaitprint_Find_Gait_Events(dat, base_filename); % Get right heel contacts

        stride_intervals = diff(right_contact_locs) / 200; % Get stride interval divided by sampling rate

        h_original(i,:) = median(bayesH(stride_intervals, 200)); % Get Hurst exponent for original timeseries

        % local row buffers for this file
        drops_used_row = NaN(1,n_drops);
        strides_row = NaN(1,n_drops);
        strides_cut_row = NaN(1,n_drops);
        strides_new_length_row = NaN(1,n_drops);
        h_random_row = NaN(1,n_drops);
        h_contiguous_row = NaN(1,n_drops);

        n_total = length(stride_intervals); % total strides for this file

        % Cut strides randomly based on %
        for j = 1:numel(drops)

            drops_used_row(j) = drops(j); % Current % of strides to be droppped

            n_keep = n_total - floor(n_total * drops(j)); % How many strides to keep

            strides_row(j) = n_total; % Number of strides

            strides_cut_row(j) = n_total - n_keep; % Number of strides cut

            stride_intervals_random = stride_intervals(sort(randperm(n_total, n_keep))); % Cut the dataframe

            strides_new_length_row(j) = length(stride_intervals_random); % Length of new time series

            h_random_row(j) = median(bayesH(stride_intervals_random, 200)); % Get H

            % Do the same for contiguous strides
            n_drop = n_total - n_keep; % Number of strides to cut

            start_idx = randi(n_total - n_drop + 1); % choose random start for contiguous block to cut

            % Create index of those to keep, then assign those that will be removed
            % a zero
            keep_mask = true(n_total, 1);
            keep_mask(start_idx:start_idx + n_drop - 1) = false;

            stride_intervals_contiguous = stride_intervals(keep_mask); % Cut based on index above

            h_contiguous_row(j) = median(bayesH(stride_intervals_contiguous, 200)); % Get H

        end

        % assign rows once per file
        strides(i,:) = strides_row;
        drops_used(i,:) = drops_used_row;
        strides_cut(i,:) = strides_cut_row;
        strides_new_length(i,:) = strides_new_length_row;
        h_random(i,:) = h_random_row;
        h_contiguous(i,:) = h_contiguous_row;

    end

end

% Combine outputs from parfor into one long table
n_files = length(my_files);
n_drops = numel(drops);

% Flatten per-(file,drop) matrices into column vectors
strides_col = strides(:);
drops_col = drops_used(:);
strides_cut_col = strides_cut(:);
strides_new_length_col = strides_new_length(:);
h_random_col = h_random(:);
h_contiguous_col = h_contiguous(:);

% Repeat trial info for each drop level
id_condition_rep = repmat(id_condition, 1, n_drops);
id_condition_col = id_condition_rep(:);

% Repeat original Hurst exponent for each drop level
h_original_rep = repmat(h_original, 1, n_drops);
h_original_col = h_original_rep(:);

% Build table in desired column order
results = table(id_condition_col, strides_col, ...
    drops_col, strides_cut_col, ...
    strides_new_length_col, h_original_col, ...
    h_random_col, h_contiguous_col, ...
    'VariableNames', {'id.condition', 'strides', ...
    'drops', 'strides.cut', ...
    'strides.new.length', 'h.original', ...
    'h.random', 'h.contiguous'});

% Remove skipped trials
valid_idx = ~isnan(h_random_col);
results = results(valid_idx,:);

writetable(results, fullfile(output_directory, 'H_Simulation_Gaitprint_Results.csv'));