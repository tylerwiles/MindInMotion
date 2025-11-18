
% TITLE: SI_Nonlinear.m
% DATE: 10/22/2025
% AUTHOR: Tyler M. Wiles, PhD
% EMAIL: twiles@ufl.edu

% DESCRIPTION:
% Calculate Hurst Exponent and Sample Entropy from stride intervals.
% R:\Ferris-Lab\share\MindInMotion\IMU_Validation\Code
% R:\Ferris-Lab\share\MindInMotion\Scripts

% NOTES/QUESTIONS:
% Check if strides were removed before the EEGLab occurs
% Follow-ups are not included.
% Which python script finds the heel contacts?
% Ask Jake or Ryan to get the gait events code. Ryan Downing. - asked Jake
% Find minimum number of strides for each chunk of time when Jake gives me
% the code and I can do them myself
% Some people explored in the uneven terrain.
% Obtain metrics for both legs and then compare between legs in R?
% Switch to step intervals?
% Everyone is cut to the same number of strides for both conditions. Do one
% for each condition? (SP/TM)
% Correlate ratings of perceived stability with hurst and entorpy?
% Opals are more reliable than loadsol probably.
% Chris Hass might have foot dominance data.

% H1022_SP_0p25_1 has weird SI Many participants do?

clear; close all; clc;

addpath(genpath('FUNCTIONS'));
addpath(genpath('DATA'));
addpath(genpath('MATLAB'));
addpath(genpath('R:\Ferris-Lab\share\MindInMotion\eeglab'));
eeglab % Load EEGLAB
my_directory = 'R:\Ferris-Lab\share\MindInMotion\Data';
% output_directory = 'C:\Users\tyler_3r9w1ip\Desktop\Github_Repositories\MindinMotion\DATA';
output_directory = 'C:\Users\twiles\Desktop\Github-Repositories\MindInMotion\DATA';

%% Find minimum number of contacts for the entire trial

my_folders = dir(fullfile(my_directory, 'H*'));
my_folders = [my_folders; dir(fullfile(my_directory, 'NH*'))];
my_folders = my_folders(~contains({my_folders.name}, '_')); % Remove *_FU, _Airram, etc
my_files = [];
files = cell(length(my_folders),1);
parfor i = 1:length(my_folders)
    disp(['Pulling: ', my_folders(i).name])
    this_path = fullfile(my_folders(i).folder, my_folders(i).name, ...
        'EEG', 'Trials');
    files_all = dir(fullfile(this_path, '*.set'));
    files_s = files_all(startsWith({files_all.name}, 'SP'));
    files_t = files_all(startsWith({files_all.name}, 'TM'));
    files{i} = [files_s; files_t];
end
my_files = vertcat(files{:});

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
    dat_start = find(strcmp(dat{:,7}, 'TrialStart'));
    dat_end = find(strcmp(dat{:,7}, 'TrialEnd'));
    dat = dat(dat_start:dat_end, :);

    % Find left and right heel contacts
    contacts_left = find(strcmp(dat{:,7}, 'LHS'));
    contacts_right = find(strcmp(dat{:,7}, 'RHS'));
    contacts = sort([contacts_left; contacts_right], 'ascend');

    % Get total number
    n_steps(i,:) = length(contacts);
    n_contacts_left(i,:) = length(contacts_left);
    n_contacts_right(i,:) = length(contacts_right);

end

n_contacts_diff = n_contacts_left - n_contacts_right;

min_contacts = table(id, treadmill, speed, terrain, trial, n_steps, n_contacts_left, n_contacts_right, n_contacts_diff, ...
    'VariableNames', {'id', 'treadmill', 'speed', 'terrain', 'trial', 'n.contacts', 'n.contacts.left', 'n.contacts.right', 'n.contacts.diff'});
writetable(min_contacts, fullfile(output_directory, 'min_contacts.csv'));

%% Hurst and Entropy

clear; close all; clc;

addpath(genpath('FUNCTIONS'));
addpath(genpath('DATA'));
addpath(genpath('MATLAB'));
addpath(genpath('R:\Ferris-Lab\share\MindInMotion\eeglab'));
eeglab % Load EEGLAB
my_directory = 'R:\Ferris-Lab\share\MindInMotion\Data';
% output_directory = 'C:\Users\tyler_3r9w1ip\Desktop\Github_Repositories\MindinMotion\DATA';
output_directory = 'C:\Users\twiles\Desktop\Github-Repositories\MindInMotion\DATA';

my_folders = dir(fullfile(my_directory, 'H*'));
my_folders = [my_folders; dir(fullfile(my_directory, 'NH*'))];
my_folders = my_folders(~contains({my_folders.name}, '_')); % Remove *_FU, _Airram, etc
my_files = [];
files = cell(length(my_folders),1);
parfor i = 1:length(my_folders)
    disp(['Pulling: ', my_folders(i).name])
    this_path = fullfile(my_folders(i).folder, my_folders(i).name, ...
        'EEG', 'Trials');
    files_all = dir(fullfile(this_path, '*.set'));
    files_s = files_all(startsWith({files_all.name}, 'SP'));
    files_t = files_all(startsWith({files_all.name}, 'TM'));
    files{i} = [files_s; files_t];
end
my_files = vertcat(files{:});

min_contacts = 65;

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
    dat_start = find(strcmp(dat{:,7}, 'TrialStart'));
    dat_end = find(strcmp(dat{:,7}, 'TrialEnd'));
    dat = dat(dat_start:dat_end, :);

    % Find left and right heel contacts
    contacts_left = find(strcmp(dat{:,7}, 'LHS'));
    contacts_right = find(strcmp(dat{:,7}, 'RHS'));
    contacts = sort([contacts_left; contacts_right], 'ascend');
    contacts_latency = table2array(dat(contacts,1));

    % If minimum number of contacts are not met then add NaN.
    if length(contacts) < min_contacts
        hursts(i,:) = NaN;
        entropies(i,:) = NaN;
    else
        % contacts_latency = contacts_latency(1:min_contacts); % Cut to first # of minimum strides
        contacts_latency = contacts_latency(max(1, numel(contacts_latency) - min_contacts + 1):end); % Cut to last # of minimum strides

        % Calculate step intervals
        step_intervals = diff(contacts_latency)/dat_temp.srate;

        % Hurst Exponent
        hurst = median(bayesH(step_intervals, 200));
        hursts(i,:) = hurst;

        % Entropy
        entropy = Samp_En(step_intervals, 2 , 0.25, std(step_intervals));
        entropies(i,:) = entropy;

        % Plot/save stride intervals to check
        % f = figure('Visible', 'off');
        % plot(step_intervals);
        % % ylim([1, 2]);
        % save_name = sprintf('%s_%s_%s_%s.png', temp_id{6}, temp_id_condition{1}, temp_id_condition{2}, temp_id_condition{3});
        % save_path = fullfile(output_directory, '\FIGURES\', save_name);
        % saveas(f, save_path);
        % close(f);
    end

end

% Export
mim_results = table(id, treadmill, speed, terrain, trial, hursts, entropies, ...
    'VariableNames', {'id', 'treadmill', 'speed', 'terrain', 'trial', 'hurst', 'entropy'});
writetable(mim_results, fullfile(output_directory, 'MIM_Nonlinear_Results.csv'));
