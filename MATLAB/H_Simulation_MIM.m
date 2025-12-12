
% TITLE: H_Simulation_MIM.m
% DATE: 11/27/2025
% AUTHOR: Tyler M. Wiles, PhD
% EMAIL: twiles@ufl.edu

% DESCRIPTION:
% Take the stride intervals from the mind in motion dataset and export them
% as individual files, per trial, to be thrown into the rscript with the
% same name. This script has to be run at school because of EEGLab.

clear; close all; clc;

addpath(genpath('FUNCTIONS'));
addpath(genpath('DATA'));
addpath(genpath('MATLAB'));
addpath(genpath('R:\Ferris-Lab\share\MindInMotion\eeglab'));
eeglab % Load EEGLAB
my_directory = 'R:\Ferris-Lab\share\MindInMotion\Data';
output_directory = 'R:\Ferris-Lab\twiles\H_Simulation\MIM';

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

parfor i = 1:length(my_files)

    % Getting filenames that will be useful for the rest of the code
    base_filename = my_files(i).name;
    full_filename = fullfile(my_files(i).folder, '\', base_filename);

    % Extracting filenames with conditions
    filename_split = strsplit(base_filename, '_');
    temp_id = strsplit(full_filename, '\');
    temp_id_condition = strsplit(base_filename, '_');

    id = temp_id{6};
    treadmill = temp_id_condition{1};
    speed = temp_id_condition{2};
    terrain = 'flat';
    trial = temp_id_condition{3}(1);
    output_name = [id, '_', treadmill, '_', speed, '_', terrain, '_', trial, '.csv'];

    dat_temp = pop_loadset(full_filename); % Load in the EEGLab structure
    dat = struct2table(dat_temp.event); % Just get time series of gait events

    % Cut timeseries to start and end of trial
    dat_start = find(strcmp(dat.type, 'TrialStart'));
    dat_end = find(strcmp(dat.type, 'TrialEnd'));
    dat = dat(dat_start:dat_end, :);

    % Find left and right heel contacts
    contacts_right = find(strcmp(dat.type, 'RHS'));
    contacts_latency = table2array(dat(contacts_right, 1));

    % Calculate step or stride intervals
    stride_intervals = diff(contacts_latency) / dat_temp.srate;

    writematrix(stride_intervals, fullfile(output_directory, output_name));

end
