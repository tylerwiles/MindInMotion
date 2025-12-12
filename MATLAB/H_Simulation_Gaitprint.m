
% TITLE: H_Simulation_Gaitprint.m
% DATE: 11/27/2025
% AUTHOR: Tyler M. Wiles, PhD
% EMAIL: twiles@ufl.edu

% DESCRIPTION:
% Take the stride intervals from the gaitprint dataset and export them as
% individual files, per trial, to be thrown into the rscript with the same
% name. This script has to be run locally on desktop.

clear; close all; clc;

addpath(genpath('FUNCTIONS'));
addpath(genpath('DATA'));
addpath(genpath('MATLAB'));
my_directory = 'E:\Research\Gaitprint\CLEAN DATA';
output_directory = 'E:\Research\Gaitprint\STRIDE_INTERVALS';
my_files = dir(fullfile(my_directory, 'S*'));
my_files = my_files(~contains({my_files.name}, 'G05'));

ref = readtable("C:\Users\tyler_3r9w1ip\Desktop\Github_Repositories\MindInMotion\DATA\H_Simulation\Gaitprint_trial_characteristics.csv");
ref(:,19) = strcat(ref{:,1}, {'_'}, ref{:,2}, {'_'}, ref{:,3}, {'_'}, ref{:,4}, {'_'}, ref{:,5});

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

        dat = readmatrix(full_filename); % Read in data

        [right_contact_locs, ~, ~, ~] = Gaitprint_Find_Gait_Events(dat, base_filename); % Get right heel contacts

        stride_intervals = diff(right_contact_locs) / 200; % Get stride interval divided by sampling rate

        writematrix(stride_intervals, fullfile(output_directory, base_filename));

    end

end
