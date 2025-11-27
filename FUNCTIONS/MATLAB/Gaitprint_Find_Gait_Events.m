function [right_contact_locs, right_toeoff_locs, left_contact_locs, left_toeoff_locs] = Gaitprint_Find_Gait_Events(dat, base_filename)

% INPUTS:
% dat = Entire Gaitprint trial matrix (~48,000 x 321)
% base_filename = trial name (i.e., S001_G01_D01_B01_T01.csv)

% OUTPUT:
% Vector for left/right heel contact/toe off. Each row is the timestamp
% where the event occurred, NOT time in seconds.

%%
% TITLE: Gaitprint_Find_Gait_Events.m
% DATE: February 5, 2024
% AUTHOR: Tyler M. Wiles, MS
% EMAIL: tylerwiles@unomaha.edu

% DESCRIPTION:
% This function cleans Gaitprint trials where there are slight difference
% in the number of gait events between the legs (~ +- 5 strides or so). The
% function will return the cleaned contacts and toe offs. % Sometimes
% missing heel contacts are synthesized by adding average amount of frames
% between contacts to the first contact prior to the missing contact.

% EXAMPLE: Below is an example ofhow to use the function, then use the gait
% events to cut gaitprint data to the last # of strides, then make the gait
% events align correctly with the cut data to be used for other functions
% like spatiotemporal calculation. 
% 
% base_filename = my_files(i).name;
% 
% dat = readmatrix(full_filename); % Read in data
% 
% [right_contact_locs, right_toeoff_locs, left_contact_locs, left_toeoff_locs] = Gaitprint_Find_Gait_Events(dat, base_filename);
% 
% % Fix everyone to their final 174 right strides
% dat = dat(right_contact_locs(length(right_contact_locs)-174):end, :);
% 
% % Onset of the first stride
% contact_threshold = right_contact_locs(length(right_contact_locs)-174);
% 
% % Cut gait events to the onset of first stride
% right_contact_locs(right_contact_locs(:) < contact_threshold+1, :) = [];
% right_toeoff_locs(right_toeoff_locs(:) < contact_threshold+1, :) = [];
% left_contact_locs(left_contact_locs(:) < contact_threshold+1, :) = [];
% left_toeoff_locs(left_toeoff_locs(:) < contact_threshold+1, :) = [];
% 
% % Align time with the data
% right_contact_locs = right_contact_locs - (contact_threshold-1);
% right_toeoff_locs = right_toeoff_locs - (contact_threshold-1);
% left_contact_locs = left_contact_locs - (contact_threshold-1);
% left_toeoff_locs = left_toeoff_locs - (contact_threshold-1);
% 
% % Extracting filenames to create figures
% png_filename = filename_split{1};
% png = append(png_filename, '.png');
% 
% [spattemp_temp] = Spatiotemporal_Calculation(base_filename, dat, right_contact_locs, right_toeoff_locs, left_contact_locs, left_toeoff_locs, 200, png_filename, png, 0, 0);

% Copyright 2023, Tyler M. Wiles

% Redistribution and use of this script, with or without
% modification, is permitted provided this copyright notice,
% the original authors name and the following disclaimer remains.

% DISCLAIMER: It is the user's responsibility to check the code is returning
% the appropriate results before reporting any scientific findings.
% The author will not accept any responsibility for incorrect results
% directly relating to user error.

%% Finding right/left heel contact/toeoff locations

[~, right_contact_locs] = findpeaks(dat(:,245),'MinPeakWidth',70, 'MinPeakDistance',50,'MinPeakHeight',0.5);
[~, right_toeoff_locs] = findpeaks(-dat(:,245),'MinPeakWidth',50);

[~, left_contact_locs] = findpeaks(dat(:,154),'MinPeakWidth',70, 'MinPeakDistance',50, 'MinPeakHeight',0.5);
[~, left_toeoff_locs] = findpeaks(-dat(:,154),'MinPeakWidth',50);

%% Trial specific fixes

if base_filename == "S014_G01_D02_B03_T02.csv"
    right_contact_locs(137) = 33232; % Incorrect heel contact
    right_toeoff_locs = [right_toeoff_locs(1:136); 33124; right_toeoff_locs(137:end)]; % Incorrect toeoff
end

if base_filename == "S015_G01_D01_B01_T03.csv" % Missing left toeoff
    left_toeoff_locs = [left_toeoff_locs(1:135); 31016; left_toeoff_locs(136:end)];
end

if base_filename == "S015_G01_D01_B02_T02.csv" % Missing left toeoff
    left_toeoff_locs = [left_toeoff_locs(1:188); 41523; left_toeoff_locs(189:end)];
end

if base_filename == "S015_G01_D02_B01_T02.csv"
    left_toeoff_locs = [left_toeoff_locs(1:67); 14972; left_toeoff_locs(68:end)]; % Missing toeoff
    left_contact_locs(87) = 19184; % Inaccurate heel contact
    left_toeoff_locs = [left_toeoff_locs(1:193); 43068; left_toeoff_locs(194:end)]; % Missing toeoff
end

if base_filename == "S015_G01_D02_B01_T02.csv"
    left_toeoff_locs = [left_toeoff_locs(1:67); 14972; left_toeoff_locs(68:end)]; % Missing toeoff
    left_contact_locs(87) = 19184; % Inaccurate heel contact
    left_toeoff_locs = [left_toeoff_locs(1:193); 43068; left_toeoff_locs(194:end)]; % Missing toeoff
end

if base_filename == "S016_G01_D01_B01_T01.csv"
    % Missing left heel contact filled by average of preceding strides
    left_contact_locs = [left_contact_locs(1:52); 10812; left_contact_locs(53:end)];
    % Same for left toeoff
    left_toeoff_locs = [left_toeoff_locs(1:53); 10940; left_toeoff_locs(54:end)];
    left_contact_locs = [left_contact_locs(:); 47890]; % Missing final left contact
end

if base_filename == "S016_G01_D01_B03_T03.csv"
    top = left_contact_locs(1:64); mid = left_contact_locs(65:74); bot = left_contact_locs(75:end);
    left_contact_locs = [top; 13417; mid; 15668; bot]; % Fill missing contact 1
    top = left_toeoff_locs(1:65); mid = left_toeoff_locs(66:75); bot = left_toeoff_locs(76:end);
    left_toeoff_locs = [top; 13545; mid; 15796; bot]; % Fill missing contact 2
end

if base_filename == "S019_G03_D01_B01_T02.csv"
    left_toeoff_locs = [left_toeoff_locs(1:210); 41477; left_toeoff_locs(211:end)];
end

if base_filename == "S024_G01_D02_B01_T03.csv" % Extra right toeoff
    right_toeoff_locs(end) = [];
end

if base_filename == "S028_G01_D02_B02_T03.csv" % Missing right toeoff
    right_toeoff_locs = [right_toeoff_locs(1:124); 22753; right_toeoff_locs(125:end)];
end

if base_filename == "S032_G01_D01_B02_T02.csv" % Missing right toeoff
    right_toeoff_locs = [right_toeoff_locs(1:68); 27076; right_toeoff_locs(69:end)];
end

if base_filename == "S038_G01_D02_B01_T02.csv" % Missing left toeoff
    left_toeoff_locs = [left_toeoff_locs(1:199); 41175; left_toeoff_locs(200:end)];
end

if base_filename == "S048_G02_D01_B02_T03.csv"
    right_contact_locs(1) = []; % First contact is false
end

if base_filename == "S049_G02_D01_B01_T01.csv"
    left_contact_locs(1) = []; % False first toeoff/contact
    left_contact_locs(103) = 28220; % Fix timing of odd heel contact
    left_toeoff_locs(1) = [];
end

if base_filename == "S064_G02_D01_B03_T01.csv" % Missing left toeoff
    left_toeoff_locs = [left_toeoff_locs(1:169); 38405; left_toeoff_locs(170:end)];
end

if base_filename == "S080_G02_D01_B01_T02.csv"
    left_contact_locs = [ceil(left_contact_locs(1)-mean(diff(left_contact_locs))); left_contact_locs];
    top = left_toeoff_locs(1); bot = left_toeoff_locs(2:end);
    left_toeoff_locs = [top; mean([top(end) bot(1)]); bot];
end

if base_filename == "S083_G02_D02_B02_T03.csv"
    top = right_contact_locs(1:207); bot = right_contact_locs(208:end);
    right_contact_locs = [top; mean([top(end) bot(1)]); bot];
    top = right_toeoff_locs(1:208); bot = right_toeoff_locs(209:end);
    right_toeoff_locs = [top; mean([top(end) bot(1)]); bot];
end

if base_filename == "S102_G03_D01_B01_T02.csv" % Missing same left contacts/toeoffs
    left_contact_locs = [left_contact_locs(1:69); 15251; left_contact_locs(70:end)];
    left_toeoff_locs = [left_toeoff_locs(1:70); 15386; left_toeoff_locs(71:end)];
    left_contact_locs = [left_contact_locs(1:180); 38901; left_contact_locs(181:end)];
    left_toeoff_locs = [left_toeoff_locs(1:181); 39039; left_toeoff_locs(182:end)];
end

if base_filename == "S102_G03_D01_B02_T01.csv"
    % Missing same left contacts
    left_contact_locs = [left_contact_locs(1:170); 35662; left_contact_locs(171:end)];
    left_contact_locs = [left_contact_locs(1:172); 36080; left_contact_locs(173:end)];
    left_contact_locs = [left_contact_locs(1:180); 37773; left_contact_locs(181:end)];
    left_toeoff_locs = [left_toeoff_locs(1:181); 37914; left_toeoff_locs(182:end)]; % Missing toeoff
    left_contact_locs = [left_contact_locs; 47886]; % Missing final contact
end

if base_filename == "S102_G03_D01_B02_T03.csv" % Missing same left contact/toeoff
    left_contact_locs = [left_contact_locs(1:104); 21644; left_contact_locs(105:end)];
    left_toeoff_locs = [left_toeoff_locs(1:105); 21771; left_toeoff_locs(106:end)];
end

if base_filename == "S102_G03_D01_B03_T02.csv" % Missing left contacts and toeoff
    left_contact_locs = [left_contact_locs(1:48); 10245; left_contact_locs(49:end)];
    left_contact_locs = [left_contact_locs(1:63); 13354; left_contact_locs(64:end)];
    left_toeoff_locs = [left_toeoff_locs(1:49); 10378; left_toeoff_locs(50:end)];
end

if base_filename == "S102_G03_D01_B03_T03.csv" % Missing left contacts and toeoffs
    left_contact_locs = [left_contact_locs(1:158); 32364; left_contact_locs(159:end)];
    left_contact_locs = [left_contact_locs(1:220); 45210; left_contact_locs(221:end)];
    left_contact_locs = [left_contact_locs(1:229); 47096; left_contact_locs(230:end)];
    left_toeoff_locs = [left_toeoff_locs(1:159); 32495; left_toeoff_locs(160:end)];
    left_toeoff_locs = [left_toeoff_locs(1:221); 45341; left_toeoff_locs(222:end)];
    left_toeoff_locs = [left_toeoff_locs(1:230); 47227; left_toeoff_locs(231:end)];
end

if base_filename == "S127_G03_D02_B03_T03.csv" % Missing left contacts and toeoffs
    left_contact_locs = [left_contact_locs(1:196); 40973; left_contact_locs(197:end)];
    left_contact_locs = [left_contact_locs(1:197); 41181; left_contact_locs(198:end)];
    left_contact_locs = [left_contact_locs(1:198); 41389; left_contact_locs(199:end)];
    left_contact_locs = [left_contact_locs(1:199); 41597; left_contact_locs(200:end)];
    left_toeoff_locs = [left_toeoff_locs(1:197); 41101; left_toeoff_locs(198:end)];
    left_toeoff_locs = [left_toeoff_locs(1:198); 41309; left_toeoff_locs(199:end)];
    left_toeoff_locs = [left_toeoff_locs(1:199); 41517; left_toeoff_locs(200:end)];
    left_toeoff_locs = [left_toeoff_locs(1:200); 41725; left_toeoff_locs(201:end)];
end

if base_filename == "S128_G03_D02_B01_T03.csv"
    left_contact_locs = [left_contact_locs(1:166); 31487; left_contact_locs(167:end)];
    left_contact_locs = [left_contact_locs(1:235); 44482; left_contact_locs(236:end)];
    left_contact_locs = [left_contact_locs(1:249); 46981; left_contact_locs(250:end)];
end

if base_filename == "S133_G03_D01_B01_T03.csv"
    left_toeoff_locs = [left_toeoff_locs(2:end)];
end

if base_filename == "S166_G05_D01_B01_T03.csv"
    left_toeoff_locs = [left_toeoff_locs(1:47); 11168; 11400; left_toeoff_locs(48:end)];
end

if base_filename == "S166_G05_D01_B02_T02.csv"
    left_toeoff_locs = [left_toeoff_locs(1:42); 9885; left_toeoff_locs(43:end)];
    left_toeoff_locs = [left_toeoff_locs(1:107); 24889; left_toeoff_locs(108:end)];
    left_toeoff_locs = [left_toeoff_locs(1:117); 27219; left_toeoff_locs(118:end)];
end

if base_filename == "S166_G05_D02_B01_T02.csv"
    right_contact_locs = [right_contact_locs(1:34); 8021; right_contact_locs(35:end)];
    right_contact_locs = [right_contact_locs(1:43); 10029; right_contact_locs(44:end)];
    right_contact_locs = [right_contact_locs(1:198); 45149; right_contact_locs(199:end)];
    left_toeoff_locs = [left_toeoff_locs(1:31); 7173; left_toeoff_locs(32:end)];
    left_toeoff_locs = [left_toeoff_locs(1:151); 34265; left_toeoff_locs(152:end)];
end

if base_filename == "S166_G05_D02_B01_T03.csv"
    left_toeoff_locs = [left_toeoff_locs(1:106); 24379; left_toeoff_locs(107:end)];
    left_toeoff_locs = [left_toeoff_locs(1:144); 32976; left_toeoff_locs(145:end)];
    left_toeoff_locs = [left_toeoff_locs(1:188); 43109; left_toeoff_locs(189:end)];
end

if base_filename == "S166_G05_D02_B03_T01.csv"
    left_toeoff_locs = [left_toeoff_locs(1:179); 40887; left_toeoff_locs(180:end)];
    left_toeoff_locs = [left_toeoff_locs(1:191); 43602; left_toeoff_locs(192:end)];
end

if base_filename == "S171_G02_D01_B02_T02.csv"
    left_toeoff_locs = [left_toeoff_locs(1:152); 33093; left_toeoff_locs(153:end)];
end
