%% HOMEWORK 2 - Imaging for neuroscience A.A. 2023-2024

% GROUP 4 %
% Coppée Aurélien:             2128487      aurelienmichaelm.coppee@studenti.unipd.it    
% De Rivo Valentino:           2125791      valentino.derivo@studenti.unipd.it           
% Donfack Mogou Glwadis Linda: 2084261      glwadislinda.donfackmogou@studenti.unipd.it  
% Friscic Tamara:              2130631      tamara.friscic@studenti.unipd.it             
% Massignan Annaluce:          2109290      annaluce.massignan@studenti.unipd.it         
% Rossato Letizia:             2102366      letizia.rossato.1@studenti.unipd.it          
% Tomasi Sara:                 2095500      sara.tomasi.4@studenti.unipd.it              

%% PATHS TO BE SET

clear
close all
clc

path_folder = pwd;

%% Utils


%% Setting folders

%create the results folder
output_file_prefx = fullfile(path_folder,'Results');

if not(exist(output_file_prefx, 'dir'))
    mkdir(output_file_prefx);
end


%% Plot the 3D configuration (sources, detectors, channels)

load ('CCW1.nirs','-mat');  
nCh = size(SD.MeasList,1)/2; % number of channels (one for each wavelength)

figure(1)
plot3(SD.SrcPos(:,1),SD.SrcPos(:,2),SD.SrcPos(:,3),'.r','MarkerSize',10) % sources x, y, z coordinates
hold on;
plot3(SD.DetPos(:,1),SD.DetPos(:,2),SD.DetPos(:,3),'.b','MarkerSize',10) % detectors x, y, z coordinates
for iCh = 1:nCh
    src = SD.SrcPos(SD.MeasList(iCh,1),:); % specific channel 
    det = SD.DetPos(SD.MeasList(iCh,2),:);
    plot3([src(1) det(1)],[src(2) det(2)],[src(3) det(3)],'g') % x y z of starting point and ending point 
end
title('3D array configuration')
legend('Source', 'Detector')

%% Compute the source-detector distance for each channel and plot all distances with a histogram

distCh = zeros(nCh,1);
for iCh = 1:nCh
    src = SD.SrcPos(SD.MeasList(iCh,1),:);
    det = SD.DetPos(SD.MeasList(iCh,2),:);
    distCh(iCh) = sqrt(sum((src-det).^2));
end

figure;
histogram(distCh,20)   
xlabel('SD distance [mm]')
ylabel('N of channels')

%% Identify 'bad' channels as those channels with signal-to-noise ratio(SNR) lower than 20.

% The output of this step should be a column vector with 0 for channels to be removed and 1 
% for channels to be kept. This vector should be placed in the SD.MeasListAct field.

dRange = 0;

% Since no information about the signal range is provided, we modified the
% removeNoisyChannels function to accept 0 as the dRange argument. If
% dRange is 0, then dRange is set as [-Inf, Inf]. This is done because we
% may decide to set a range.

%{
if dRange == 0
    dRange(1) = -Inf;        
    dRange(2) = Inf;
end
%}

% We also tried to use the range of the lab to check the differences

% dRange = [0.001 2.5];

SNRrange = 20;
remCh = removeNoisyChannels(d,dRange,SNRrange);

SD.MeasListAct = remCh;

%% Plot the 3D array highlighting the bad channels with a different color

figure;
plot3(SD.SrcPos(:,1),SD.SrcPos(:,2),SD.SrcPos(:,3),'.r','MarkerSize',10)
hold on;
plot3(SD.DetPos(:,1),SD.DetPos(:,2),SD.DetPos(:,3),'.b','MarkerSize',10)
for iCh = 1:nCh
        src = SD.SrcPos(SD.MeasList(iCh,1),:);
        det = SD.DetPos(SD.MeasList(iCh,2),:);
    if remCh(iCh) == 1
        plot3([src(1) det(1)],[src(2) det(2)],[src(3) det(3)],'g')
    else
        plot3([src(1) det(1)],[src(2) det(2)],[src(3) det(3)],'y')
    end
end
title('3D array configuration')
legend('Source', 'Detector', 'Good channel', '', 'Bad channel')     % Rivedere

%% Plot the 3D array keeping in the plot only the good channels

figure;
plot3(SD.SrcPos(:,1),SD.SrcPos(:,2),SD.SrcPos(:,3),'.r','MarkerSize',10)
hold on;
plot3(SD.DetPos(:,1),SD.DetPos(:,2),SD.DetPos(:,3),'.b','MarkerSize',10)
for iCh = 1:nCh
     src = SD.SrcPos(SD.MeasList(iCh,1),:);
        det = SD.DetPos(SD.MeasList(iCh,2),:);
    if remCh(iCh) == 1
        plot3([src(1) det(1)],[src(2) det(2)],[src(3) det(3)],'g')
    end
end
title('3D array configuration (good bros only)')
legend('Source', 'Detector')




%% Check the distances now (NOT IN THE ASSIGNMENT)

distCh_tmp = zeros(nCh,1);

for iCh = 1:nCh
    if remCh(iCh) == 1
        src = SD.SrcPos(SD.MeasList(iCh,1),:);
        det = SD.DetPos(SD.MeasList(iCh,2),:);
        distCh_tmp(iCh) = sqrt(sum((src-det).^2));    
    end
end

figure;
histogram(distCh_tmp)   
xlabel('SD distance [mm]')
ylabel('N of channels')

%% PREPROCESSING
% a. Conversion to optical density changes

meanValue = mean(d);
dodConv = -log(abs(d)./meanValue);

% b. Motion correction

% CHOOSE THE BEST MOTION CORRECTION TECHNIQUE (5)
% Number of channels available
nCh = size(SD.MeasList,1)/2;

% Plot intensity data at first wavelength to have a look at artifacts
figure;
plot(t,d(:,1:nCh))
xlabel('Time [s]')
ylabel('Intensity [A.U.]')
xlim([t(1) t(end)])
title('Wavelength 1')

%% Wavelet motion correction--> there are artifacts but without changes of the baseline
SD.MeasListAct = remCh;
iqr = 0.5;
% Run wavelet motion correction
dodWavelet = hmrMotionCorrectWavelet(dodConv,SD,iqr);

% Plot original optical density data at first wavelength considering only good channels
dodConvGood = dodConv(:,remCh==1);
figure;
plot(t,dodConvGood(:,1:end/2)) % just plot first wavelength channels
xlabel('Time [s]')
ylabel('Optical density [A.U.]')
xlim([t(1) t(end)]) % This can be used so that the figure will display exactly these time points and not leave white empty spaces between axis and start/end of the signal
title('Wavelength 1')

% Plot wavelet corrected optical density data at first wavelength considering only good channels
dodWavGood = dodWavelet(:,remCh==1);
figure;
plot(t,dodWavGood(:,1:end/2))
xlabel('Time [s]')
ylabel('Optical density corrected [A.U.]')
xlim([t(1) t(end)])
title('Wavelength 1')

% Compare uncorrected vs. wavelet corrected data at each good channel
for iCh = 1:nCh
    if remCh(iCh) == 1 % Plot only if it is a good channel
        figure;
        plot(t,dodConv(:,iCh))
        hold on;
        plot(t,dodWavelet(:,iCh))
        title(num2str(iCh))
        legend('dod','dod Wavelet corrected')
        xlim([t(1) t(end)])
        pause
        close
    end
end
%% c. Band-pass filtering with cut-off frequency 0.01 and 0.5 Hz-->noise
% removal

lowerCutOff = 0.01;
higherCutOff = 0.5;
fs = 1/(t(2)-t(1)); % compute sampling frequency
dodFilt = hmrBandpassFilt(dodWavelet,fs,lowerCutOff,higherCutOff);

% Plot filtered optical density data at one selected good channel (here channel 1)
dodFiltGood = dodFilt(:,remCh==1);
figure;
plot(t,dodWavGood(:,1))
hold on;
plot(t,dodFiltGood(:,1))
xlabel('Time [s]')
ylabel('Optical density [A.U.]')
xlim([t(1) t(end)])
title('Wavelength 1')
legend('original','band-pass filtered')

%% d. Computation of the average optical density hemodynamic response for each channel
%and condition in a time range of -2 to 36 seconds from stimulus onset with the block
%average approach


tRange = [-2 36]; % range of timimg around stimulus to define a trial
tHRF = tRange(1):1/fs:tRange(2); % time vector for the hemodynamic response (and trials)

dcAvg = zeros(length(tHRF),size(dc,2),size(dc,3),size(s,2));