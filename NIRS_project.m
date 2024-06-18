%% HOMEWORK 2 - Imaging for neuroscience A.A. 2023-2024

% GROUP 4 %
% Coppée Aurélien:             2128487      aurelienmichaelm.coppee@studenti.unipd.it    
% De Rivo Valentino:           2125791      valentino.derivo@studenti.unipd.it           
% Donfack Mogou Glwadis Linda: 2084261      glwadislinda.donfackmogou@studenti.unipd.it  
% Friscic Tamara:              2130631      tamara.friscic@studenti.unipd.it             
% Massignan Annaluce:          2109290      annaluce.massignan@studenti.unipd.it         
% Rossato Letizia:             2102366      letizia.rossato.1@studenti.unipd.it          
% Tomasi Sara:                 2095500      sara.tomasi.4@studenti.unipd.it              

%%
clear; close all; clc;

%% Utils
% add all the utils you need in the folder where your code is
addpath(genpath(pwd))

%% Setting folders

% create the results folder
output_file_prefx = fullfile(pwd,'Results');

if not(exist(output_file_prefx, 'dir'))
    mkdir(output_file_prefx);
end

%% 1. Plot the 3D configuration (sources, detectors, channels)

load ('CCW1.nirs','-mat');  
nCh = size(SD.MeasList,1)/2; %number of channels

figure(1)
plot3(SD.SrcPos(:,1),SD.SrcPos(:,2),SD.SrcPos(:,3),'.r','MarkerSize',10) %sources x, y, z coordinates
hold on;
plot3(SD.DetPos(:,1),SD.DetPos(:,2),SD.DetPos(:,3),'.b','MarkerSize',10) %detectors x, y, z coordinates
for iCh = 1:nCh
    src = SD.SrcPos(SD.MeasList(iCh,1),:);  
    det = SD.DetPos(SD.MeasList(iCh,2),:);
    plot3([src(1) det(1)],[src(2) det(2)],[src(3) det(3)],'g')  
end
title('3D array configuration')
legend('Source', 'Detector','Channel')

%% 2. Source-detector distance for each channel and histogram

distCh = zeros(nCh,1);
for iCh = 1:nCh
    src = SD.SrcPos(SD.MeasList(iCh,1),:);
    det = SD.DetPos(SD.MeasList(iCh,2),:);
    distCh(iCh) = sqrt(sum((src-det).^2));
end

figure(2);
histogram(distCh,20)   
xlabel('SD distance [mm]')
ylabel('N of channels')
title('Source-detector distances')

%% 3. Identify 'bad' channels

dRange = 0;

SNRrange = 20;
remCh = removeNoisyChannels(d,dRange,SNRrange);

SD.MeasListAct = remCh;

%% (3.) Plot the 3D array highlighting the bad channels with a different color

figure(3);
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
legend('Source', 'Detector', 'Good channel', '', 'Bad channel')

%% (3.) Plot the 3D array keeping in the plot only the good channels

figure(4);
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

%% Check the distances now 

distCh_tmp = zeros(nCh,1);

for iCh = 1:nCh
    if remCh(iCh) == 1
        src = SD.SrcPos(SD.MeasList(iCh,1),:);
        det = SD.DetPos(SD.MeasList(iCh,2),:);
        distCh_tmp(iCh) = sqrt(sum((src-det).^2));    
    end
end

figure(5);
histogram(distCh_tmp(distCh_tmp>0),20)
xlabel('SD distance [mm]')
ylabel('N of channels')

%% 4. PREPROCESSING
% 4a. Conversion to optical density changes

meanValue = mean(d);
dodConv = -log(abs(d)./meanValue);

% 4b. Motion correction:
%% 5. CHOOSE THE BEST MOTION CORRECTION TECHNIQUE

% Plot intensity data at first wavelength to have a look at artifacts
figure(6);
plot(t,d(:,1:nCh))
xlabel('Time [s]')
ylabel('Intensity [A.U.]')
xlim([t(1) t(end)])
title('Wavelength 1')

%% Wavelet motion correction--> there are artifacts but without big changes of the baseline

iqr = 0.5;
% Run wavelet motion correction - necessary only in the first run
%dodWavelet = hmrMotionCorrectWavelet(dodConv,SD,iqr);
%save("dodWavelet.mat", "dodWavelet")
load('dodWavelet.mat')

% Plot original optical density data at first wavelength considering only good channels
dodConvGood = dodConv(:,remCh==1);

figure(7);
plot(t,dodConvGood(:,1:end/2)) 
xlabel('Time [s]')
ylabel('Optical density [A.U.]')
xlim([t(1) t(end)]), ylim([-0.2 0.2])
title('Wavelength 1')

% Plot wavelet corrected optical density data at first wavelength considering only good channels
dodWavGood = dodWavelet(:,remCh==1);

figure(8);
plot(t,dodWavGood(:,1:end/2))
xlabel('Time [s]')
ylabel('Optical density corrected [A.U.]')
xlim([t(1) t(end)]), ylim([-0.2 0.2])
title('Wavelength 1')

% Compare uncorrected vs. wavelet corrected data at each good channel
visualize = 1;
if visualize ==1
    for iCh = 1:nCh
        if remCh(iCh) == 1 % Plot only if it is a good channel
            figure(9);
            plot(t,dodConv(:,iCh))
            hold on;
            plot(t,dodWavelet(:,iCh))
            title(num2str(iCh))
            legend('dod','dod Wavelet corrected')
            xlim([t(1) t(end)]), ylim([-0.2 0.2])
            pause(0.5)
            hold off;
        end
    end
end

%%
% 4c. Band-pass filtering

lowerCutOff = 0.01;
higherCutOff = 0.5;
fs = 1/(t(2)-t(1)); %sampling frequency
dodFilt = hmrBandpassFilt(dodWavelet,fs,lowerCutOff,higherCutOff);

% Plot filtered optical density data at one selected good channel (here channel 1)
dodFiltGood = dodFilt(:,remCh==1);

figure(10);
plot(t,dodWavGood(:,1))
hold on;
plot(t,dodFiltGood(:,1))
xlabel('Time [s]')
ylabel('Optical density [A.U.]')
xlim([t(1) t(end)])
title('Channel 1 - Wavelength 1')
legend('dod Wavelet corrected','band-pass filtered')

%%
% 4d. Average optical density hemodynamic response for each channel
%     and condition with the block average approach

tRange = [-2 36];
sRange = fix(tRange*fs);
tHRF = tRange(1):1/fs:tRange(2)-1/fs;
dodAvg = zeros(length(tHRF),size(dodFilt,2),size(s,2));

for iS = 1:size(s,2) % for each condition
    stimulusTiming = find(s(:,iS)==1);
    ytrial = zeros(length(tHRF),size(dodFilt,2),length(stimulusTiming));
    
    nTrial = 0;
    for iT = 1:length(stimulusTiming)
        if (stimulusTiming(iT)+sRange(1))>=1 && (stimulusTiming(iT)+sRange(2))<=size(dodFilt,1) 
            nTrial = nTrial + 1;
            ytrial(:,:,nTrial) = dodFilt(stimulusTiming(iT)+[sRange(1):sRange(2)],:,:); % extract the trial from the dc data
        end
    end
    
    dodAvg(:,:,iS) = mean(ytrial(:,:,1:nTrial),3);

    % Correct for the baseline
    for ii = 1:size(dodAvg,2) % for each channel 
        foom = mean(dodAvg(1:-sRange(1),ii,iS),1); %baseline = avg of the signal in the -2:0 seconds time range
        dodAvg(:,ii,iS) = dodAvg(:,ii,iS) - foom;  %subtract the baseline from the average
    end
end

% Plot the average hemodynamic response (good channels, 1st wavelet)
visualize = 1;
if visualize == 1
    for iCh = 1:nCh
        if remCh(iCh)==1
            figure(11);
            plot(tHRF,squeeze(dodAvg(:,iCh,1)),'r','LineWidth',2)
            hold on;
            plot(tHRF,squeeze(dodAvg(:,iCh,2)),'b','LineWidth',2)
            legend('Do nothing','Look at the pattern')
            title(num2str(iCh))
            xlabel('Time [s]')
            ylabel('Optical density [A.U.]')
            xlim([tHRF(1) tHRF(end)])
            ylim([-0.03 0.03]) 
            pause(0.5)
            hold off;
        end
    end
end

%% 6. Test all motion correction approaches

%% (6.) Spline motion correction

% Detect motion artifacts in signal
tMotion = 0.5;
tMask = 2; 
SDThresh = 8;
AmpThresh = 0.3; 

tIncMan = ones(length(t),1); 
[tInc,tIncCh] = hmrMotionArtifactByChannel(dodConv, fs, SD, tIncMan, tMotion, tMask, SDThresh, AmpThresh);

% Visualize detected motion artifacts in good channels at first wavelength
dodConvGood = dodConv(:,remCh==1);
figure(12);
plot(t,dodConvGood(:,1:end/2))
hold on;
for i = 1:length(tInc) % here we use tInc since we are displaying all channels and we want to see all time points that have been detected as motion in at least one channel
    if tInc(i) == 0 % For each sample identified as motion artifact
        lh = plot([t(i) t(i)],[-0.5 0.2],'r','LineWidth',0.5); % Plot vertical red line
        lh.Color = [lh.Color 0.05]; % the first three entries define color, the fourth one transparency
    end
end
title('First Wavelength')
xlim([t(1) t(end)])
xlabel('Time [s]')
ylabel('\DeltaOD [A.U.]')

% Motion artifacts plot with intensity data
figure(13);
plot(t,d(:,1:nCh))
hold on;
for i = 1:length(tInc) % here we use tInc since we are displaying all channels and we want to see all time points that have been detected as motion in at least one channel
    if tInc(i) == 0 % For each sample identified as motion artifact
        lh = plot([t(i) t(i)],[0 0.3],'r','LineWidth',0.5); % Plot vertical red line
        lh.Color = [lh.Color 0.05]; % the first three entries define color, the fourth one transparency
    end
end
title('First Wavelength')
xlim([t(1) t(end)])
xlabel('Time [s]')
ylabel('Intensity [A.U.]')

% Spline interpolation
p = 0.99;
dodSpline = hmrMotionCorrectSpline(dodConv,t,SD,tIncCh,p);

% Plot original optical density data in all good channels of first wavelength
figure(14);
plot(t,dodConvGood(:,1:end/2))
xlabel('Time [s]')
ylabel('Optical density [A.U.]')
xlim([t(1) t(end)])
title('Wavelength 1')

% Plot spline corrected optical density data in all good channels of first wavelength
dodSplineGood = dodSpline(remCh == 1);
figure(15);
plot(t,dodSpline(:,1:end/2))
xlabel('Time [s]')
ylabel('Optical density corrected [A.U.]')
xlim([t(1) t(end)])
title('Wavelength 1')

% Compare uncorrected vs. spline corrected data at each good channel with
% superimposed vertical lines for where artifacts were detected
visualize = 1;
if visualize ==1
    for iCh = 1:nCh
        if remCh(iCh) == 1 % display only good channels
            figure(16);
            plot(t,dodConv(:,iCh))
            hold on;
            plot(t,dodSpline(:,iCh))
            title(num2str(iCh))
            xlim([t(1) t(end)]), ylim([-0.2 0.2])
            hold on;
            for i = 1:size(tIncCh,1) 
                if tIncCh(i,iCh) == 0 % here we use tIncCh since we are displaying channels one at a time and we are interested in evaluating whether spline was able to correct the artifacts detected specifically in each single channel
                    lh = plot([t(i) t(i)],[-0.5 0.2],'r','LineWidth',0.5);
                    lh.Color = [lh.Color 0.05]; % the first three entries define color, the fourth one transparency
                end
            end
            legend('dod','dod Spline corrected')
            pause(0.5)
            hold off;            
        end
    end
end

%% (6.) tPCA motion correction
varThresh = 0.97; 
nIter = 5; 

[dodPCA,tIncPCAafter,svs,nSV,tIncPCAbefore] = hmrMotionCorrectPCArecurse(dodConv,fs,SD,tIncMan,tMotion,tMask,SDThresh,AmpThresh,varThresh,nIter);

% Plot original optical density data in all good channels of first wavelength
figure(17);
plot(t,dodConvGood(:,1:end/2))
xlabel('Time [s]')
ylabel('Optical density [A.U.]')
xlim([t(1) t(end)])
title('Wavelength 1')

% Plot tPCA corrected optical density data in all good channels of first wavelength
dodPCAGood = dodPCA(:,remCh == 1);
figure(18);
plot(t,dodPCAGood(:,1:end/2))
xlabel('Time [s]')
ylabel('Optical density corrected [A.U.]')
xlim([t(1) t(end)])
title('Wavelength 1')

% Compare uncorrected vs. tPCA corrected data at each good channel with
% superimposed vertical lines for where artifacts were detected
visualize=1;
if visualize ==1
    for iCh = 1:nCh
        if remCh(iCh) == 1 % only good channels
            figure(19);
            plot(t,dodConv(:,iCh))
            hold on;
            plot(t,dodPCA(:,iCh))
            title(num2str(iCh))
            xlim([t(1) t(end)]), ylim([-0.2 0.2])
            for i = 1:length(tIncPCAbefore) % here we use tIncPCAbefore since we are displaying all channels and we want to know on which segments of contaminated data tPCA worked on (remember that tPCA is a multi-channel approach working on all channels)
                if tIncPCAbefore(i) == 0
                    lh = plot([t(i) t(i)],[-0.5 0.2],'r','LineWidth',0.5);
                    lh.Color = [lh.Color 0.05]; % the first three entries define color, the fourth one transparency
                end
            end
            legend('dod','dod tPCA corrected')
            pause(0.5)
            hold off;
        end
    end
end

%% (6.) Comparison between approaches

for iCh = 1:nCh
    if remCh(iCh) == 1 % display only good channels 
        figure(20);
        plot(t,dodConv(:,iCh))
        hold on;
        plot(t,dodSpline(:,iCh))
        plot(t,dodWavelet(:,iCh))
        plot(t,dodPCA(:,iCh))
        title(num2str(iCh))
        xlim([t(1) t(end)]), ylim([-0.2 0.2])
        legend('dod','dod Spline','dod Wavelet','dod PCA')
        pause(0.8)
        hold off;
    end
end

% farei il plot di 4 canali significativi giusto per far vedere che quella
% che abbiamo scelto è in effetti la migliore


%% 7 Display whole array sensitivity for the first wavelength

% Load data
load(fullfile(pwd, 'MNI','HeadVolumeMesh.mat'))
load(fullfile(pwd, 'MNI','GMSurfaceMesh.mat'))
load(fullfile(pwd, 'MNI','ScalpSurfaceMesh.mat'))
load(fullfile(pwd, 'MNI','TissueMask.mat'))

% Load txt file with 10-5 positions
fid = fopen(fullfile(pwd, 'MNI','10-5_Points.txt'),'r');
tmp = textscan(fid,'%s %f %f %f','CollectOutput',1);
fclose(fid);
tenFive = tmp{2};
tenFiveLabel = tmp{1};

% Load txt file with cranial landmarks coordinates
fid = fopen(fullfile(pwd, 'MNI','LandmarkPoints.txt'),'r');
tmp = textscan(fid,'%s %f %f %f','CollectOutput',1);
fclose(fid);
cranialL = tmp{2};
cranialLLabel = tmp{1};

% Load Jacobian
load('CCW.jac','-mat')

% whole array sensitivity for first wavelength
HeadVolumeMesh.node(:,4) = (sum(J{1}.vol));

% Display whole array sensitivity (1st wavelength) on GM volume mesh
figure(21);
plotmesh(HeadVolumeMesh.node,HeadVolumeMesh.elem(HeadVolumeMesh.elem(:,5)==4,1:4))
caxis([-3 0])
hold on;
plot3(SD.SrcPos(:,1),SD.SrcPos(:,2),SD.SrcPos(:,3),'.r','MarkerSize',10)
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
title('Whole array sensitivity on GM volume mesh with all channels')

% Remove bad channels from Jacobian
for i = 1:length(SD.Lambda)
    tmp = J{i}.vol;
    JCropped{i} = tmp(SD.MeasListAct(SD.MeasList(:,4)==i)==1,:);
end

HeadVolumeMesh.node(:,4) = (sum(JCropped{1}));

% Display whole array sensitivity (1st wavelength) on GM volume mesh with only good channels
figure(22);
plotmesh(HeadVolumeMesh.node,HeadVolumeMesh.elem(HeadVolumeMesh.elem(:,5)==4,1:4))
caxis([-3 0])
hold on;
plot3(SD.SrcPos(:,1),SD.SrcPos(:,2),SD.SrcPos(:,3),'.r','MarkerSize',10)
plot3(SD.DetPos(:,1),SD.DetPos(:,2),SD.DetPos(:,3),'.b','MarkerSize',10)
for iCh = 1:nCh
    src = SD.SrcPos(SD.MeasList(iCh,1),:);
    det = SD.DetPos(SD.MeasList(iCh,2),:);
    if remCh(iCh) == 1
        plot3([src(1) det(1)],[src(2) det(2)],[src(3) det(3)],'g')
    end
end
title('Whole array sensitivity on GM volume mesh - only good channels')


%% 8 Reconstruction of HbO and HbR images

% Compute inverse of Jacobian
lambda1 = 0.1;
invJ = cell(length(SD.Lambda),1);

for i = 1:length(SD.Lambda) %for each Jacobian
    Jtmp = JCropped{i};
    JJT = Jtmp*Jtmp';
    S=svd(JJT);
    invJ{i} = Jtmp'/(JJT + eye(length(JJT))*(lambda1*max(S)));
end


% Data to reconstruct are optical density changes compared to a baseline.
% In our case the baseline is 0, since we have already correct for the baseline
% in section 4d, therefore we want to reconstruct 0-our data
datarecon = -dodAvg;

% Inizialize matrices and load useful stuff
nNodeVol = size(HeadVolumeMesh.node,1);  %The node count of the volume mesh
nNodeGM = size(GMSurfaceMesh.node,1); %The node count of the GM mesh
nFrames = size(datarecon,1); % Number of samples to reconstruct
load('vol2gm.mat')
wavelengths = SD.Lambda; % wavelengths of the system
nWavs = length(SD.Lambda); % n of wavelengths
nCond = size(s,2); % number of condition

% Initialize final results matrices
hbo.vol = zeros(nFrames,nNodeVol,nCond);
hbr.vol = zeros(nFrames,nNodeVol,nCond);
hbo.gm = zeros(nFrames,nNodeGM,nCond);
hbr.gm = zeros(nFrames,nNodeGM,nCond);

% Obtain specific absorption coefficients
Eall = [];
for i = 1:nWavs
    Etmp = GetExtinctions(wavelengths(i));
    Etmp = Etmp(1:2); %HbO and HbR only
    Eall = [Eall; Etmp./1e7]; %This will be nWavs x 2;
end

% For each condition, perform reconstruction
for cond = 1:nCond
    
    % For each frame
    for frame = 1:nFrames
        
        % Reconstruct absorption changes
        muaImageAll = zeros(nWavs,nNodeVol);
        for wav = 1:nWavs
            dataTmp = squeeze(datarecon(frame,SD.MeasList(:,4)==wav & SD.MeasListAct==1,cond));
            invJtmp = invJ{wav};
            tmp = invJtmp * dataTmp';
            muaImageAll(wav,:) = tmp; %This will be nWavs * nNode
        end
        
        % Convert to concentration changes
        hbo_tmpVol = (Eall(2,2)*muaImageAll(1,:) - Eall(1,2)*muaImageAll(2,:))/(Eall(1,1)*Eall(2,2)-Eall(1,2)*Eall(2,1));
        hbr_tmpVol = (muaImageAll(2,:)-Eall(1,2)*hbo_tmpVol)/Eall(2,2);        
        
        % Map to GM surface mesh
        hbo_tmpGM = (vol2gm*hbo_tmpVol');
        hbr_tmpGM = (vol2gm*hbr_tmpVol');
        
        % Book-keeping and saving
        hbo.vol(frame,:,cond) = hbo_tmpVol;
        hbr.vol(frame,:,cond) = hbr_tmpVol;
        hbo.gm(frame,:,cond) = hbo_tmpGM;
        hbr.gm(frame,:,cond) = hbr_tmpGM;
        
    end
end

% Plot reconstructed images at different time points for each condition
tRecon = [0 10 15 20]; % time point 0 10  15  and 20 seconds
baseline = abs(tRange(1)); % two seconds of baseline
sRecon = fix(tRecon*fs)+fix(baseline*fs); % Convert to samples
load greyJet % load colormap to make better images
for iCond = 1:nCond
    for iT = 1:length(sRecon)
        % Assign image to fourth column of node
        GMSurfaceMesh.node(:,4) = hbo.gm(sRecon(iT),:,iCond);
        figure(22+iCond);
        set(gcf, 'Position', get(0, 'Screensize'));
        subplot(2,2,iT),
        plotmesh(GMSurfaceMesh.node,GMSurfaceMesh.face)
        caxis([-0.05 0.05]) % Set the limit of the colorbar
        view([0 90]) % Set the view angle
        title(['HbO cond ' num2str(iCond) ' t = ' num2str(tRecon(iT)) ' s'])
        colormap(greyJet) % set the loaded colormap
        hb = colorbar;
        hb.Label.String = {'\DeltaHbO [\muM]'}; % assign label to colorbar
        axis off % remove axis
        
        GMSurfaceMesh.node(:,4) = hbr.gm(sRecon(iT),:,iCond);
        figure(22+nCond+iCond);
        set(gcf, 'Position', get(0, 'Screensize'));
        subplot(2,2,iT),
        plotmesh(GMSurfaceMesh.node,GMSurfaceMesh.face)
        view([0 90])
        caxis([-0.05 0.05])
        title(['HbR cond ' num2str(iCond) ' t = ' num2str(tRecon(iT)) ' s'])
        colormap(greyJet)
        hb = colorbar;
        hb.Label.String = {'\DeltaHbR [\muM]'};
        axis off
    end
end





