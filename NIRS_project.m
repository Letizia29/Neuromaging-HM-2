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

addpath(genpath(pwd))
% Utils_path = ''; % insert your Utils path
% addpath(genpath(Utils_path))

%% Setting folders

% create the results folder
output_file_prefx = fullfile(path_folder,'Results');

if not(exist(output_file_prefx, 'dir'))
    mkdir(output_file_prefx);
end


%% 1. Plot the 3D configuration (sources, detectors, channels)

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

%% 2. Compute the source-detector distance for each channel and plot all distances with a histogram

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
title('Source-detector distances')

%% 3. Identify 'bad' channels as those channels with signal-to-noise ratio(SNR) lower than 20.

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

%% 4. PREPROCESSING
% a. Conversion to optical density changes

meanValue = mean(d);
dodConv = -log(abs(d)./meanValue);

% b. Motion correction:
%% 5. CHOOSE THE BEST MOTION CORRECTION TECHNIQUE
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
% Run wavelet motion correction - necessary only in the first run
dodWavelet = hmrMotionCorrectWavelet(dodConv,SD,iqr);
save("dodWavelet.mat", "dodWavelet")
% load('dodWavelet.mat')

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
        pause(0.5)
        close
    end
end

%%
% c. Band-pass filtering with cut-off frequency 0.01 and 0.5 Hz-->noise
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

%% Convert to concentration changes
% Compute DPF
alpha = 223.3;
beta = 0.05624;
gamma = 0.8493;
delta = -5.723*10^(-7);
epsilon = 0.001245;
csi = -0.9025;

lambda = SD.Lambda; % to get the wavelengths from your probe information
age = 28; % Age of the participant

DPF = zeros(1,length(lambda));
% for iL = 1:length(lambda)
%     DPF(1,iL) = alpha + beta * age^gamma + delta * lambda(iL)^3 + epsilon * lambda(iL)^2 + csi * lambda(iL);
% end
%or
DPF = alpha + beta * age^gamma + delta * lambda.^3 + epsilon * lambda.^2 + csi * lambda;

% Convert to concentration changes
dc = hmrOD2Conc(dodFilt,SD,DPF);

% Plot concentration changes data for the HbO chromophore for all good channels
figure;
plot(t,squeeze(dc(:,1,remCh(1:nCh)==1))) % squeeze can be used to remove the dimension from the matrix that we have selected and therefore set to 1, to make again the matrix 2D
xlabel('Time [s]')
ylabel('Concentration changes [M]')
xlim([t(1) t(end)])
title('HbO')

%%
% d. Computation of the average optical density hemodynamic response for each channel
% and condition in a time range of -2 to 36 seconds from stimulus onset with the block
% average approach

tRange = [-2 36]; % range of timing around stimulus to define a trial
sRange = fix(tRange*fs); % convert the time in seconds to samples
tHRF = tRange(1):1/fs:tRange(2)-1/fs; % time vector for the hemodynamic response (and trials)
dcAvg = zeros(length(tHRF),size(dc,2),size(dc,3),size(s,2)); % initialize the matrix that will contain our average hemodynamic response for each chromophore, channel and condition
for iS = 1:size(s,2) % for each condition
    % Get the timing of stimulus presentation for that condition
    stimulusTiming = find(s(:,iS)==1); 
    % Initialize the matrix that will contain the single trial responses
    % for that condition for all chromophores and channels
    ytrial = zeros(length(tHRF),size(dc,2),size(dc,3),length(stimulusTiming));
    
    nTrial = 0;
    for iT = 1:length(stimulusTiming) % for each stimulus presented (for eacht trial)
        if (stimulusTiming(iT)+sRange(1))>=1 && (stimulusTiming(iT)+sRange(2))<=size(dc,1) % Check that there are enough data pre and post stimulus (this is useful to check that the first stimulus is presented at least 2 seconds after the start of the acquisition and that the last stimulus has at least 36 seconds of data afterwards)
            nTrial = nTrial + 1;
            ytrial(:,:,:,nTrial) = dc(stimulusTiming(iT)+[sRange(1):sRange(2)],:,:); % extract the trial from the dc data
        end
    end
    
    % Average trials (the fourth dimension of the ytrial matrix)
    dcAvg(:,:,:,iS) = mean(ytrial(:,:,:,1:nTrial),4);
    % Correct for the baseline
    for ii = 1:size(dcAvg,3) % for each channel 
        foom = mean(dcAvg(1:-sRange(1),:,ii,iS),1); % compute baseline as average of the signal in the -2:0 seconds time range
        dcAvg(:,:,ii,iS) = dcAvg(:,:,ii,iS) - foom; % subtract the baseline from the average hemodynamic responses
    end
end

% Plot the average hemodynamic response for HbO and HbR at the selected channels
% for iCh = 1:length(nCh)
%     figure;
%     plot(tHRF,squeeze(dcAvg(:,1,chSel(iCh),1)),'r','LineWidth',2)
%     hold on;
%     plot(tHRF,squeeze(dcAvg(:,1,chSel(iCh),2)),'b','LineWidth',2)
%     plot(tHRF,squeeze(dcAvg(:,2,chSel(iCh),1)),'m--','LineWidth',2)
%     plot(tHRF,squeeze(dcAvg(:,2,chSel(iCh),2)),'c--','LineWidth',2)
%     legend('HbO right hand','HbO left hand','HbO feet','HbR right hand','HbR left hand','HbR feet')
%     title(num2str(chSel(iCh)))
%     xlabel('Time [s]')
%     ylabel('\DeltaHb [M]')
%     xlim([tHRF(1) tHRF(end)])
%     ylim([-1e-7 2e-7]) % Limit the y axis so that all plots have the same axis limits and can be compared
% end

%% 6. Test all motion correction approaches on the data and establish whether your choice based on 
% theory (the one at point 5) did provide good qualitative results.
%% Spline motion correction
% Detect motion artifacts in signal
tMotion = 0.5;
tMask = 2; 
SDThresh = 12; 
AmpThresh = 0.5; 

tIncMan = ones(length(t),1); % set it to vectors of ones (this is a vector used to remove manually parts of the data if needed)
% Motion detection technique. tIncCh is a matrix number of samples x twice n of
% channels which contains for each channel (column) the information about
% whether an artifact was present (0s) or not (1s). tInc is a vector which
% contains information on whether at that time sample in any of the channel
% was present an artifact (0s) or not (1s). tInc can therefore be obtained
% from tIncCh by setting to 0 every row that contains at least one 0. 
[tInc,tIncCh] = hmrMotionArtifactByChannel(dodConv, fs, SD, tIncMan, tMotion, tMask, SDThresh, AmpThresh);

% Visualize detected motion artifacts in good channels at first wavelength
dodConvGood = dodConv(:,remCh==1);
figure;
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

% Spline interpolation
p = 0.99;
dodSpline = hmrMotionCorrectSpline(dodConv,t,SD,tIncCh,p);

% Plot original optical density data in all good channels of first wavelength
figure;
plot(t,dodConvGood(:,1:end/2))
xlabel('Time [s]')
ylabel('Optical density [A.U.]')
xlim([t(1) t(end)])
title('Wavelength 1')

% Plot spline corrected optical density data in all good channels of first wavelength
dodSplineGood = dodSpline(remCh == 1);
figure;
plot(t,dodSpline(:,1:end/2))
xlabel('Time [s]')
ylabel('Optical density corrected [A.U.]')
xlim([t(1) t(end)])
title('Wavelength 1')

% Compare uncorrected vs. spline corrected data at each good channel with
% superimposed vertical lines for where artifacts were detected
for iCh = 1:nCh
    if remCh(iCh) == 1 % display only good channels
        figure;
        plot(t,dodConv(:,iCh))
        hold on;
        plot(t,dodSpline(:,iCh))
        title(num2str(iCh))
        %legend('dod','dod Spline corrected')
        xlim([t(1) t(end)])
        hold on;
        for i = 1:size(tIncCh,1) 
            if tIncCh(i,iCh) == 0 % here we use tIncCh since we are displaying channels one at a time and we are interested in evaluating whether spline was able to correct the artifacts detected specifically in each single channel
                lh = plot([t(i) t(i)],[-0.5 0.2],'r','LineWidth',0.5);
                lh.Color = [lh.Color 0.05]; % the first three entries define color, the fourth one transparency
            end
        end
        pause
        close
    end
end

%% tPCA motion correction
varThresh = 0.97; % % of variance to remove
nIter = 5; % n of iterations

% Looking at the help of the function, we know that one of the outputs
% (tIncPCAbefore) contains information about detected motion artifacts
% before applying the PCA
[dodPCA,tIncPCAafter,svs,nSV,tIncPCAbefore] = hmrMotionCorrectPCArecurse(dodConv,fs,SD,tIncMan,tMotion,tMask,SDThresh,AmpThresh,varThresh,nIter);

% Plot original optical density data in all good channels of first wavelength
figure;
plot(t,dodConvGood(:,1:end/2))
xlabel('Time [s]')
ylabel('Optical density [A.U.]')
xlim([t(1) t(end)])
title('Wavelength 1')

% Plot tPCA corrected optical density data in all good channels of first wavelength
dodPCAGood = dodPCA(:,remCh == 1);
figure;
plot(t,dodPCAGood(:,1:end/2))
xlabel('Time [s]')
ylabel('Optical density corrected [A.U.]')
xlim([t(1) t(end)])
title('Wavelength 1')

% Compare uncorrected vs. tPCA corrected data at each good channel with
% superimposed vertical lines for where artifacts were detected
for iCh = 1:nCh
    if remCh(iCh) == 1 % only good channels
        figure;
        plot(t,dodConv(:,iCh))
        hold on;
        plot(t,dodPCA(:,iCh))
        title(num2str(iCh))
        %legend('dod','dod tPCA corrected')
        xlim([t(1) t(end)])
        hold on;
        for i = 1:length(tIncPCAbefore) % here we use tIncPCAbefore since we are displaying all channels and we want to know on which segments of contaminated data tPCA worked on (remember that tPCA is a multi-channel approach working on all channels)
            if tIncPCAbefore(i) == 0
                lh = plot([t(i) t(i)],[-0.5 0.2],'r','LineWidth',0.5);
                lh.Color = [lh.Color 0.05]; % the first three entries define color, the fourth one transparency
            end
        end
        pause
        close
    end
end

%% Comparison between approaches
for iCh = 1:nCh
    if remCh(iCh) == 1 % display only good channels 
        figure;
        plot(t,dodConv(:,iCh))
        hold on;
        plot(t,dodSpline(:,iCh))
        plot(t,dodWavelet(:,iCh))
        plot(t,dodPCA(:,iCh))
        title(num2str(iCh))
        xlim([t(1) t(end)])
        hold on;
        for i = 1:length(tInc) 
            if tIncCh(i,iCh) == 0 % Display artifacts identified in each channel so that it is easier to evaluate the differences among approaches by looking at specific segments of data
                lh = plot([t(i) t(i)],[-0.5 0.2],'r','LineWidth',0.5);
                lh.Color = [lh.Color 0.05]; % the first three entries define color, the fourth one transparency
            end
        end
        legend('dod','dod Spline','dod Wavelet','dod PCA')
        pause
        close
    end
end


%% Part 7
% Load Jacobian

load('CCW.jac','-mat')
HeadVolumeMesh.node(:,4) = (sum(J{1}.vol));
% Display whole array sensitivity on GM volume mesh
figure;
plotmesh(HeadVolumeMesh.node,HeadVolumeMesh.elem(HeadVolumeMesh.elem(:,5)==4,1:4))
caxis([-3 0])

% Remove bad channels from Jacobian
for i = 1:length(SD.Lambda) % for the first wavelength only
    tmp = J{1}.vol;
    JCropped{i} = tmp(SD.MeasListAct(SD.MeasList(:,4)==i)==1,:);
end

HeadVolumeMesh.node(:,4) = (sum(JCropped{1}));
% Display whole array sensitivity on GM volume mesh with only good channels
figure;
plotmesh(HeadVolumeMesh.node,HeadVolumeMesh.elem(HeadVolumeMesh.elem(:,5)==4,1:4))
caxis([-3 0])


%% 8
% Compute inverse of Jacobian
lambda1 = 0.1;
invJ = cell(length(SD.Lambda),1);
for i = 1:length(SD.Lambda) %for each Jacobian
    Jtmp = JCropped{i};
    JJT = Jtmp*Jtmp';
    S=svd(JJT);
    invJ{i} = Jtmp'/(JJT + eye(length(JJT))*(lambda1*max(S)));
end
    
    % Average trials (the fourth dimension of the ytrial matrix)
    dodAvg(:,:,iS) = mean(ytrial(:,:,1:nTrial),3);
    % Correct for the baseline
    for ii = 1:size(dodAvg,2) % for each channel (but you can do it directly without the for cycle)
        foom = mean(dodAvg(1:-sRange(1),ii,iS),1); % compute baseline as average of the signal in the -2:0 seconds time range
        dodAvg(:,ii,iS) = dodAvg(:,ii,iS) - foom; % subtract the baseline from the average hemodynamic responses
    end

% Data to reconstruct are optical density changes compared to a baseline.
% In our case the baseline is 0, therefore we want to reconstruct 0-our
% data
datarecon = -dodAvg;

% Inizialize matrices and load useful stuff
nNodeVol = size(HeadVolumeMesh.node,1);  %The node count of the volume mesh
nNodeGM = size(GMSurfaceMesh.node,1); %The node count of the GM mesh
nFrames = size(datarecon,1); % Number of samples to reconstruct
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
        figure;
        plotmesh(GMSurfaceMesh.node,GMSurfaceMesh.face)
        caxis([-0.05 0.05]) % Set the limit of the colorbar
        view([0 90]) % Set the view angle
        title(['HbO cond ' num2str(iCond) ' t = ' num2str(tRecon(iT)) ' s'])
        colormap(greyJet) % set the loaded colormap
        hb = colorbar;
        hb.Label.String = {'\DeltaHbO [\muM]'}; % assign label to colorbar
        axis off % remove axis
        
        GMSurfaceMesh.node(:,4) = hbr.gm(sRecon(iT),:,iCond);
        figure;
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




