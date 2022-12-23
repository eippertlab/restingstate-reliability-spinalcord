% 2. Methods
% 2.3. Data preprocessing
% 2.3.1. Preprocessing of physiological data

% In this code, the physiological nouise metrics were calculated for the following
% manuscript:  "" and for the physiological data associated with the following dataset:
% https://openneuro.org/datasets/ds004386 --> see and subject specific
% sub-directories under derivatives/ folder

% NOTE:
% 'raw' data were acquired with a BrainAmp (BrainAmp ExG system; Brain Products GmbH, Gilching, Germany)
% EEGLAB (https://sccn.ucsd.edu/eeglab/index.php, version 2019.0;) with (https://fsl.fmrib.ox.ac.uk/eeglab/fmribplugin/
% were used to process raw data.
% in-house MATLAB scripts by Dr. Ulrike Horn were used to detect R-peaks.


% Merve Kaptan, mkaptan@cbs.mpg.de
% 22.12.2022
%%
clc
clear all
close all
%%
dataDir = '/data/pt_02098/RELIABILITY_FC/Reliability_Spinal_RestingStatefMRI/derivatives';
codeDir = '/data/pt_02098/RELIABILITY_FC/Reliability_Spinal_RestingStatefMRI/derivatives/Code';
savePath = '/data/pt_02098/RELIABILITY_FC/Reliability_Spinal_RestingStatefMRI/derivatives/Physio/';
sessions = {'auto', 'manual'};

cd(dataDir)
subjects = dir('sub-ZS*');
subjects([9,18,30]) = [];  % exclude subjects with ECG recording problems

saveMode = 1;      % save the results or not

addpath(genpath(codeDir))

if ~exist(savePath)
    mkdir(savePath)
    
end

%%
% some settings to calculate breathing period
sr = 1000; %sampling rate
order = 1; % filter order
fn = sr/2;
lo = 0.6/fn; % low pass
hi = 0.01/fn; % high pass



for sub = 1:numel(subjects)
    
    for ses = 1:numel(sessions)
        
        cd(fullfile(dataDir,subjects(sub).name,'physio'))
        
        %load preprocessed physio data
        tmp = load([sessions{ses} '_physio_table.txt']);
        
        % RR is a vector containing RR intervals in seconds.
        RR = diff(find(tmp(:,1)));
        
        SDNN(sub,ses) = std(RR); %take the standard deviation as SD heart period
        IBImean(sub,ses) = mean(RR); % take the mean as MEAN heart period
       
        % calculate breathing period 
        
        raw = tmp(:,2);
        raw = raw - mean(raw); %mean-center
        [b,a] = butter(order,[hi lo], 'bandpass');
        filt =  filtfilt(b, a, raw); %filter
        
        med_filtdata = medfilt1(filt, sr + 1); %median filter
        
        % find pos/neg zero crossings
        respstamp_ind = find(diff(sign(med_filtdata)) == -2);
        respstamp = respstamp_ind/sr;
        
        ibi = diff(respstamp);
        indx = find(ibi < 1);
        respstamp_ind(indx + 1) = [];
        ibi(indx + 1) = [];
        AvgBP(sub,ses) = mean(ibi);
        StdBP(sub,ses) = std(ibi);
        
        clear tmp RR ibi indx
        
        % calculate mean DVARS and mean REFRMS (these are calculated in the script II_III_II_REL_Preprocessing.m, after moco)
        
        cd(fullfile(dataDir,subjects(sub).name,'func'))
        
        tmp = load(['dvars_' sessions{ses} '.txt']);
        dVars(sub,ses) = mean(tmp(2:end));
        clear tmp
        
        tmp = load(['refrmsFALK_' sessions{ses} '.txt']);
        refrms(sub,ses) = mean(tmp(2:end));
        clear tmp
        
        
        
    end
end

%% calculate ICC and bootstrapped CIs

ICCrefrms =  ICC(refrms, '2', '1');
ICCrefrms_CI = bootci(100000,@myICC,refrms);

ICCdVars =  ICC(dVars, '2', '1');
ICCdVars_CI = bootci(100000,@myICC,dVars);

ICCSDNN = ICC(SDNN, '2', '1');
ICCSDNN_CI = bootci(100000,@myICC,SDNN);

ICCIBImean = ICC(IBImean, '2', '1');
ICCIBImean_CI = bootci(100000,@myICC,IBImean);

ICCAvgBP = ICC(AvgBP, '2', '1');
ICCAvgBP_CI = bootci(100000,@myICC,AvgBP);

ICCStdBP = ICC(StdBP, '2', '1');
ICCStdBP_CI = bootci(100000,@myICC,StdBP);

datasets1 = [ICCrefrms, ICCdVars, ICCIBImean, ICCSDNN, ICCAvgBP, ICCStdBP];
datasets2(1,:) = ICCrefrms_CI;
datasets2(2,:) = ICCdVars_CI;
datasets2(3,:) = ICCIBImean_CI;
datasets2(4,:) = ICCSDNN_CI;
datasets2(5,:) = ICCAvgBP_CI;
datasets2(6,:) = ICCStdBP_CI;


if saveMode
    
    save([savePath filesep 'physioRel_N' num2str(numel(subjects)) '.mat'], 'SDNN', 'IBImean',...
        'dVars', 'refrms', 'AvgBP', 'StdBP', 'datasets1','datasets2', 'subjects');
    
    
end

