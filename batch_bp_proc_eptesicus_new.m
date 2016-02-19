% 2015 10 23  Process beampattern data and save to folder
% 2016 02 15  Update for 1 marker case

username = getenv('username');
pname = ['C:\Users\',username,'\Dropbox\0_ANALYSIS\bp_processing'];  % base path
fname = 'eptesicus_20150824_file_match.xlsx';  % spreadsheet containing all matching files of different types
trial_to_proc = 11:25;   % row index of the trials to process in the spreadsheet above
chk_indiv_call = 1;     % whether to display the cut-out section for each call/channel
track_cut_idx = 1:800;  % index for frames with acoustic data
save_dir = ['C:\Users\',username,'\Dropbox\0_ANALYSIS\bp_processing\proc_output_eptesicus_new'];  % path to save processing output
if ~exist(save_dir,'dir')
    mkdir(save_dir);
end
load(['C:\Users\',username,'\Dropbox\0_CODE\beampattern_processing\bpf30.mat']);  % filter use only when detecting Rousettus clicks

for tnum = trial_to_proc
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    data.track.fs = 200;   % video frame rate [Hz]
    data.track.smooth_len = 10;  % number of points used to smooth tracks
    data.track.head_aim_est_time_diff = 50; % [ms] time span between points used for estimating head aim from bat position
    data.track.head_n_prescribe = [0,0,1];  % precribed head normal vector (only used in 1-marker case, ignored in 3-marker case)
    
    data.param.tempF = 75;  % temperature [deg F]
    data.param.humid = 50;  % humidity (relative in %)
    data.param.extract_call_len = 10;        % [ms] length of recording section isolated from around a given call
    data.param.tolernace = 2;               % tolerance for call misalignment, make it bigger when mic location uncertainty is large
    data.param.tukeywin_proportion = 0.25;  % proportional of tukeywin for call tapering
    data.param.dura_flag = 1;   % 1-use duration marked in the mic detect file (FM bats)
                                % 0-use default detection stuff (Rousettus)
    
    % Below are only used when data.param.dura_flag=1, such as for Rousettus clicks
    data.param.call_short_len = 0.5;        % [ms] desired length of extracted call
    data.param.call_portion_front = 0.2;    % proportion of extracted call before the peak
    data.param.click_th = 0.1;              % threshold for extracting click ranges
    data.param.click_bpf = bpf30;           % bandpass filter for use in click range estimation

    data.param.mic_floor_idx = [4,5,7,17,24,26,27];  % index of mics on floor, used to estimate head normal vector
                                                     % in missing marker scenarios for 3-marker case
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    data = bp_proc(data,pname,fname,tnum,chk_indiv_call,track_cut_idx);
    ff = [data.files.mic_data,'_bp_proc.mat'];
    save(fullfile(save_dir,ff),'-struct','data');
    
    clear data
end

