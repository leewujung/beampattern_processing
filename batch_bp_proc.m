% 2015 10 23  Process beampattern data and save to folder
clear

pname = 'E:\Desktop\bp_proc_main_folder';
fname = 'SPECIES_DATE_file_match.xlsx';

trial_to_proc = 1;
chk_indiv_call = 0;
track_cut_idx = 1:400;
save_dir = [pname '\proc_output'];
if ~exist(save_dir,'dir')
    mkdir(save_dir);
end
load('bpf30.mat');

for tnum = trial_to_proc
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    data=[];
    
    data.track.fs = 100;   % frame rate [Hz]
    data.track.smooth_len = 10;  % number of points used to smooth tracks
    data.track.head_aim_est_time_diff = 50; % [ms] time span between points used for estimating head aim from bat position
    
    data.param.head_aim_prescribed=0;%head aim is prescribed (using value for head_aim_prescribe and head_n_prescribe)
    
    data.track.head_aim_prescribe = [0,-1,0];  % precribed head aim vector (only used in 1-marker case, ignored in 3-marker case)
    data.track.head_n_prescribe = [0 0 1];% precribed head normal vector (only used in 1-marker case, ignored in 3-marker case)
    
    
    data.param.tempF = 75;  % temperature [deg F]
    data.param.humid = 50;  % humidity (relative in %)
    data.param.extract_call_len = 10;  % [ms]
    data.param.call_short_len = 0.5;  % desired length of extracted call [sec]
    data.param.call_portion_front = 0.2;         % proportion of extracted call before the peak
    data.param.tolerance = 2;
    data.param.tukeywin_proportion = 0.25;  % proportional of tukeywin for call tapering
    data.param.dura_flag = 1;   % 1-use duration marked in the mic detect file (FM bats)
                                % 0-use default detection stuff (Rousettus)
    data.param.axis_orient= [1 2 3]; %axis orientation ([3 1 2] for small space beampatter 2015 data)
    data.param.zero_bat2mic_angle=0;%mics are pointed at bat at all times
    
    % Below are only used when data.param.dura_flag=0, such as for Rousettus clicks
    data.param.call_short_len = 0.5;        % [ms] desired length of extracted call
    data.param.call_portion_front = 0.2;    % proportion of extracted call before the peak
    data.param.click_th = 0.1;              % threshold for extracting click ranges
    data.param.click_bpf = bpf30;           % bandpass filter for use in click range estimation
    
    
    data.param.PSD_type='pwelch'; %power spectral density calculation
%     data.param.PSD_type='FFT';
    
    data.param.mic_floor_idx = [4,5,7,17,24,26,27];  % index of mics on floor, used to estimate head normal vector, should be in counterclockwise order according to the direction "up"
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    data = bp_proc(data,pname,fname,tnum,chk_indiv_call,track_cut_idx);
    ff = [data.files.mic_data,'_bp_proc.mat'];
    save(fullfile(save_dir,ff),'-struct','data');
end

