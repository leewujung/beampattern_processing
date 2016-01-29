% 2015 10 23  Process beampattern data and save to folder

clear;
addpath('F:\repositories\beampattern_processing\')

chk_indiv_call = 0;
track_cut_idx = 1:800;

pname='..\';
fname=['analysis\eptesicus_20150824_file_match.xlsx'];
save_dir=['..\proc_output\'];

% bandpass filter for use in click range estimation, only used when dura_flag=1
load('bpf30.mat');

% trial_to_proc = 26:34; %LB62
% trial_to_proc = 12:25; %LB53 
% trial_to_proc = 1:10; %LB88
trial_to_proc = [1:10 12:25 26:34]; %all
for tnum = trial_to_proc
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  data.track.fs = 200;   % frame rate [Hz]
  data.param.tempF = 75;  % temperature [deg F]
  data.param.humid = 50;  % humidity (relative in %)
  data.param.extract_call_len = 10;  % [ms]
  data.param.call_short_len = 0.5;  % desired length of extracted call [sec]
  data.param.call_portion_front = 0.2;         % proportion of extracted call before the peak
  data.param.tolernace = 2;
  data.param.tukeywin_proportion = 0.25;  % proportional of tukeywin for call tapering
  data.param.dura_flag = 1;   % 1-use duration marked in the mic detect file (FM bats)
  % 0-use default detection stuff (Rousettus)
  data.param.click_th = 0.1;  % threshold for extracting click ranges, only used when dura_flag=1
  data.param.click_bpf = bpf30;     % bandpass filter for use in click range estimation, only used when dura_flag=1
  data.param.mic_floor_idx = [4,5,7,17,24,26,27];  % index of mics on floor, used to estimate head normal vector
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  data = bp_proc(data,pname,fname,tnum,chk_indiv_call,track_cut_idx);
  ff = [data.files.mic_data,'_bp_proc.mat'];
  save(fullfile(save_dir,ff),'-struct','data');
  
  clear data
end

