% 2015 10 23  Process beampattern data and save to folder

username = getenv('username');
pname = ['C:\Users\',username,'\Dropbox\0_ANALYSIS\bp_processing_large_room'];
fname = 'rousettus_20141204_file_match.xlsx';
trial_to_proc = 2:5;
chk_indiv_call = 1;
track_cut_idx = [];
save_dir = ['C:\Users\',username,'\Dropbox\0_ANALYSIS\bp_processing_large_room\proc_output'];
if ~exist(save_dir,'dir')
    mkdir(save_dir);
end

for tnum = trial_to_proc
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    data.track.fs = 250;   % frame rate [Hz]
    data.param.tempF = 75;  % temperature [deg F]
    data.param.humid = 50;  % humidity (relative in %)
    data.param.extract_call_len = 10;  % [ms]
    data.param.call_short_len = 0.5;  % desired length of extracted call [sec]
    data.param.call_portion_front = 0.2;         % proportion of extracted call before the peak
    data.param.tolernace = 3.5;
    data.param.tukeywin_proportion = 0.25;  % proportional of tukeywin for call tapering
    data.param.dura_flag = 0;  % 1-use duration marked in the mic detect file (FM bats)
    % 0-use default detection stuff (Rousettus)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    data = bp_proc(data,pname,fname,tnum,chk_indiv_call,track_cut_idx);
    
    bat_pos_fname = data.files.bat_pos;
    C = strsplit(data.files.bat_pos,'_');
    ff = [data.files.mic_data,'_',C{end},'_bp_proc.mat'];
    save(fullfile(save_dir,ff),'-struct','data');
    
    clear data
end

