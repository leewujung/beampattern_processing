% 2017 03 06  Copy checked good call idx to new results

clear
usrn = getenv('username');
base_path = 'F:\Dropbox\0_ANALYSIS\bp_processing';
bat_proc_path_checked = './proc_output_rousettus_checked';
bat_proc_path_new = './proc_output_rousettus_new_checked';
bat_proc_file = dir(fullfile(base_path,bat_proc_path_checked,'rousettus_20150825_*.mat'));
freq_wanted = 35e3;
goodcall_angle_range = [-60 60];
mic_area_frac = 0.4;

for iF = 1:length(bat_proc_file)

    data_checked = load(fullfile(base_path,bat_proc_path_checked,bat_proc_file(iF).name));
    data = load(fullfile(base_path,bat_proc_path_new,bat_proc_file(iF).name));
    
    data.proc.chk_good_call = data_checked.proc.chk_good_call;

    save(fullfile(base_path,bat_proc_path_new,bat_proc_file(iF).name),'-struct','data');  % update chk_bad_call
    
    clear data
end