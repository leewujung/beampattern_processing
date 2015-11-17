% 2015 11 12  Determine good calls from processed files

clear
usrn = getenv('username');
base_path = ['C:\Users\',usrn,'\Dropbox\0_ANALYSIS\bp_processing'];
bat_proc_path = './proc_output';
bat_proc_file = dir(fullfile(base_path,bat_proc_path,'rousettus_20150825_*.mat'));
freq = 35e3;
goodcall_angle_range = [-60 60];
mic_area_frac = 0.4;

for iF = 1:length(bat_proc_file)

    data = load(fullfile(base_path,bat_proc_path,bat_proc_file(iF).name));
    
    for iC = 1:length(data.proc.chk_good_call)  % go through all processed calls
        [~,fidx] = min(abs(freq-data.proc.call_freq_vec{iC}));
        call_dB = squeeze(data.proc.call_psd_dB_comp_re20uPa_withbp{iC}(fidx,:));
        mic_to_bat_angle = squeeze(data.proc.mic_to_bat_angle(iC,:,:));
        data.proc.chk_good_call(iC) = isgoodcall(mic_to_bat_angle(:,1)/pi*180,mic_to_bat_angle(:,2)/pi*180,call_dB,goodcall_angle_range,mic_area_frac);
    end

    save(fullfile(base_path,bat_proc_path,bat_proc_file(iF).name),'-struct','data');  % update chk_bad_call
    
    clear data
end