% 2015 11 12  Determine good calls from processed files

clear
usrn = getenv('username');
base_path = ['C:\Users\',usrn,'\Dropbox\0_ANALYSIS\bp_processing'];
bat_proc_path = './proc_output';
bat_proc_file = dir(fullfile(base_path,bat_proc_path,'rousettus_20150825_*.mat'));
freq_wanted = 35e3;
goodcall_angle_range = [-60 60];
mic_area_frac = 0.4;

for iF = 1:length(bat_proc_file)

    data = load(fullfile(base_path,bat_proc_path,bat_proc_file(iF).name));
    
    for iC = 1:length(data.proc.chk_good_call)  % go through all processed calls
        call_dB = nan(data.mic_data.num_ch_in_file,1);
        for iM=1:data.mic_data.num_ch_in_file
            freq = data.proc.call_freq_vec{iC,iM};
            [~,fidx] = min(abs(freq-freq_wanted));
            call_dB(iM) = data.proc.call_psd_dB_comp_re20uPa_withbp{iC,iM}(fidx);
        end
        mic_to_bat_angle = squeeze(data.proc.mic_to_bat_angle(iC,:,:));
        if ~isnan(mic_to_bat_angle(:,1))==0  % if all positions are NaN
            data.proc.chk_good_call(iC) = 0;
        else
            data.proc.chk_good_call(iC) = isgoodcall(mic_to_bat_angle(:,1)/pi*180,mic_to_bat_angle(:,2)/pi*180,call_dB,goodcall_angle_range,mic_area_frac);
        end
    end

    save(fullfile(base_path,bat_proc_path,bat_proc_file(iF).name),'-struct','data');  % update chk_bad_call
    
    clear data
end