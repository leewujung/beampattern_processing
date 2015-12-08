function [call_dB,az,el,ch_include_idx] = get_call_azel_dB_data(data,freq_wanted,call_num)

% Get call info
call_dB = nan(1,data.mic_data.num_ch_in_file);
for iM=1:data.mic_data.num_ch_in_file
    freq = data.proc.call_freq_vec{call_num,iM};
    [~,fidx] = min(abs(freq-freq_wanted));
    call_dB(iM) = data.proc.call_psd_dB_comp_re20uPa_withbp{call_num,iM}(fidx);
end

% Get angle info
mic_to_bat_angle = squeeze(data.proc.mic_to_bat_angle(call_num,:,:));
az = mic_to_bat_angle(:,1);
el = mic_to_bat_angle(:,2);

% Check for channels to be excluded
if isempty(data.proc.ch_ex{call_num})
    ch_ex_manual = [];
else
    ch_ex_manual = data.proc.ch_ex{call_num};
end
ch_ex_sig = find(isnan(call_dB));  % low quality channel from call extraction function

% Check for mics without location data
ch_good_loc = ~isnan(data.mic_loc(:,1))';

% Index of mics to be included
ch_include_idx = ~ismember(1:data.mic_data.num_ch_in_file,union(ch_ex_manual,ch_ex_sig)) & ch_good_loc;

