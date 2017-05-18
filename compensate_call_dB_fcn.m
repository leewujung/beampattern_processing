function data = compensate_call_dB_fcn(data)
% Compensate for various factors for all channel
% 2015 10 21  Version: used outside of GUI
% 2015 10 22  change data structure for extracted short call section from
%             3D mtx to cell to accommodate different call length

proc_call_num = length(data.mic_data.call_idx_w_track);

for iC = 1:proc_call_num

mic_to_bat_dist = data.proc.mic_to_bat_dist(iC,:);     % distance between bat and mic
bat_to_mic_angle = data.proc.bat_to_mic_angle(iC,:);   % angle to compensate for mic beampattern
% call_psd_raw_dB = data.proc.call_psd_raw_dB{iC};  % call spectrum in dB scale
call_psd_raw_dB = cell2mat(data.proc.call_psd_raw_dB(iC,:)')';  % call spectrum in dB scale
call_freq = data.proc.call_freq_vec{iC,1};  % frequency vector of the call spectrum
num_ch = data.mic_data.num_ch_in_file;  % number of channels in file
d0 = 0.1;  % [m] reference distance from bat

% Transmission loss: air absorption and spreading loss
[alpha, alpha_iso,~,~] = ...
    air_absorption_vec(call_freq,data.param.tempC,data.param.humid);  % atmospheric aborption [dB/m]

air_attn_dB = repmat(alpha,length(mic_to_bat_dist),1).*repmat(mic_to_bat_dist'-d0,1,length(alpha));  % [dB]
spreading_loss_dB = repmat(20*log10(mic_to_bat_dist'/d0),1,length(alpha));  % [dB]
TL_dB = (air_attn_dB+spreading_loss_dB)';  % [dB]

% Mic sensitivity
mic_sens_dB = interp1(data.mic_sens.freq,...
                data.mic_sens.sens,...
                call_freq);  % interpolate mic sensitivity vector

% Mic beampattern: have to loop because interp2 is needed but bp is in 3D mtx
bp_compensation = nan(size(call_psd_raw_dB));
if isempty(bp_compensation)  % **call duration marking error**
    % Fake save data
    data.param.d0 = d0;
    data.param.alpha(iC,:) = cell(1,num_ch);
    data.param.alpha_iso(iC,:) = cell(1,num_ch);
    data.proc.air_attn_dB(iC,:) = cell(1,num_ch);
    data.proc.spreading_loss_dB(iC,:) = cell(1,num_ch);
    data.proc.TL_dB(iC,:) = cell(1,num_ch);
    data.proc.mic_bp_compensation_dB(iC,:) = cell(1,num_ch);
    data.proc.mic_sens_dB(iC,:) = cell(1,num_ch);
    data.proc.call_psd_dB_comp_nobp(iC,:) = cell(1,num_ch);
    data.proc.call_psd_dB_comp_withbp(iC,:) = cell(1,num_ch);
    data.proc.call_psd_dB_comp_re20uPa_nobp(iC,:) = cell(1,num_ch);
    data.proc.call_psd_dB_comp_re20uPa_withbp(iC,:) = cell(1,num_ch);
    data.proc.call_p2p_SPL_comp_re20uPa(iC,:) = NaN(1,num_ch);
    continue
end
[X,Y] = meshgrid(data.mic_bp.theta/180*pi,data.mic_bp.freq);  % convert from [deg] to [rad]
for iM=1:length(bat_to_mic_angle)
    [XI,YI] = meshgrid(bat_to_mic_angle(iM),call_freq);
    bp_compensation(:,iM) = interp2(X,Y,data.mic_bp.bp(:,:,iM),XI,YI);  % interpolate mic beampattern
end

% Broadband call SPL and PSD
call_psd_dB_comp_nobp = call_psd_raw_dB + TL_dB - mic_sens_dB;
call_psd_dB_comp_withbp = call_psd_raw_dB + TL_dB - mic_sens_dB - bp_compensation;
gain_dB_fac = repmat(data.mic_gain(:)',size(call_psd_dB_comp_nobp,1),1);
call_psd_dB_comp_re20uPa_nobp = call_psd_dB_comp_nobp + 20*log10(1/20e-6) - gain_dB_fac;
call_psd_dB_comp_re20uPa_withbp = call_psd_dB_comp_withbp + 20*log10(1/20e-6) - gain_dB_fac;

% Call SPL using p2p voltage
ch_data = cell2mat(data.proc.call_align_short(iC,:)');
ch_p2p = max(ch_data,[],2) + ...
  abs(min(cell2mat(data.proc.call_align_short(iC,:)'),[],2));
call_p2p_ch_dB = 20*log10(ch_p2p');

%using pwelch to calc. the max freq in the data
A=pwelch(ch_data',128,120,call_freq,data.mic_data.fs);
[~,max_freq]=max(A); %where the peak should be...
max_freq_idx = sub2ind(size(A), max_freq, 1:size(A,2));

TL_dB_ch = TL_dB(max_freq_idx);
mic_sens_dB_mean_ch = mic_sens_dB(max_freq_idx);
call_p2p_SPL_comp_re20uPa = call_p2p_ch_dB + TL_dB_ch - mic_sens_dB_mean_ch +...
  20*log10(1/20e-6) - data.mic_gain';

% call SPL using RMS
freqs_RMS=data.proc.call_rms_fcenter{iC,1};
TL_dB_RMS_freq = interp1(call_freq,TL_dB,freqs_RMS);
mic_sens_dB_RMS_freq = interp1(call_freq,mic_sens_dB,freqs_RMS);

call_rms_dB = cell2mat(data.proc.call_rms_dB(iC,:)')';
call_RMS_SPL_comp_re20uPa=nan(length(freqs_RMS),size(data.proc.call_rms_fcenter,2));
for iF = 1:length(freqs_RMS)
  call_RMS_SPL_comp_re20uPa(iF,:) = ...
    call_rms_dB(iF,:) + TL_dB_RMS_freq(iF,:) - mic_sens_dB_RMS_freq(iF,:) +...
    20*log10(1/20e-6) - data.mic_gain';
end

% Save data
call_len = size(call_psd_raw_dB,2);
data.param.d0 = d0;
data.param.alpha(iC,:) = num2cell(repmat(alpha',1,call_len)',2);
data.param.alpha_iso(iC,:) = num2cell(repmat(alpha_iso',1,call_len)',2);
data.proc.air_attn_dB(iC,:) = num2cell(air_attn_dB,2);
data.proc.spreading_loss_dB(iC,:) = num2cell(spreading_loss_dB,2);
data.proc.TL_dB(iC,:) = num2cell(TL_dB',2);
data.proc.mic_bp_compensation_dB(iC,:) = num2cell(bp_compensation',2);
data.proc.mic_sens_dB(iC,:) = num2cell(mic_sens_dB',2);
data.proc.call_psd_dB_comp_nobp(iC,:) = num2cell(call_psd_dB_comp_nobp',2);
data.proc.call_psd_dB_comp_withbp(iC,:) = num2cell(call_psd_dB_comp_withbp',2);
data.proc.call_psd_dB_comp_re20uPa_nobp(iC,:) = num2cell(call_psd_dB_comp_re20uPa_nobp',2);
data.proc.call_psd_dB_comp_re20uPa_withbp(iC,:) = num2cell(call_psd_dB_comp_re20uPa_withbp',2);
data.proc.call_p2p_SPL_comp_re20uPa(iC,:) = call_p2p_SPL_comp_re20uPa;
data.proc.call_RMS_SPL_comp_re20uPa(iC,:) =  num2cell(call_RMS_SPL_comp_re20uPa',2);

end