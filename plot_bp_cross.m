function plot_bp_cross(handles,freq_wanted)
% Plot interpolated beampattern in polar plot for cross configuration
% Interpolation superimposed on actual measurements

% Load data
data = getappdata(0,'data');
gui_op = getappdata(0,'gui_op');

if isempty(data.mic_vh)
    msgbox('This is 2D config only data');
    gui_op.mic_config = 'rb_2d';
    setappdata(0,'gui_op',gui_op);
    set(handles.config_radio_grp,'SelectedObject',handles.rb_2d);
    return
end

mic_to_bat_angle = squeeze(data.proc.mic_to_bat_angle_x(gui_op.current_call_idx,:,:));
if nargin == 1 || isempty(freq_wanted) %allowing optional passing in of freq
  freq_wanted = str2double(get(handles.edit_bp_freq,'String'))*1e3;  % beampattern frequency [Hz]
end

call_dB = nan(1,data.mic_data.num_ch_in_file);
for iM=1:data.mic_data.num_ch_in_file
    freq = data.proc.call_freq_vec{gui_op.current_call_idx,iM};
    [~,fidx] = min(abs(freq-freq_wanted));
    call_dB(iM) = data.proc.call_psd_dB_comp_re20uPa_withbp{gui_op.current_call_idx,iM}(fidx);
end


% Check for channels to be excluded
if isempty(data.proc.ch_ex{gui_op.current_call_idx})
    ch_ex_manual = [];
else
    ch_ex_manual = data.proc.ch_ex{gui_op.current_call_idx};
end
ch_ex_sig = find(isnan(call_dB));  % low quality channel from call extraction function

% Check for mics without location data
ch_good_loc = ~isnan(data.mic_loc(:,1))';

% Channels to be considered
notnanidx = ~ismember(1:data.mic_data.num_ch_in_file,union(ch_ex_manual,ch_ex_sig)) & ch_good_loc;

% Azimuth and elevation stuff
az_idx = ceil(data.mic_vh(:))==1 & notnanidx(:);
el_idx = floor(data.mic_vh(:))==0 & notnanidx(:);
az = mic_to_bat_angle(az_idx,1);  % azimuth
el = mic_to_bat_angle(el_idx,2);  % elevation
if strcmp(gui_op.linlog,'rb_lin')  % log/linear display
    az_amp = 10.^(call_dB(az_idx)/10);
    az_amp = az_amp/min(az_amp);
    el_amp = 10.^(call_dB(el_idx)/10);
    el_amp = el_amp/min(el_amp);
else
    az_amp = call_dB(az_idx);
    el_amp = call_dB(el_idx);
    az_amp = az_amp-max(az_amp)+30;
    az_amp(az_amp<0) = NaN;
    el_amp = el_amp-max(el_amp)+30;
    el_amp(el_amp<0) = NaN;
end

az_plot = [az,az_amp(:),find(az_idx)];
[az_plot,sort_index] = sortrows(az_plot,1);
el_plot = [el,el_amp(:),find(el_idx)];
el_plot = sortrows(el_plot,1);

mic_num_seq = 1:length(az_idx);
if sum(el_idx)
    mic_el = mic_num_seq(el_idx);
    [el_x,el_y] = pol2cart(el_plot(:,1),el_plot(:,2));
else
    mic_el = [];
    el_x = [];
    el_y = [];
end
if sum(az_idx)
    mic_az = mic_num_seq(az_idx);
    [az_x,az_y] = pol2cart(az_plot(:,1),az_plot(:,2));
else
    mic_az = [];
    az_x = [];
    az_y = [];
end

% Plot elevation polar circle
axes(handles.axes_bp_contour);
cla(handles.axes_bp_contour,'reset');
if ~isempty(el_x)
    pp = polar(el_plot(:,1),el_plot(:,2),'.-');
    set(pp,'markersize',20);
    hold on
    text(el_x,el_y,num2str(el_plot(:,3)),'color','r');
else
    polar(0,0,'.-');
end
title('Elevation');
% ff = findall(gca,'type','text');
% t = strtrim(get(ff,'String'));
% t = t(13:end);
% for r = 1:length(rho_labels)
%    set(ff(strcmp(t,rho_labels{r})),'String',rho_labels2{r})
% end


% Plot azimuth polar circle
axes(handles.axes_bp);
cla(handles.axes_bp,'reset');
if ~isempty(az_x)
    pp = polar(az_plot(:,1),az_plot(:,2),'.-');  % flip az +/-
    set(pp,'markersize',20);
    hold on
    text(az_x,az_y,num2str(mic_az(sort_index)'),'color','r');  % flip az +/-
else
    polar(0,0,'.-');
end
title('Azimuth');




