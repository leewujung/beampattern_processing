function plot_bp_multi_freq(handles)
% Plot interpolated beampattern on azimuth-elevation plane

% Load data
data = getappdata(0,'data');
gui_op = getappdata(0,'gui_op');
mic_to_bat_angle = squeeze(data.proc.mic_to_bat_angle(gui_op.current_call_idx,:,:));

start_freq = str2double(get(handles.bp_freq_anim_start,'String'))*1e3;
end_freq = str2double(get(handles.bp_freq_anim_end,'String'))*1e3;

freq_skip=5e3;

freqs=start_freq:freq_skip:end_freq;
cols=parula(length(freqs));

axes(handles.axes_bp_contour);
cla(handles.axes_bp_contour,'reset');
axesm eckert4;
framem('fedgecolor',200*ones(1,3)/255,'flonlimit',[-180 180]);
gridm('gcolor',190*ones(1,3)/255,'glinestyle','-');
axis off
hold on

% Check for channels to be excluded
if isempty(data.proc.ch_ex{gui_op.current_call_idx})
  ch_ex_manual = [];
else
  ch_ex_manual = data.proc.ch_ex{gui_op.current_call_idx};
end

% Check for mics without location data
ch_good_loc = ~isnan(data.mic_loc(:,1))';

for ff=1:length(freqs)
  freq_wanted = freqs(ff);
  
  call_dB = nan(1,data.mic_data.num_ch_in_file);
  for iM=1:data.mic_data.num_ch_in_file
    freq = data.proc.call_freq_vec{gui_op.current_call_idx,iM};
    [~,fidx] = min(abs(freq-freq_wanted));
    call_dB(iM) = data.proc.call_psd_dB_comp_re20uPa_withbp{gui_op.current_call_idx,iM}(fidx);
  end
  
  ch_ex_sig = find(isnan(call_dB));  % low quality channel from call extraction function
    
  % Interpolation
  mic_num = 1:data.mic_data.num_ch_in_file;
  angle_notnanidx = ~ismember(1:data.mic_data.num_ch_in_file,...
    union(ch_ex_manual,ch_ex_sig)) & ch_good_loc;
  az = mic_to_bat_angle(angle_notnanidx,1);
  el = mic_to_bat_angle(angle_notnanidx,2);
  
  maxref = max(call_dB(angle_notnanidx));
  [azq,elq] = meshgrid(min(az):pi/180:max(az),min(el):pi/180:max(el));
  if strcmp(gui_op.interp,'rb_natural')  % natural neighbor interpolation
    vq = griddata(az,el,call_dB(angle_notnanidx),azq,elq,'natural');
  elseif strcmp(gui_op.interp,'rb_rbf')  % radial basis function interpolation
    vq = rbfinterp([azq(:)';elq(:)'],rbfcreate([az(:)';el(:)'],...
      call_dB(angle_notnanidx),'RBFFunction','multiquadrics'));
    vq = reshape(vq,size(azq));
  end
  vq_norm = vq-maxref;
  
  % Find indices within measured polygon
  k = boundary(az,el,0);  % outer boundary of all measured points
  [in,on] = inpolygon(azq,elq,az(k),el(k));
  in_smpl_poly = in|on;
  clear in on
  vq(~in_smpl_poly) = NaN;
  vq_norm(~in_smpl_poly) = NaN;
  
  if abs(freq_wanted - str2double(get(handles.edit_bp_freq,'String'))*1e3)...
      <freq_skip
    % Update peak call amplitude
    set(handles.edit_peak_db,'String',sprintf('%2.1f',maxref));
  end
  
  contourm(elq/pi*180,azq/pi*180,vq_norm,-3,'linewidth',2,'color',cols(ff,:));
end
textm(el/pi*180,az/pi*180,num2str(mic_num(angle_notnanidx)'),...
  'horizontalalignment','center','fontsize',10);
tightmap
title('-3 dB contour across freqs')
hold off
caxis([freqs(1) freqs(end)]);
colormap(cols);
colorbar('ticks',freqs,'ticklabels',num2str(freqs'./1e3),...
  'Location','southoutside');
