function plot_bp_multi_freq(handles)
% Plot interpolated beampattern on azimuth-elevation plane

% Load data
data = getappdata(0,'data');
gui_op = getappdata(0,'gui_op');

start_freq = str2double(get(handles.bp_freq_anim_start,'String'))*1e3;
end_freq = str2double(get(handles.bp_freq_anim_end,'String'))*1e3;

freq_skip=5e3;

freqs=start_freq:freq_skip:end_freq;
cols=parula(length(freqs));

% Check for channels to be excluded
if isempty(data.proc.ch_ex{gui_op.current_call_idx})
  ch_ex_manual = [];
else
  ch_ex_manual = data.proc.ch_ex{gui_op.current_call_idx};
end

% Check for mics without location data
ch_good_loc = ~isnan(data.mic_loc(:,1))';

if ~isempty(data.mic_vh) %cross config
  axes(handles.axes_bp_contour);
  cla(handles.axes_bp_contour,'reset');
  
  axes(handles.axes_bp);
  cla(handles.axes_bp,'reset');
  
  mic_to_bat_angle = squeeze(data.proc.mic_to_bat_angle_x(gui_op.current_call_idx,:,:));
else
  axes(handles.axes_bp_contour);
  cla(handles.axes_bp_contour,'reset');
  axesm eckert4;
  framem('fedgecolor',200*ones(1,3)/255,'flonlimit',[-180 180]);
  gridm('gcolor',190*ones(1,3)/255,'glinestyle','-');
  axis off
  hold on
  
  mic_to_bat_angle = squeeze(data.proc.mic_to_bat_angle(gui_op.current_call_idx,:,:));
end

for ff=1:length(freqs)
  freq_wanted = freqs(ff);
  
  call_dB = nan(1,data.mic_data.num_ch_in_file);
  for iM=1:data.mic_data.num_ch_in_file
    if strcmp(gui_op.linlog,'rb_RMS') %use RMS
      freq = data.proc.call_rms_fcenter{gui_op.current_call_idx,iM};
      [~,fidx] = min(abs(freq-freq_wanted));
      call_dB(iM) = data.proc.call_RMS_SPL_comp_re20uPa{gui_op.current_call_idx,iM}(fidx);
    else %use the PSD
      freq = data.proc.call_freq_vec{gui_op.current_call_idx,iM};
      [~,fidx] = min(abs(freq-freq_wanted));
      call_dB(iM) = data.proc.call_psd_dB_comp_re20uPa_withbp{gui_op.current_call_idx,iM}(fidx);
    end
  end
  
  ch_ex_sig = find(isnan(call_dB));  % low quality channel from call extraction function
  
  if ~isempty(data.mic_vh) %cross config
    % Channels to be considered
    notnanidx = ~ismember(1:data.mic_data.num_ch_in_file,union(ch_ex_manual,ch_ex_sig))...
      & ch_good_loc;
    
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
    if ~isempty(el_x)
      pp = polar(el_plot(:,1),el_plot(:,2),'.-');
      set(pp,'markersize',20,'color',cols(ff,:));
      if ff==1
        text(el_x,el_y,num2str(el_plot(:,3)),'color','r');
      end
    else
      polar(0,0,'.-');
    end
    hold on
    title('Elevation');
    % ff = findall(gca,'type','text');
    % t = strtrim(get(ff,'String'));
    % t = t(13:end);
    % for r = 1:length(rho_labels)
    %    set(ff(strcmp(t,rho_labels{r})),'String',rho_labels2{r})
    % end
    
    
    % Plot azimuth polar circle
    axes(handles.axes_bp);
    if ~isempty(az_x)
      pp = polar(az_plot(:,1),az_plot(:,2),'.-');  % flip az +/-
      set(pp,'markersize',20,'color',cols(ff,:));
      if ff==1
        text(az_x,az_y,num2str(mic_az(sort_index)'),'color','r');  % flip az +/-
      end
    else
      polar(0,0,'.-');
    end
    hold on
    title('Azimuth');
    
  else %it's a 2D pattern
    % Interpolation
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
end

if isempty(data.mic_vh)
  mic_num = 1:data.mic_data.num_ch_in_file;
  textm(el/pi*180,az/pi*180,num2str(mic_num(angle_notnanidx)'),...
    'horizontalalignment','center','fontsize',10);
  tightmap
  title('-3 dB contour across freqs')
end

caxis([freqs(1) freqs(end)]);
colormap(cols);
cbar=colorbar('ticks',freqs,'ticklabels',num2str(freqs'./1e3),...
  'Location','south');
cbar.Label.String = 'Freq (kHz)';
cbar_pos = cbar.Position;
cbar.Position = cbar_pos + [0 -.05 0 0];
