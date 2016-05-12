function plot_bat_mic_vector()
% Plot bat track and various vectors on view track GUI

if isappdata(0,'track_gui_handles')
    % load data and GUI handle
    data = getappdata(0,'data');
    gui_op = getappdata(0,'gui_op');
    track_gui_handles = getappdata(0,'track_gui_handles');

    call_loc_on_track_idx_all = data.track.call_loc_idx_on_track_interp;  % call emission location in terms of idx of interpolated track
    bat_loc_at_call = data.proc.bat_loc_at_call(gui_op.current_call_idx,:);
    bat_loc_at_call_rep = repmat(bat_loc_at_call,size(data.mic_loc,1),1);
    mic_to_bat_vec = squeeze(data.proc.mic_to_bat_vec(gui_op.current_call_idx,:,:));
    bat_headaim = data.head_aim.head_aim_int(call_loc_on_track_idx_all(gui_op.current_call_idx),:);
    bat_headnorm = data.head_normal.head_normal_int(call_loc_on_track_idx_all(gui_op.current_call_idx),:);
    bat_rr = cross(bat_headaim,bat_headnorm);
    
    % plot
    axes(track_gui_handles.axes2);
    corder = get(gca,'colororder');  % get line colors
    set(track_gui_handles.h_curr_call,...
        'Xdata',data.track.track_interp(data.track.call_loc_idx_on_track_interp(gui_op.current_call_idx),1),...
        'Ydata',data.track.track_interp(data.track.call_loc_idx_on_track_interp(gui_op.current_call_idx),2),...
        'Zdata',data.track.track_interp(data.track.call_loc_idx_on_track_interp(gui_op.current_call_idx),3));
    % Comment out below because it's too messy to see all the vectors, but
    % it's good for checking purposes
%     if isfield(track_gui_handles,'h_mic_to_bat_v')  % mic to bat vector
%         set(track_gui_handles.h_mic_to_bat_v,...
%             'XData',bat_loc_at_call_rep(:,1),...
%             'YData',bat_loc_at_call_rep(:,2),...
%             'ZData',bat_loc_at_call_rep(:,3),...
%             'UData',mic_to_bat_vec(:,1)/2,...
%             'VData',mic_to_bat_vec(:,2)/2,...
%             'WData',mic_to_bat_vec(:,3)/2);
%     else
%         h_mic_to_bat_v = quiver3(bat_loc_at_call_rep(:,1),bat_loc_at_call_rep(:,2),bat_loc_at_call_rep(:,3),...
%             mic_to_bat_vec(:,1)/2,mic_to_bat_vec(:,2)/2,mic_to_bat_vec(:,3)/2,'color',corder(5,:));
%         track_gui_handles.h_mic_to_bat_v = h_mic_to_bat_v;
%     end
    if isfield(track_gui_handles,'h_bat_aim_v')  % bat head aim
        set(track_gui_handles.h_bat_aim_v,...
            'XData',bat_loc_at_call(1),...
            'YData',bat_loc_at_call(2),...
            'ZData',bat_loc_at_call(3),...
            'UData',bat_headaim(:,1),...
            'VData',bat_headaim(:,2),...
            'WData',bat_headaim(:,3));
    else
        h_bat_aim_v = quiver3(bat_loc_at_call(1),bat_loc_at_call(2),bat_loc_at_call(3),...
            bat_headaim(1)*2,bat_headaim(2)*2,bat_headaim(3),'linewidth',2,'color',corder(2,:));
        track_gui_handles.h_bat_aim_v = h_bat_aim_v;
    end
    if isfield(track_gui_handles,'h_bat_norm_v')  % bat head normal vector
        set(track_gui_handles.h_bat_norm_v,...
            'XData',bat_loc_at_call(1),...
            'YData',bat_loc_at_call(2),...
            'ZData',bat_loc_at_call(3),...
            'UData',bat_headnorm(:,1),...
            'VData',bat_headnorm(:,2),...
            'WData',bat_headnorm(:,3));
    else
        h_bat_norm_v = quiver3(bat_loc_at_call(1),bat_loc_at_call(2),bat_loc_at_call(3),...
            bat_headnorm(1)*2,bat_headnorm(2)*2,bat_headnorm(3),'linewidth',2,'color',corder(3,:));
        track_gui_handles.h_bat_norm_v = h_bat_norm_v;
    end
    if isfield(track_gui_handles,'h_bat_rr_v')  % vector to the right of the bat
        set(track_gui_handles.h_bat_rr_v,...
            'XData',bat_loc_at_call(1),...
            'YData',bat_loc_at_call(2),...
            'ZData',bat_loc_at_call(3),...
            'UData',bat_rr(:,1),...
            'VData',bat_rr(:,2),...
            'WData',bat_rr(:,3));
    else
        h_bat_rr_v = quiver3(bat_loc_at_call(1),bat_loc_at_call(2),bat_loc_at_call(3),...
            bat_rr(1),bat_rr(2),bat_rr(3),'linewidth',2,'color',corder(4,:));
        track_gui_handles.h_bat_rr_v = h_bat_rr_v;
    end

    %removing previous handles for beampattern
    if isfield(track_gui_handles,'h_beam_dir')
      h_beam_dir = track_gui_handles.h_beam_dir;
      delete(h_beam_dir);
    end
    if isfield(track_gui_handles,'h_beam_patt')
      h_beam_patt = track_gui_handles.h_beam_patt;
      delete(h_beam_patt);
    end
    if isfield(track_gui_handles,'h_beam_text')
      h_beam_text = track_gui_handles.h_beam_text;
      delete(h_beam_text);
    end
    %plot beampattern
    if get(track_gui_handles.checkbox_top_view_bp,'value')
      view(2)
      
      mic_num = 1:data.mic_data.num_ch_in_file;
      freq_tag=findobj('tag','edit_bp_freq');
      freq_desired=str2double(get(freq_tag,'String'));
      flatten_deg=30;
      beam_pattern_length=.6;
      bat=bat_loc_at_call;
      plot_mic_nums=1;
      
      call_I = nan(1,data.mic_data.num_ch_in_file);
      rb_log = findobj('tag','rb_log');
      for iM=mic_num
        freq = data.proc.call_freq_vec{gui_op.current_call_idx,iM};
        [~,fidx] = min(abs(freq-freq_desired*1e3));
        if get(rb_log,'value')
          call_I(iM) = data.proc.call_psd_dB_comp_re20uPa_withbp{gui_op.current_call_idx,iM}(fidx);
        else
          call_I(iM) = ...
            10.^(data.proc.call_psd_dB_comp_re20uPa_withbp{gui_op.current_call_idx,iM}(fidx)/10);
        end
      end
      
      ch_ex = data.proc.ch_ex;
      ch_ex_sig = find(isnan(call_I)); % low quality channel from call extraction function
      ch_good_loc = ~isnan(data.mic_loc(:,1))';  % Check for mics without location data
      
      angle_notnanidx = ~ismember(mic_num,ch_ex_sig) & ch_good_loc & ~ismember(mic_num,ch_ex{gui_op.current_call_idx});
      
      
      mic_vec=data.mic_loc(angle_notnanidx,:)-...
        repmat(bat,size(data.mic_loc(angle_notnanidx,:),1),1);
      [az,el] =...
        cart2sph(mic_vec(:,1),mic_vec(:,2),mic_vec(:,3));
      
      indx_flat=find( el <= flatten_deg/180*pi & el >= -flatten_deg/180*pi ); %pull out indexes within 1 degree el in the interpolated data
      angle_midline = az(indx_flat);
            
      vq = call_I(angle_notnanidx);
      vq_norm=vq-max(vq);
      vq2d = vq_norm(indx_flat);
      vq2d = vq2d - min(vq2d);
      norm_vq2d = vq2d'./max(abs(vq2d'));
      
      %from here: http://www.mathworks.com/matlabcentral/fileexchange/35122-gaussian-fit
      [sigma, thdirfit]=gaussfit(angle_midline,vq2d');
      
      [bx,by]=pol2cart(thdirfit,beam_pattern_length*.8);
      h_beam_dir = plot([bat(1) bat(1)+bx],[bat(2) bat(2)+ by],...
            '-','linewidth',2,'color',[55 126 184]./255);
          
      %plot bp
      [~,isort]=sort(angle_midline);
      [Bx,By]=pol2cart(angle_midline(isort),norm_vq2d(isort).*beam_pattern_length);
      h_beam_patt = plot(repmat(bat(1),size(Bx))+Bx,repmat(bat(2),size(By))+By,...
        '-+','linewidth',1,'color',[.4 .4 .4]);
      if plot_mic_nums
        mics_data=find(angle_notnanidx);
        mics_sorted=mics_data(isort);
        h_beam_text = text(repmat(bat(1),size(Bx))+Bx,repmat(bat(2),size(By))+By,...
          num2str(mics_sorted'));
      end
      
      track_gui_handles.h_beam_dir = h_beam_dir;
      track_gui_handles.h_beam_patt = h_beam_patt;
      track_gui_handles.h_beam_text = h_beam_text;
    else
      view(3)
    end
    
    
    setappdata(0,'track_gui_handles',track_gui_handles);
end