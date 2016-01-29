function data = get_time_series_around_call_fcn(data)
% Cut out section with call and align according to bat-to-mic distances
% 2015 10 21  Version: used outside of GUI

extract_len_pt = data.param.extract_call_len_pt;  % sample points
extract_len_idx = data.param.extract_call_len_idx;  % index

proc_call_num = length(data.mic_data.call_idx_w_track);

for iC = 1:proc_call_num
  % Extract signal ===============
  curr_call_global_idx = data.mic_data.call_idx_w_track(iC);  % idx of current call among all detected calls
  if ~isempty(find(data.mic_data.call(curr_call_global_idx).locs + extract_len_idx<0,1))
    data.proc.call_loc_on_track_interp(iC,:) = nan;
    data.proc.call_receive_time(iC,:) = nan;
    data.proc.call_align(iC,:,:) = zeros(extract_len_pt,data.mic_data.num_ch_in_file);
    data.proc.call_no_align(iC,:,:) = zeros(extract_len_pt,data.mic_data.num_ch_in_file);
    data.proc.call_align_se_idx(iC,:,:) = zeros(data.mic_data.num_ch_in_file,2);
    data.proc.call_no_align_se_idx(iC,:,:) = zeros(data.mic_data.num_ch_in_file,2);
    continue
  end
  call_mtx_simple = data.mic_data.sig(data.mic_data.call(...
    curr_call_global_idx).locs + extract_len_idx,:);
  
  % adjust call detection time to peak of the channel with max signal amplitude
  [max_pk_val,max_pk_idx_in_ch] = max(call_mtx_simple,[],1);
  chch = [(1:data.mic_data.num_ch_in_file)',max_pk_val(:)];
  chch(isnan(data.mic_loc(:,1)),:) = [];  % delete those channels without locations
  chch = sortrows(chch,2);  % sort according to pk value in each channel
  max_ch_idx = chch(end,1);  % take top channel with max pk value
  current_call_detect_idx = max_pk_idx_in_ch(max_ch_idx) + ...
    data.mic_data.call(curr_call_global_idx).locs + extract_len_idx(1) -1;
  current_call_detect_time = data.mic_data.sig_t(current_call_detect_idx);
  
  % find the index of bat trajectory that produced the call (interpolated track)
  [~,call_emission_idx_in_traj] = min(abs(current_call_detect_time - ...
    data.time_of_call_at_mic(:,max_ch_idx)));
  
  % extract and align signals
  call_align = zeros(extract_len_pt,data.mic_data.num_ch_in_file);
  call_no_align = zeros(extract_len_pt,data.mic_data.num_ch_in_file);
  call_align_se_idx = zeros(data.mic_data.num_ch_in_file,2);
  call_no_align_se_idx = zeros(data.mic_data.num_ch_in_file,2);
  for iM=1:data.mic_data.num_ch_in_file  % loop through all mics
    if isnan(data.mic_loc(iM,1))
      call_align(:,iM) = nan(extract_len_pt,1);
      call_align_se_idx(iM,:) = nan(1,2);
      call_no_align(:,iM) = nan(extract_len_pt,1);
      call_no_align_se_idx(iM,:) = nan(1,2);
    else
      call_receive_time = data.time_of_call_at_mic(call_emission_idx_in_traj,iM);
      [~,call_receive_idx] = min(abs(data.mic_data.sig_t - call_receive_time));
      want_idx = call_receive_idx + extract_len_idx([1 end]);
      if want_idx(1)<1
        want_idx = [1 length(extract_len_idx)];
      end
      if want_idx(2)>length(data.mic_data.sig_t)
        want_idx = length(data.mic_data.sig_t)+[-length(extract_len_idx)+1 0];
      end
      call_align(:,iM) = data.mic_data.sig(want_idx(1):want_idx(2), iM);
      call_align_se_idx(iM,:) = want_idx;
      call_no_align(:,iM) = data.mic_data.sig(data.mic_data.call(...
        curr_call_global_idx).locs + extract_len_idx, iM);
      call_no_align_se_idx(iM,:) = data.mic_data.call(curr_call_global_idx).locs...
        + extract_len_idx([1 end]);
    end
  end
  
  
  % Save data ===================
  
  data.proc.call_loc_on_track_interp(iC,:) = call_emission_idx_in_traj;  % actual call emission location on track (interpolated)
  data.proc.call_receive_time(iC,:) = call_receive_time;  % expected call received time for each channel
  
  data.proc.call_align(iC,:,:) = call_align;  % raw time series without correcting for bat-to-mic distance
  data.proc.call_no_align(iC,:,:) = call_no_align;  % time series corrected for bat-to-mic distance
  
  data.proc.call_align_se_idx(iC,:,:) = call_align_se_idx;  % start and end index of extracted section
  data.proc.call_no_align_se_idx(iC,:,:) = call_no_align_se_idx;
  
end

