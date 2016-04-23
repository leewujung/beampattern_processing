function data = bp_proc(data,pname,fname,tnum,chk_indiv_call,track_cut_idx)
% Beampattern processing main code
% all functions transplated from beampattern_gui_v6.m
% 
% pname   path to the info matching file
% fname   filename of the info matching file
% tnum    trial number to be processed
% chk_indiv_call  1-plot and check peak detection for each call
%
% Wu-Jung Lee | leewujung@gmail.com
% 2015 10 21  
% 2015 10 23  Combine all mic related info into one function
% 2015 10 27  Change the directory structure so that raw mic data and
%             detection results can be loaded from different directories
% 2015 10 28  Enable reading call_start_idx and call_end_idx in the folder
% 2015 11 12  Plot and check the peak detection for each call

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Manual params
data.path.base_dir = pname;
data.files.match_fname = fname;
data.param.tempC = (data.param.tempF-32)*9/5;  % temperature [deg C]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Set default path
[~,~, data.param.c, data.param.c_iso] = ...    % sounds speed [m/s]
    air_absorption_vec(1e3,data.param.tempC,data.param.humid);
data.path.bat_pos = './bat_pos';
data.path.mic_data = './mic_data';
data.path.mic_detect = './mic_detect';
data.path.mic_info = './mic_info';
data.path.mic_bp = './mic_bp';
data.path.mic_sens = './mic_sens';
data.path.proc_output = './proc_output';
% cd(data.path.base_dir);


%% Load all related files
[~,tt,~] = xlsread(fullfile(data.path.base_dir,data.files.match_fname));
tt(1,:) = [];
data.match_seq = tt;
data.trial_num = tnum;
data.files.bat_pos = data.match_seq{tnum,1};
data.files.mic_data = data.match_seq{tnum,2};
data.files.mic_detect = [strtok(data.files.mic_data,'.'),'_detect.mat'];  % xxx_detect.mat
data.files.mic_info = data.match_seq{tnum,3};
data.files.mic_sens = data.match_seq{tnum,4};
data.files.mic_bp = data.match_seq{tnum,5};

disp('------------------------------------------------------');
disp(['Processing ',data.files.mic_data]);

% Load data and info
disp('Loading all data and related info...');
data = load_mic_info(data);
data = load_mic_bp_sens(data);
data = load_bat_pos(data,track_cut_idx);
data = load_mic_data(data);


%% Gather angle and call data

% Calculate all angle info
disp('Calculating all angle info...');
data = time_delay_btwn_bat_mic(data);
data = mic2bat(data);
data = bat2mic(data);

% Reserve space for GUI checking
data.proc.chk_good_call = zeros((length(data.mic_data.call_idx_w_track)),1);  % set all to bad call
if isempty(data.mic_data.call_idx_w_track)
    disp('No calls with bat position info in this file!');
    return
end
data.proc.ch_ex{length(data.mic_data.call_idx_w_track)} = [];  % channels to be excluded

% Call amplitude calculation and compensation
disp('Calculating call amplitude...');
data = get_time_series_around_call_fcn(data);
data = get_call_fcn(data,chk_indiv_call);
data = compensate_call_dB_fcn(data);

% Remove raw mic signals and long sections from structure
data.proc = rmfield(data.proc,'call_align');
data.proc = rmfield(data.proc,'call_no_align');
data.mic_data = rmfield(data.mic_data,'sig');
data = rmfield(data,'match_seq');


function data = time_delay_btwn_bat_mic(data)
% Find actual call emission time given trajectory and call received time

bat_traj = data.track.track_interp;
mic_loc = data.mic_loc;
bat_traj_time = data.track.track_interp_time;

% Calculation
time_from_bat_to_mic = zeros(size(bat_traj,1),size(mic_loc,1));
for iT=1:size(bat_traj,1)
    for iCH=1:data.mic_data.num_ch_in_file
            time_from_bat_to_mic(iT,iCH) = norm(bat_traj(iT,:)-mic_loc(iCH,:))/data.param.c;
    end
end
if size(bat_traj_time,1)==1
    bat_traj_time = bat_traj_time';
end
time_of_call_at_mic = time_from_bat_to_mic + repmat(bat_traj_time,1,size(time_from_bat_to_mic,2));

data.time_from_bat_to_mic = time_from_bat_to_mic;  % time delay between the bat and the mic
data.time_of_call_at_mic = time_of_call_at_mic;  % received time of call on channel iCH if bat emits call at location idx iT


function data =  mic2bat(data)
% Calculate the azimuth and elevation angle of each mic with respective to
% the bat, and then plot the vector figure and the az-el location figure

proc_call_num = length(data.mic_data.call_idx_w_track);

for iC=1:proc_call_num

% Project the mics onto sphere around bat
curr_call_loc_idx_on_track = data.track.call_loc_idx_on_track_interp(iC);
bat_loc_at_call = data.track.track_interp(curr_call_loc_idx_on_track,:);
bat_loc_at_call_rep = repmat(bat_loc_at_call,size(data.mic_loc,1),1);
mic_to_bat_vec = data.mic_loc - bat_loc_at_call_rep;
mic_to_bat_dist = diag(sqrt(mic_to_bat_vec*mic_to_bat_vec'));
mic_to_bat_vec = mic_to_bat_vec./repmat(mic_to_bat_dist,1,3);

aim_v = data.head_aim.head_aim_int(curr_call_loc_idx_on_track,:);  % head aim at call
norm_v = data.head_normal.head_normal_int(curr_call_loc_idx_on_track,:);  % head normal at call

% Calculate the azimuth and elevation angle of each mic and object
[mic2bat_2d,mic2bat_x] = find_mic_az_el_to_bat_fcn(mic_to_bat_vec,aim_v,norm_v);

% Save data
mic_to_bat_dist = mic_to_bat_dist(:)';  % make this into row vector

data.proc.bat_loc_at_call(iC,:) = bat_loc_at_call(1,:);  % bat location at call emission (interpolated track)
data.proc.mic_to_bat_dist(iC,:) = mic_to_bat_dist;  % distance from bat to mic [m]
data.proc.mic_to_bat_vec(iC,:,:) = mic_to_bat_vec;  % vector direction from bat to mic
data.proc.mic_to_bat_angle(iC,:,:) = mic2bat_2d;  % angle of each mic from bat's perspective [azimuth, elevation] --> 2d config
data.proc.mic_to_bat_angle_x(iC,:,:) = mic2bat_x;  % angle of each mic from bat's perspective [azimuth, elevation] --> cross config
data.proc.source_head_aim(iC) = data.track.marker_indicator(curr_call_loc_idx_on_track);  % source of head aim vector
                                                                                          % 1-from marker
                                                                                          % 0-from smoothed track
end


function data = bat2mic(data)
% Calculate the angle between the mic axis and the vector from mic to bat

proc_call_num = length(data.mic_data.call_idx_w_track);

for iC=1:proc_call_num

bat_loc_at_call_rep = repmat(data.proc.bat_loc_at_call(iC,:),size(data.mic_loc,1),1);

% Find angle from the bat to each mic
bat_to_mic_vec = bat_loc_at_call_rep-data.mic_loc;
bat_to_mic_vec = bat_to_mic_vec./repmat(sqrt(diag(bat_to_mic_vec*bat_to_mic_vec')),1,3);
bat_to_mic_angle = acos(diag(bat_to_mic_vec*data.mic_vec'));

% Save data
data.proc.bat_to_mic_angle(iC,:) = bat_to_mic_angle;
end


function data = load_mic_info(data)
A = load(fullfile(data.path.base_dir,data.path.mic_info,data.files.mic_info));
data.mic_loc = A.mic_loc(:,[3 1 2]);  % permute to get the x-y-z coordinate right
data.mic_vec = A.mic_vec(:,[3 1 2]);
data.mic_vh = A.mic_vh;
data.mic_gain = A.mic_gain;
clear A


function data = load_mic_bp_sens(data)
data.mic_sens = load(fullfile(data.path.base_dir,data.path.mic_sens,data.files.mic_sens));
data.mic_bp = load(fullfile(data.path.base_dir,data.path.mic_bp,data.files.mic_bp));


function data = load_bat_pos(data,cut_idx)
bat = load(fullfile(data.path.base_dir,data.path.bat_pos,data.files.bat_pos));
sm_len = data.track.smooth_len;  % track smoothing length
diff_len = round(data.track.head_aim_est_time_diff*data.track.fs/1e3);  % pt difference for estimating head aim from track

if length(bat.bat_pos)==1  % 1 marker on head/bat position
    % Raw tracks
    pos = cell2mat(bat.bat_pos);
    if isempty(cut_idx)  % extract only the part with acoustic data
        cut_idx = 1:size(pos,1);
    end
    track = pos(cut_idx,:);
    track_t = -fliplr(0:size(track,1)-1)/data.track.fs;  % Time stamp of the track [sec]

    % Find segments and smoothing
    seg_idx = find_seg(track,sm_len);  % continuous segment
    track_sm = sm_track(track(:,[3 1 2]),sm_len,seg_idx);  % smoothing
    
    % Get head aim and head normal
    head_aim = nan(size(track_sm));
    head_aim(1:end-diff_len+1,:) = [track_sm(diff_len:end,1:2)-track_sm(1:end-diff_len+1,1:2),zeros(size(track_sm,1)-diff_len+1,1)];
    head_aim = norm_mtx_vec(head_aim);
    
    head_n = nan(size(track_sm));
    notnanidx_head_n = ~isnan(head_aim(:,1));
    head_n(notnanidx_head_n,:) = repmat(data.track.head_n_prescribe,sum(notnanidx_head_n),1);
    
    % Interpolate to finer resolution for aligning calls
    track_int_t_interval = 1e-3;  % interpolate to 1ms interval
    track_int_t = track_t(1):track_int_t_interval:track_t(end);
    track_int = int_track(track_t,track_sm,track_int_t);
    head_aim_int = int_track(track_t,head_aim,track_int_t);
    head_n_int = int_track(track_t,head_n,track_int_t);
    
    % Indicator of where head aim/normal comes from
    marker_indic = zeros(size(track_int,1),1);  % 0-head aim derived from track
    
    % Save raw and smoothed marker positions
    data.track.marked_pos = pos(:,[3 1 2]);
    
    
elseif length(bat.bat_pos)==3  % 3 markers on the head
    % Raw tracks
    pos = cell2mat(bat.bat_pos);
    pos = reshape(pos,length(pos),3,[]);
    if isempty(cut_idx)  % extract only the part with acoustic data
        cut_idx = 1:size(pos,1);
    end
    pos = pos(cut_idx,:,:);
    track = nanmean(pos,3);  % raw bat track: mean of three points
    % track_all3 = mean(pos,3);  % raw bat track with all 3 points presents
    % track = track(:,[3 1 2]);  % change axis sequence to corresponding the ground reference
    track_t = -fliplr(0:size(track,1)-1)/data.track.fs;  % Time stamp of the track [sec]
    
    % Find segments
    seg_idx = find_seg(track,sm_len);
    % seg_idx_all3 = find_seg(track_all3,sm_len);
    seg_idx_tip = find_seg(pos(:,:,1),sm_len);
    seg_idx_left = find_seg(pos(:,:,2),sm_len);
    seg_idx_right = find_seg(pos(:,:,3),sm_len);
    
    % Smooth track and marker pos
    track_sm = sm_track(track(:,[3 1 2]),sm_len,seg_idx);
    tip = sm_track(pos(:,[3 1 2],1),sm_len,seg_idx_tip);
    left = sm_track(pos(:,[3 1 2],2),sm_len,seg_idx_left);
    right = sm_track(pos(:,[3 1 2],3),sm_len,seg_idx_right);
    
    head_aim = norm_mtx_vec(tip-(left+right)/2);
    head_n = norm_mtx_vec(cross(tip-left,tip-right));
    
    % Interpolate to finer resolution for aligning calls
    track_int_t_interval = 1e-3;  % interpolate to 1ms interval
    track_int_t = track_t(1):track_int_t_interval:track_t(end);
    track_int = int_track(track_t,track_sm,track_int_t);
    head_aim_int = int_track(track_t,head_aim,track_int_t);
    head_n_int = int_track(track_t,head_n,track_int_t);
    
    % Indicator of where head aim/normal comes from
    marker_indic = zeros(size(track_int,1),1);
    marker_indic(~isnan(head_aim_int(:,1))) = 1;  % 1-head aim derived from interpolated marker locations
                                                  % 0-head aim derived from track
    
    % Fake head aim from smoothed track
    head_aim_fake = [track_sm(sm_len:end,1:2)-track_sm(1:end-sm_len+1,1:2),zeros(size(track_sm,1)-sm_len+1,1)];
    head_aim_fake = norm_mtx_vec(head_aim_fake);
    head_aim_fake_int = int_track(track_t(1:end-sm_len+1),head_aim_fake,track_int_t);
    
    % Fake head normal from mics on the floor
    A = data.mic_loc(data.param.mic_floor_idx,:);
    A0 = bsxfun(@minus,A,mean(A,1)); % Subtract "mean" point
    [~,~,V] = svd(A0,0);
    norm_vec = norm_mtx_vec(V(:,3)');
    
    % Fill in gaps for head aim and head normal
    nanidx = isnan(head_aim_int(:,1));
    head_aim_int(nanidx,:) = head_aim_fake_int(nanidx,:);
    head_n_int(nanidx,:) = repmat(norm_vec,sum(nanidx),1);
    
    % Save raw and smoothed marker positions
    data.track.marked_pos = pos(:,[3 1 2],:);
    data.track.tip_smooth = tip;
    data.track.left_smooth = left;
    data.track.right_smooth = right;
end

% Save trajectory
data.track.track_raw = track;
data.track.track_raw_time = track_t;
data.track.track_smooth = track_sm;
data.track.track_interp = track_int;
data.track.track_interp_time = track_int_t;

data.track.marker_indicator = marker_indic;  % Indicator of where head aim/normal comes from

data.head_aim.head_aim_smooth = head_aim;
data.head_aim.head_aim_int = head_aim_int;

data.head_normal.head_normal_smooth = head_n;
data.head_normal.head_normal_int = head_n_int;



function seg_idx = find_seg(pos,sm_len)
% Find continous segment in the trajectories with small gaps filled
notnan = ~isnan(pos(:,1));
idx_nan = find(diff(notnan)~=0)+1;

g = normpdf(-sm_len:sm_len,0,sm_len);
w = conv(double(notnan),g,'same');
idx = find(diff(w~=0)~=0)+1;
idx_up = find(diff(w~=0)>0)+1;
idx_dn = find(diff(w~=0)<0)+1;

if ~isempty(idx)
    [~,iconv] = min(abs(repmat(idx',length(idx_nan),1)-repmat(idx_nan,1,length(idx))),[],1);
    if isempty(idx_up) && ~isempty(idx_dn)
        seg_idx = [1;idx_nan(iconv)];
    elseif ~isempty(idx_up) && isempty(idx_dn)
        seg_idx = [idx_nan(iconv);size(pos,1)];
    else
        if idx_up(1)~=idx(1)  % if the first up edge not at the beginning
            seg_idx = [1;idx_nan(iconv)];
        else
            seg_idx = idx_nan(iconv);
        end
    end
else
    seg_idx = [1,size(pos,1)];
end
if mod(length(seg_idx),2)~=0  % if the last down edge not at the end
    seg_idx = [seg_idx;length(notnan)];
end
seg_idx = reshape(seg_idx,2,[])';


function v_int = int_track(x,v,x_int)
v_int(:,1) = interp1(x,v(:,1),x_int);
v_int(:,2) = interp1(x,v(:,2),x_int);
v_int(:,3) = interp1(x,v(:,3),x_int);


function v_sm = sm_track(v,sm_len,seg_idx)
v_sm = nan(size(v));
for iS=1:size(seg_idx,1)
    v_sm(seg_idx(iS,1):seg_idx(iS,2),1) = smooth(v(seg_idx(iS,1):seg_idx(iS,2),1),sm_len);
    v_sm(seg_idx(iS,1):seg_idx(iS,2),2) = smooth(v(seg_idx(iS,1):seg_idx(iS,2),2),sm_len);
    v_sm(seg_idx(iS,1):seg_idx(iS,2),3) = smooth(v(seg_idx(iS,1):seg_idx(iS,2),3),sm_len);
end

function mtx_v_norm = norm_mtx_vec(mtx_v)
dd = diag(sqrt(mtx_v*mtx_v'));
mtx_v_norm = mtx_v./repmat(dd,1,3);


function data = load_mic_data(data)
A = load(fullfile(data.path.base_dir,data.path.mic_data,data.files.mic_data));
mic_data = load(fullfile(data.path.base_dir,data.path.mic_detect,data.files.mic_detect));
mic_data.sig = A.sig;
mic_data.fs = A.fs;
clear A
mic_data.sig_t = -fliplr(0:size(mic_data.sig,1)-1)/mic_data.fs;  % time stamps for mic signals [sec]
data.mic_data = mic_data;

% Have fs and can calculate these params
data.param.extract_call_len_pt = round(data.param.extract_call_len*1e-3*data.mic_data.fs);  % sample points
data.param.extract_call_len_idx = -round((data.param.extract_call_len_pt+1)/2)+(1:data.param.extract_call_len_pt);

% Get call idx info on selected track segment
data = get_call_on_seg_stuff(data);


function data = get_call_on_seg_stuff(data)
% Get call and idx info

% Adjust invalid index --> these will be deleted later in this function
call_locs_ini = [data.mic_data.call.locs];
call_locs_ini(call_locs_ini<1) = 1;
call_time = data.mic_data.sig_t(call_locs_ini);
% call_time = data.mic_data.sig_t([data.mic_data.call.locs]);

nan_se_call_idx = isnan([data.mic_data.call.call_start_idx])|...  % Delete calls if start/end haven't been marked
                  isnan([data.mic_data.call.call_end_idx]);
[~,track_interp_time_idx] = min(abs(repmat(call_time,length(data.track.track_interp_time),1)-...
                                    repmat(data.track.track_interp_time',1,length(call_time))),[],1);
notnan_track_idx = isnan(data.head_aim.head_aim_int(track_interp_time_idx,1));   % Delete calls without track info
good_call_idx = find(~(nan_se_call_idx(:)|notnan_track_idx(:)));
call_loc_idx_on_track_interp = track_interp_time_idx(good_call_idx);

% Only include calls with enough flanking region (front and back) to be extracted
call_locs = [data.mic_data.call(good_call_idx).locs];
call_idx_extract_len = (call_locs+data.param.extract_call_len_idx(1))>1 &...
                       (call_locs+data.param.extract_call_len_idx(end))<size(data.mic_data.sig,1);
good_call_idx = good_call_idx(call_idx_extract_len);

% Delete calls in which the marked channel is the one without mic locations
ch_sel = [data.mic_data.call(good_call_idx).channel_marked];
nanidx_mic_loc = find(isnan(data.mic_loc(:,1)));
good_call_idx = good_call_idx(~ismember(ch_sel,nanidx_mic_loc));

% call_time(nan_se_call_idx) = [];
% [~,track_interp_time_idx] = min(abs(repmat(call_time,length(data.track.track_interp_time),1)-...
%                                     repmat(data.track.track_interp_time',1,length(call_time))),[],1);
% notnan_track_idx = find(~isnan(data.head_aim.head_aim_int(track_interp_time_idx,1)));
% call_loc_idx_on_track_interp = track_interp_time_idx(notnan_track_idx);

% Save data
data.mic_data.call_idx_w_track = good_call_idx;  % idx of calls in mic_data.call within the selected track
data.track.call_loc_idx_on_track_interp = call_loc_idx_on_track_interp;  % call emission location in terms of idx of interpolated track

