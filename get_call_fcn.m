function data = get_call_fcn(data,varargin)
% Extract short section around call
% 
% dura_flag  whether to use marked call start/end index or not
%
% Wu-Jung Lee | leewujung@gmail.com
%
% 2015 10 21  Version: used outside of GUI
% 2015 10 22  change data structure for extracted short call section from
%             3D mtx to cell to accommodate different call length
% 2015 10 28  incorporate reading in call duration
% 2015 11 12  modify plot_opt

if isfield(data.param,'dura_flag')  % duration mark flag set
    dura_flag = data.param.dura_flag;
else
    dura_flag = 0;  % no duration mark flag
end

if nargin==2
    plot_opt = varargin{1};
else
    plot_opt = 0;
end

if plot_opt
    fig_chk = figure('position',[300 50 550 700]);
    corder = get(gca,'colororder');
end

proc_call_num = length(data.mic_data.call_idx_w_track);  % number of processed calls


for iC = 1:proc_call_num

    % Load params
    tolerance = data.param.tolernace;  % tolerance for finding calls around peak of max channel
    tukeywin_prop = data.param.tukeywin_proportion;  % tukey window taper porportion
    num_ch = data.mic_data.num_ch_in_file;  % number of channels in file
    call_long = squeeze(data.proc.call_align(iC,:,:));
    if dura_flag==1
        sidx_in_long = call_sidx-data.proc.call_align_se_idx(iC,ch_sel,1)+1;  % call start idx in extracted portion
        eidx_in_long = call_eidx-data.proc.call_align_se_idx(iC,ch_sel,1)+1;  % call end idx in extracted portion
        mark_ch = data.mic_data.call_idx_w_track(iC);  % channel selected when call was marked
        call_template = call_long(sidx_in_long:eidx_in_long,mark_ch);  % curved out template for call
        call_template_len_pt = length(call_template);
    else
        call_len = data.param.call_short_len;
        call_portion_front = data.param.call_portion_front;
        call_len_pt = round(call_len*1e-3*data.mic_data.fs);
        call_len_idx = -round(call_len_pt*call_portion_front)+(1:call_len_pt);
    end
    
    if dura_flag==1  % Use duration marked in the mic detect file ================================
        % Xcorr to find call in all channel
        ch_xcorr = nan(size(call_long,1)*2-1,num_ch);
        ch_xcorr_env = nan(size(call_long,1)*2-1,num_ch);
        for iM = 1:num_ch
            if ~isnan(data.mic_loc(iM,1))  % if mic location available
                [ch_xcorr(:,iM),xcorr_lags] = xcorr(call_long(:,iM),call_template);
                ch_xcorr_env(:,iM) = abs(hilbert(ch_xcorr(:,iM)));
            end
        end
        [mm,mm_idx] = max(ch_xcorr_env,[],1);  % max of each channel
        [maxenv_top,max_ch_idx] = max(mm);  % max of all channels
        
        % Find tolerated section
        tol_len_pt_half = round(tolerance*1e-3*data.mic_data.fs/2);
        chk_se_idx = mm_idx(max_ch_idx) + tol_len_pt_half*[-1 1];  % index of section to be checked for a peak
        if chk_se_idx(1)<=0
            chk_se_idx(1) = 1;
        end
        if chk_se_idx(2)>size(call_long,1)*2-1
            chk_se_idx(2) = size(call_long,1)*2-1;
        end
        chk_se_idx = chk_se_idx(1):chk_se_idx(2);
        
        % Xcorr across all channels to extract call section
        call_short = nan(call_template_len_pt,num_ch);
        call_short_se_idx = nan(num_ch,2);
        for iM=1:num_ch
            if ~isnan(data.mic_loc(iM,1))  % if mic location available
                maxenv = max(ch_xcorr_env(:,iM));
                if maxenv>maxenv_top*0.5  % if strong signal consider secondary arrival
                    [~,ch_xcorr_pk_idx_tmp] = findpeaks(ch_xcorr_env(chk_se_idx,iM),'SortStr','descend','MinPeakDistance',50,'MinPeakHeight',maxenv*0.5);
                    ch_xcorr_pk_idx_tmp = min(ch_xcorr_pk_idx_tmp);  % take the first arrival in case there is stronger echo
                else  % if very weak signal don't consider second arrival
                    [~,ch_xcorr_pk_idx_tmp] = findpeaks(ch_xcorr_env(chk_se_idx,iM),'SortStr','descend','MinPeakDistance',50,'NPeak',1);
                end
                ch_xcorr_pk_idx_tmp = ch_xcorr_pk_idx_tmp+chk_se_idx(1)-1;
                want_idx = xcorr_lags(ch_xcorr_pk_idx_tmp(1))+[0,call_template_len_pt-1];
                if want_idx(1)<1
                    want_idx(1) = 1;
                end
                if want_idx(2)>size(call_long,1)
                    want_idx(2)=size(call_long,1);
                end
                call_short_se_idx(iM,:) = want_idx;
                want_idx = want_idx(1):want_idx(2);
                call_short(:,iM) = [call_long(want_idx,iM);zeros(call_template_len_pt-length(want_idx),1)];
            end
        end
        
    else  % Use default detection stuff ================================
        % Find location idx in max channel
        pk_loc_ch = [zeros(data.mic_data.num_ch_in_file,2),(1:data.mic_data.num_ch_in_file)'];
        [pk_loc_ch(:,1),pk_loc_ch(:,2)] = max(call_long,[],1);  % peak and loc of each channel
        pk_loc_ch(isnan(pk_loc_ch)) = 0;
        pk_loc_ch = flipud(sortrows(pk_loc_ch,1));
        [~,aa_idx] = min(pk_loc_ch(pk_loc_ch(:,1)>=pk_loc_ch(1,1)*0.8,2));  % find the earliest of all loud arrivals
        max_ch_idx = pk_loc_ch(aa_idx,3);  % index of max channel
        pk_idx = pk_loc_ch(aa_idx,2);  % index of peak in the max channel

        % Section used as template
        want_idx = pk_idx + call_len_idx([1 end]);
        if want_idx(1)<1
            want_idx(1) = 1;
        end
        if want_idx(2)>size(call_long,1)
            want_idx(2) = length(call_max_ch_env);
        end
        want_idx = want_idx(1):want_idx(2);
        
        % Use call in max channel as xcorr template
        call_template = call_long(want_idx,max_ch_idx);  % call in max channel as xcorr template
        [~,call_template_pk_shift] = max(smooth(call_template.^2,50));  % shift need to be added after xcorr
        
        % Xcorr max ch and tolrated shift section
        [ch_xcorr_max,ch_lags] = xcorr(call_long(:,max_ch_idx),call_template);
        [maxenv_top,max_ch_pk_idx] = max(abs(hilbert(ch_xcorr_max)));
        tol_len_pt_half = round(tolerance*1e-3*data.mic_data.fs/2);
        chk_se_idx = max_ch_pk_idx + tol_len_pt_half*[-1 1];  % index of section to be checked for a peak
        if chk_se_idx(1)<=0
            chk_se_idx(1) = 1;
        end
        if chk_se_idx(2)>size(call_long,1)*2-1
            chk_se_idx(2) = size(call_long,1)*2-1;
        end
        chk_se_idx = chk_se_idx(1):chk_se_idx(2);
        
        % Xcorr across all channels
        ch_xcorr = zeros(size(call_long,1)*2-1,num_ch);
        ch_xcorr_env = zeros(size(call_long,1)*2-1,num_ch);
        ch_xcorr_pk_idx = zeros(size(call_long,2),1);
        ch_mm = zeros(size(call_long,2),1);
        for iM=1:num_ch  % going through all channels
            if isnan(data.mic_loc(iM,1))  % if mic location not available
                ch_xcorr_pk_idx(iM) = NaN;
            else
                [ch_xcorr(:,iM),~] = xcorr(call_long(:,iM),call_template);
                ch_xcorr_env(:,iM) = abs(hilbert(ch_xcorr(:,iM)));
                maxenv = max(ch_xcorr_env(:,iM));
                if maxenv>maxenv_top*0.5  % if strong signal consider secondary arrival
                    [ch_mm_tmp,ch_xcorr_pk_idx_tmp] = findpeaks(ch_xcorr_env(chk_se_idx,iM),'SortStr','descend','MinPeakDistance',50,'MinPeakHeight',maxenv*0.5);
                    ch_xcorr_pk_idx_tmp = min(ch_xcorr_pk_idx_tmp);  % take the first arrival in case there is stronger echo
                else  % if very weak signal don't consider second arrival
                    [ch_mm_tmp,ch_xcorr_pk_idx_tmp] = findpeaks(ch_xcorr_env(chk_se_idx,iM),'SortStr','descend','MinPeakDistance',50,'NPeak',1);
                end
                ch_xcorr_pk_idx_tmp = ch_xcorr_pk_idx_tmp+chk_se_idx(1)-1;
                if isempty(ch_xcorr_pk_idx_tmp)
                    ch_xcorr_pk_idx(iM) = NaN;
                    ch_mm(iM) = NaN;
                    fprintf('%s: Problem in call#%d:%d channel#%d\n',datestr(now,'HH:MM AM'),iC,data.mic_data.call_idx_w_track(iC),iM);
                else
                    ch_xcorr_pk_idx(iM) = ch_xcorr_pk_idx_tmp(1);
                    ch_mm(iM) = ch_mm_tmp(1);
                end
            end
        end
        

        if plot_opt
            % find shift_gap between channels
            shift_gap = max(max(ch_xcorr_env));  % vertical shift gaps for display all channels together
            shift_gap = max([floor(shift_gap/0.1)*0.1 0.1]);
            shift_gap = shift_gap/2;
            tstamp = (0:size(ch_xcorr_env,1)-1)/data.mic_data.fs*1e3;
            
            % plot
            figure(fig_chk);
            cla
            plot(tstamp,ch_xcorr_env+repmat((1:num_ch)*shift_gap,length(ch_xcorr_env),1),'color',corder(1,:));
            hold on
            notnanidx = ~isnan(ch_xcorr_pk_idx);
            plot(tstamp(ch_xcorr_pk_idx(notnanidx)),ch_mm(notnanidx)+(find(notnanidx))*shift_gap,'r.','markersize',10);
            plot(tstamp(chk_se_idx(1))*[1 1],[0,shift_gap*(num_ch+2)],'k');  % boundary of peak detection
            plot(tstamp(chk_se_idx(end))*[1 1],[0,shift_gap*(num_ch+2)],'k');
            if any(~notnanidx)  % highlight channels failing peak detection within boundary
                for iN = find(~notnanidx)'
                    plot(tstamp,ch_xcorr_env(:,iN)+iN*shift_gap,'linewidth',2,'color',corder(2,:));
                end
            end
            title(sprintf('Call#%d:%d on track',iC,data.mic_data.call_idx_w_track(iC)));
            ylim([0,shift_gap*(num_ch+2)]);
            set(gca,'ytick',(1:num_ch)*shift_gap,'yticklabel',1:num_ch);
            xlabel('Time (ms)');
            ylabel('Channel number');
            pause(0.5)
            hold off
        end
        
        % Extract call according to xcorr peak location
        call_short = nan(call_len_pt,num_ch);
        call_short_se_idx = nan(num_ch,2);
        for iM=1:num_ch
            if isnan(ch_xcorr_pk_idx(iM))
                call_short(:,iM) = nan(call_len_pt,1);
            else
                match_idx = max([ch_lags(ch_xcorr_pk_idx(iM))+call_template_pk_shift 1]);
                click_idx = find_click_range(call_long(:,iM),call_long(match_idx,iM),match_idx,...
                                             data.param.click_th,data.param.click_bpf);
                click_idx(1) = max([click_idx(1) match_idx+call_len_idx(1)-1]);
                if isnan(click_idx(2))||click_idx(2)>match_idx+call_len_idx(end)-1
                    click_idx(2) = match_idx+call_len_idx(end)-1;
                end
                if click_idx(1)<1
                    click_idx(1) = 1;
                end
                if click_idx(2)>size(call_long,1)
                    click_idx(2)=size(call_long,1);
                end
                call_short_se_idx(iM,:) = click_idx;
                want_idx = click_idx(1):click_idx(2);
                call_short(:,iM) = [call_long(want_idx,iM);zeros(call_len_pt-length(want_idx),1)];
            end
        end
        
    end

% Calculate fft ================================
w = tukeywin(size(call_short,1),tukeywin_prop);
call_short_taper = call_short.*repmat(w,1,size(call_short,2));  % taper call for spectrum estimation
call_fft = fft(call_short_taper);
call_freq_vec = linspace(0,data.mic_data.fs/2,round((size(call_fft,1)+1)/2));
call_fft_len = length(call_freq_vec);
% call_psd = 2*abs(call_fft(1:call_fft_len,:)).^2/(length(call_fft)*data.mic_data.fs);
call_psd = 2*abs(call_fft(1:call_fft_len,:)).^2;
call_psd_dB = 10*log10(call_psd);


% Save data ====================================
data.proc.call_align_short(iC,:) = num2cell(call_short',2);
data.proc.call_align_short_se_idx(iC,:,:) = call_short_se_idx;
data.proc.call_fft(iC,:) = num2cell(call_fft',2);  % call spectrum
data.proc.call_freq_vec(iC,:) = num2cell(repmat(call_freq_vec',1,size(call_short,2))',2);  % frequency vector for call spectrum
data.proc.call_psd_raw_linear(iC,:) = num2cell(call_psd',2);  % spectrum of extracted calls, linear scale
data.proc.call_psd_raw_dB(iC,:) = num2cell(call_psd_dB',2);   % spectrum of extracted calls, dB scale


clearvars -except proc_call_num plot_opt fig_chk corder dura_flag data

end

if plot_opt
    close(fig_chk);
end
