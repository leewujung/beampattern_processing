function plot_time_series_in_gui(handles)
% 2015 10 22  Updated to accommodate different data structure

data = getappdata(0,'data');
gui_op = getappdata(0,'gui_op');

call_no_align = zeros(data.param.extract_call_len_pt,data.mic_data.num_ch_in_file);  % length of signal x number of mics
call_align = zeros(data.param.extract_call_len_pt,data.mic_data.num_ch_in_file);  % length of signal x number of mics
for iM=1:length(data.mic_loc)
    if isnan(data.mic_loc(iM,1))
        call_no_align(:,iM) = nan(data.param.extract_call_len_pt,1);
        call_align(:,iM) = nan(data.param.extract_call_len_pt,1);
    else
        idx = data.proc.call_no_align_se_idx(gui_op.current_call_idx,iM,:);
        call_no_align(:,iM) = data.mic_data.sig(idx(1):idx(2),iM);
        idx = data.proc.call_align_se_idx(gui_op.current_call_idx,iM,:);
        call_align(:,iM) = data.mic_data.sig(idx(1):idx(2),iM);
    end
end

shift_gap = max(max(call_no_align));  % vertical shift gaps for display all channels together
shift_gap = ceil(shift_gap/0.01)*0.01;

axes(handles.axes_time_series);
corder = get(gca,'colororder');
plot((0:data.param.extract_call_len_pt-1)/data.mic_data.fs*1e3,...  % raw time series
     call_no_align + repmat((1:length(data.mic_loc))*shift_gap, data.param.extract_call_len_pt, 1),...
     'color',ones(1,3)*190/255);
hold on
plot((0:data.param.extract_call_len_pt-1)/data.mic_data.fs*1e3,...  % shifted time series
     call_align + repmat((1:length(data.mic_loc))*shift_gap, data.param.extract_call_len_pt, 1),...
     'color',corder(1,:));
for iM=1:length(data.mic_loc)
    if ~isnan(data.mic_loc(iM,1))
        want_idx = data.proc.call_align_short_se_idx(gui_op.current_call_idx,iM,:);
        want_idx = want_idx(1):want_idx(2);
        plot((want_idx-1)/data.mic_data.fs*1e3,...  % extracted short call
            call_align(want_idx,iM)+iM*shift_gap,...
            'color',corder(2,:));
    end
end
hold off
set(gca,'ytick',(1:length(data.mic_loc))*shift_gap,...
        'yticklabel',{num2str((1:length(data.mic_loc))')});
ylabel('Mic digitization sequence','fontsize',12);
xlabel('Time (ms)','fontsize',12);
title('Mic recording');
ylim([0 shift_gap*(length(data.mic_loc)+1)]);