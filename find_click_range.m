function click_idx = find_click_range(sig,sig_max,sig_max_loc,th,numfilt)

sig_filt = filtfilt(numfilt,1,sig);
sm_env = smooth(abs(hilbert(sig_filt)),20);
% [sig_max,sig_max_loc] = max(sm_env);
indic = sm_env>sig_max*th;  % all locations above threshold
indic_keyidx = find(diff(indic)~=0)+1;
if ~isempty(indic_keyidx)
    click_idx(1) = max([max(indic_keyidx(indic_keyidx<sig_max_loc)) NaN]);
    click_idx(2) = min([min(indic_keyidx(indic_keyidx>sig_max_loc)) NaN]);
else
    click_idx = nan(1,2);
end

