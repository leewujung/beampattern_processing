function click_idx = find_click_range(sig,th,numfilt)

sig_filt = filtfilt(numfilt,1,sig);
sm_env = smooth(abs(hilbert(sig_filt)),30);
[mm,mmidx] = max(sm_env);
indic = sm_env>mm*th;  % all locations above threshold
indic_keyidx = find(diff(indic)~=0)+1;
click_idx(1) = max(indic_keyidx(indic_keyidx<mmidx));
click_idx(2) = min(indic_keyidx(indic_keyidx>mmidx));


