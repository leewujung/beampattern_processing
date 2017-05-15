function varargout = edit_call_section(varargin)
% EDIT_CALL_SECTION MATLAB code for edit_call_section.fig
%      EDIT_CALL_SECTION, by itself, creates a new EDIT_CALL_SECTION or raises the existing
%      singleton*.
%
%      H = EDIT_CALL_SECTION returns the handle to a new EDIT_CALL_SECTION or the handle to
%      the existing singleton*.
%
%      EDIT_CALL_SECTION('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in EDIT_CALL_SECTION.M with the given input arguments.
%
%      EDIT_CALL_SECTION('Property','Value',...) creates a new EDIT_CALL_SECTION or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before edit_call_section_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to edit_call_section_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help edit_call_section

% Last Modified by GUIDE v2.5 17-Nov-2015 15:51:10

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @edit_call_section_OpeningFcn, ...
                   'gui_OutputFcn',  @edit_call_section_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before edit_call_section is made visible.
function edit_call_section_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to edit_call_section (see VARARGIN)

% Choose default command line output for edit_call_section
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% Set default structure
if ~isappdata(0,'gui_call_op');
    gui_call_op.curr_ch = 1;
    setappdata(0,'gui_call_op',gui_call_op);
else
    gui_call_op = getappdata(0,'gui_call_op');
end

% Load call_align
data = getappdata(0,'data');
gui_op = getappdata(0,'gui_op');
call_align = zeros(data.param.extract_call_len_pt,data.mic_data.num_ch_in_file);  % length of signal x number of mics
for iM=1:length(data.mic_loc)
    if isnan(data.mic_loc(iM,1))
        call_align(:,iM) = nan(data.param.extract_call_len_pt,1);
    else
        idx = data.proc.call_align_se_idx(gui_op.current_call_idx,iM,:);
        call_align(:,iM) = data.mic_data.sig(idx(1):idx(2),iM);
    end
end
gui_call_op.call_align = call_align;
gui_call_op.se_idx = squeeze(data.proc.call_align_short_se_idx(gui_op.current_call_idx,:,:));
gui_call_op.fs = data.mic_data.fs;
gui_call_op.time = (0:size(call_align,1)-1)/data.mic_data.fs;
gui_call_op.num_ch_in_file = data.mic_data.num_ch_in_file;
setappdata(0,'gui_call_op',gui_call_op);

plot_everything(handles,gui_op.current_call_idx);

% Set figure handle to appdata
if ~isappdata(0,'edit_call_gui_handles')
    setappdata(0,'edit_call_gui_handles',handles);
end

% UIWAIT makes edit_call_section wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = edit_call_section_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in button_start.
function button_start_Callback(hObject, eventdata, handles)
% hObject    handle to button_start (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
data = getappdata(0,'data');
gui_call_op = getappdata(0,'gui_call_op');
gui_op = getappdata(0,'gui_op');
hh_main = getappdata(0,'bp_check_gui_handles');

% move start line
[xsel,~] = ginput(1);   % time
gui_call_op.se_idx(gui_call_op.curr_ch,1) = round(xsel/1e3*gui_call_op.fs);  % convert time to index

% update figure
axes(handles.axes_spectrogram);
set(gui_call_op.sline_spec,'XData',[1 1]*gui_call_op.se_idx(gui_call_op.curr_ch,1)/gui_call_op.fs*1e3);
axes(handles.axes_time_series);
set(gui_call_op.sline_time,'XData',[1 1]*gui_call_op.se_idx(gui_call_op.curr_ch,1)/gui_call_op.fs*1e3);

% Update call_dB data
data.proc.call_align_short_se_idx(gui_op.current_call_idx,gui_call_op.curr_ch,1) = gui_call_op.se_idx(gui_call_op.curr_ch,1);
data = update_call_time_series_data(data,gui_op,gui_call_op);
data = update_call_dB_data(data,gui_op,gui_call_op);

% Update main GUI figure
reload_bp_timeseries(hh_main,gui_op.mic_config);

% Store data & upload main GUI figure
setappdata(0,'gui_call_op',gui_call_op);
setappdata(0,'data',data);
reload_bp_timeseries(hh_main,gui_op.mic_config);



% --- Executes on button press in button_end.
function button_end_Callback(hObject, eventdata, handles)
% hObject    handle to button_end (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
data = getappdata(0,'data');
gui_call_op = getappdata(0,'gui_call_op');
gui_op = getappdata(0,'gui_op');
hh_main = getappdata(0,'bp_check_gui_handles');

% move end line
[xsel,~] = ginput(1);   % time
gui_call_op.se_idx(gui_call_op.curr_ch,2) = round(xsel/1e3*gui_call_op.fs);  % convert time to index

% update figure
axes(handles.axes_spectrogram);
set(gui_call_op.eline_spec,'XData',[1 1]*gui_call_op.se_idx(gui_call_op.curr_ch,2)/gui_call_op.fs*1e3);
axes(handles.axes_time_series);
set(gui_call_op.eline_time,'XData',[1 1]*gui_call_op.se_idx(gui_call_op.curr_ch,2)/gui_call_op.fs*1e3);

% Update call_dB data
data.proc.call_align_short_se_idx(gui_op.current_call_idx,gui_call_op.curr_ch,2) = gui_call_op.se_idx(gui_call_op.curr_ch,2);
data = update_call_time_series_data(data,gui_op,gui_call_op);
data = update_call_dB_data(data,gui_op,gui_call_op);

% Store data & upload main GUI figure
setappdata(0,'gui_call_op',gui_call_op);
setappdata(0,'data',data);
reload_bp_timeseries(hh_main,gui_op.mic_config);



% --- Executes on button press in button_delete_ch.
function button_delete_ch_Callback(hObject, eventdata, handles)
% hObject    handle to button_delete_ch (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
data = getappdata(0,'data');
gui_call_op = getappdata(0,'gui_call_op');
gui_op = getappdata(0,'gui_op');
hh_main = getappdata(0,'bp_check_gui_handles');

% Update data
gui_call_op.se_idx(gui_call_op.curr_ch,:) = nan(1,2);  % convert time to index
data.proc.call_align_short_se_idx(gui_op.current_call_idx,gui_call_op.curr_ch,:) = nan(1,2);

% update figure
axes(handles.axes_spectrogram);
set(gui_call_op.sline_spec,'XData',[1 1]*gui_call_op.se_idx(gui_call_op.curr_ch,1)/gui_call_op.fs*1e3);
set(gui_call_op.eline_spec,'XData',[1 1]*gui_call_op.se_idx(gui_call_op.curr_ch,2)/gui_call_op.fs*1e3);
axes(handles.axes_time_series);
set(gui_call_op.sline_time,'XData',[1 1]*gui_call_op.se_idx(gui_call_op.curr_ch,1)/gui_call_op.fs*1e3);
set(gui_call_op.eline_time,'XData',[1 1]*gui_call_op.se_idx(gui_call_op.curr_ch,2)/gui_call_op.fs*1e3);

% Update call_dB data
data = update_call_time_series_data(data,gui_op,gui_call_op);
data = update_call_dB_data(data,gui_op,gui_call_op);

% Store data & upload main GUI figure
setappdata(0,'gui_call_op',gui_call_op);
setappdata(0,'data',data);
reload_bp_timeseries(hh_main,gui_op.mic_config);


% --- Executes on button press in button_both.
function button_both_Callback(hObject, eventdata, handles)
% hObject    handle to button_both (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
data = getappdata(0,'data');
gui_call_op = getappdata(0,'gui_call_op');
gui_op = getappdata(0,'gui_op');
hh_main = getappdata(0,'bp_check_gui_handles');

% move end line
[xsel,~] = ginput(2);   % time
gui_call_op.se_idx(gui_call_op.curr_ch,:) = round(xsel/1e3*gui_call_op.fs);  % convert time to index

% update figure
axes(handles.axes_spectrogram);
set(gui_call_op.sline_spec,'XData',[1 1]*gui_call_op.se_idx(gui_call_op.curr_ch,1)/gui_call_op.fs*1e3);
set(gui_call_op.eline_spec,'XData',[1 1]*gui_call_op.se_idx(gui_call_op.curr_ch,2)/gui_call_op.fs*1e3);
axes(handles.axes_time_series);
set(gui_call_op.sline_time,'XData',[1 1]*gui_call_op.se_idx(gui_call_op.curr_ch,1)/gui_call_op.fs*1e3);
set(gui_call_op.eline_time,'XData',[1 1]*gui_call_op.se_idx(gui_call_op.curr_ch,2)/gui_call_op.fs*1e3);

% Update call_dB data
data.proc.call_align_short_se_idx(gui_op.current_call_idx,gui_call_op.curr_ch,:) = gui_call_op.se_idx(gui_call_op.curr_ch,:);
data = update_call_time_series_data(data,gui_op,gui_call_op);
data = update_call_dB_data(data,gui_op,gui_call_op);

% Store data & upload main GUI figure
setappdata(0,'gui_call_op',gui_call_op);
setappdata(0,'data',data);
reload_bp_timeseries(hh_main,gui_op.mic_config);



function reload_bp_timeseries(hh,mic_config)
% Reload bp plot and recording time series in main GUI
plot_time_series_in_gui(hh);  % display time series of first call
if strcmp(mic_config,'rb_cross');
    plot_bp_cross(hh);  % display beampattern
else
    plot_bp_2d(hh);  % display beampattern
end



function data = update_call_time_series_data(data,gui_op,gui_call_op)
% Update call_align_short time series
tukeywin_prop = data.param.tukeywin_proportion;  % tukey window taper porportion
iC = gui_op.current_call_idx;
iM = gui_call_op.curr_ch;
se_idx = gui_call_op.se_idx(iM,:);  % new segment start/end index
if any(isnan(se_idx))
    call_short = NaN;
else
    call_short = gui_call_op.call_align(se_idx(1):se_idx(2),gui_call_op.curr_ch)';  % extract new short call segment
end
w = tukeywin(length(call_short),tukeywin_prop);
call_short_taper = call_short.*w';  % taper call for spectrum estimation


call_freq_vec = linspace(0,data.mic_data.fs/2,round((size(call_short_taper,1)+1)/2));
if ~isfield(data.param,'PSD_type') || strcmp(data.param.PSD_type,'FFT')
    % Calculate fft ================================
    call_fft = fft(call_short_taper);
    call_fft_len = length(call_freq_vec);
    call_psd = 2*abs(call_fft(1:call_fft_len,:)).^2/(length(call_fft)*data.mic_data.fs);
    % call_psd = 2*abs(call_fft(1:call_fft_len,:)).^2;
    call_psd_dB = 10*log10(call_psd);
elseif strcmp(data.param.PSD_type,'pwelch')
    call_psd=pwelch(call_short_taper,128,120,call_freq_vec,data.mic_data.fs);
    call_psd_dB = 10*log10(call_psd);
end

%calculate RMS ==================================
call_rms = sqrt(mean(call_short_taper.^2));

data.proc.call_align_short{iC,iM} = call_short;
data.proc.call_align_short_se_idx(iC,iM,:) = se_idx;  % update start/end idx for this channel
data.proc.call_fft{iC,iM} = call_fft;  % call spectrum
data.proc.call_freq_vec{iC,iM} = call_freq_vec;  % frequency vector for call spectrum
data.proc.call_psd_raw_linear{iC,iM} = call_psd;  % spectrum of extracted calls, linear scale
data.proc.call_psd_raw_dB{iC,iM} = call_psd_dB;   % spectrum of extracted calls, dB scale
data.proc.call_rms(iC,iM) = call_rms;



function data = update_call_dB_data(data,gui_op,gui_call_op)
% Update call_dB data
iC = gui_op.current_call_idx;
iM = gui_call_op.curr_ch;
mic_to_bat_dist = data.proc.mic_to_bat_dist(iC,iM);     % distance between bat and mic
bat_to_mic_angle = data.proc.bat_to_mic_angle(iC,iM);   % angle to compensate for mic beampattern
call_psd_raw_dB = data.proc.call_psd_raw_dB{iC,iM};  % call spectrum in dB scale
call_freq = data.proc.call_freq_vec{iC,iM};  % frequency vector of the call spectrum

% Transmission loss: air absorption and spreading loss
[alpha, alpha_iso,~,~] = ...
    air_absorption_vec(call_freq,data.param.tempC,data.param.humid);  % atmospheric aborption [dB/m]

% Attenuation and spreading loss
air_attn_dB = alpha*(mic_to_bat_dist'-data.param.d0);  % [dB]
spreading_loss_dB = 20*log10(mic_to_bat_dist'/data.param.d0);  % [dB]
TL_dB = air_attn_dB+spreading_loss_dB';  % [dB]

% Mic sensitivity
mic_sens_dB = interp1(data.mic_sens.freq,data.mic_sens.sens(:,iM),call_freq);  % interpolate mic sensitivity vector

% Mic beampattern: have to loop because interp2 is needed but bp is in 3D mtx
[X,Y] = meshgrid(data.mic_bp.theta,data.mic_bp.freq);
[XI,YI] = meshgrid(bat_to_mic_angle,call_freq);
bp_compensation = interp2(X,Y,data.mic_bp.bp(:,:,iM),XI,YI)';  % interpolate mic beampattern

% Broadband call SPL and PSD
call_psd_dB_comp_nobp = call_psd_raw_dB + TL_dB - mic_sens_dB;
call_psd_dB_comp_withbp = call_psd_raw_dB + TL_dB - mic_sens_dB - bp_compensation;
gain_dB_fac = data.mic_gain(iM);
call_psd_dB_comp_re20uPa_nobp = call_psd_dB_comp_nobp + 20*log10(1/20e-6) - gain_dB_fac;
call_psd_dB_comp_re20uPa_withbp = call_psd_dB_comp_withbp + 20*log10(1/20e-6) - gain_dB_fac;

% Call SPL using p2p voltage
ch_data = cell2mat(data.proc.call_align_short(iC,:)');
ch_p2p = max(ch_data,[],2) + ...
  abs(min(cell2mat(data.proc.call_align_short(iC,:)'),[],2));
call_p2p_ch_dB = 20*log10(ch_p2p');

%using pwelch to calc. the max freq in the data
A=pwelch(ch_data',128,120,call_freq,data.mic_data.fs);
[~,max_freq]=max(A); %where the peak should be...
max_freq_idx = sub2ind(size(A), max_freq, 1:size(A,2));

TL_dB_ch = TL_dB(max_freq_idx);
mic_sens_dB_mean_ch = mic_sens_dB(max_freq_idx);
call_p2p_SPL_comp_re20uPa = call_p2p_ch_dB + TL_dB_ch - mic_sens_dB_mean_ch +...
  20*log10(1/20e-6) - data.mic_gain';

% Save data
data.param.alpha{iC,iM} = alpha;
data.param.alpha_iso{iC,iM} = alpha_iso;
data.proc.air_attn_dB{iC,iM} = air_attn_dB;
data.proc.spreading_loss_dB{iC,iM} = spreading_loss_dB;
data.proc.TL_dB{iC,iM} = TL_dB;
data.proc.mic_bp_compensation_dB{iC,iM} = bp_compensation;
data.proc.mic_sens_dB{iC,iM} = mic_sens_dB;
data.proc.call_psd_dB_comp_nobp{iC,iM} = call_psd_dB_comp_nobp;
data.proc.call_psd_dB_comp_withbp{iC,iM} = call_psd_dB_comp_withbp;
data.proc.call_psd_dB_comp_re20uPa_nobp{iC,iM} = call_psd_dB_comp_re20uPa_nobp;
data.proc.call_psd_dB_comp_re20uPa_withbp{iC,iM} = call_psd_dB_comp_re20uPa_withbp;
data.proc.call_p2p_SPL_comp_re20uPa(iC,:) = call_p2p_SPL_comp_re20uPa;




% --- Executes on button press in button_ch_next.
function button_ch_next_Callback(hObject, eventdata, handles)
% hObject    handle to button_ch_next (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
gui_call_op = getappdata(0,'gui_call_op');
gui_op = getappdata(0,'gui_op');
gui_call_op.curr_ch = min(gui_call_op.curr_ch+1,gui_call_op.num_ch_in_file);
set(handles.edit_ch,'String',num2str(gui_call_op.curr_ch));

setappdata(0,'gui_call_op',gui_call_op);

plot_everything(handles,gui_op.current_call_idx);  % update spectrogram and time series


% --- Executes on button press in butto_ch_previous.
function butto_ch_previous_Callback(hObject, eventdata, handles)
% hObject    handle to butto_ch_previous (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
gui_call_op = getappdata(0,'gui_call_op');
gui_op = getappdata(0,'gui_op');
gui_call_op.curr_ch = max(gui_call_op.curr_ch-1,1);
set(handles.edit_ch,'String',num2str(gui_call_op.curr_ch));

setappdata(0,'gui_call_op',gui_call_op);

plot_everything(handles,gui_op.current_call_idx);  % update spectrogram and time series


function edit_ch_Callback(hObject, eventdata, handles)
% hObject    handle to edit_ch (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
gui_call_op = getappdata(0,'gui_call_op');
gui_op = getappdata(0,'gui_op');
tmp = mod(str2double(get(hObject,'String')),gui_call_op.num_ch_in_file);
if tmp==0
    gui_call_op.curr_ch = gui_call_op.num_ch_in_file;
else
    gui_call_op.curr_ch = tmp;
end
setappdata(0,'gui_call_op',gui_call_op);

plot_everything(handles,gui_op.current_call_idx);  % update spectrogram and time series



function plot_everything(handles,curr_call_idx)
gui_call_op = getappdata(0,'gui_call_op');
[~,F,T,P] = spectrogram(gui_call_op.call_align(:,gui_call_op.curr_ch),128,120,128,gui_call_op.fs);
P = 10*log10(abs(P));
xlim_curr = T([1 end])*1e3;

axes(handles.axes_spectrogram);
imagesc(T*1e3,F/1e3,P);
axis xy
hold on
gui_call_op.sline_spec = plot([1 1]*gui_call_op.se_idx(gui_call_op.curr_ch,1)/gui_call_op.fs*1e3,[0 gui_call_op.fs/1e3/2],'r','linewidth',2);
gui_call_op.eline_spec = plot([1 1]*gui_call_op.se_idx(gui_call_op.curr_ch,2)/gui_call_op.fs*1e3,[0 gui_call_op.fs/1e3/2],'r','linewidth',2);
hold off
xlim(xlim_curr);
set(gca,'xticklabel','');
ylabel('Frequency (kHz)');
title(sprintf('Call #%d',curr_call_idx),'fontsize',16);

axes(handles.axes_time_series);
cla
plot(gui_call_op.time*1e3,gui_call_op.call_align(:,gui_call_op.curr_ch));
yybnd = ceil(max(gui_call_op.call_align(:,gui_call_op.curr_ch))/0.01)*0.01*[-1 1];
hold on
gui_call_op.sline_time = plot([1 1]*gui_call_op.se_idx(gui_call_op.curr_ch,1)/gui_call_op.fs*1e3,[-5 5],'r','linewidth',2);
gui_call_op.eline_time = plot([1 1]*gui_call_op.se_idx(gui_call_op.curr_ch,2)/gui_call_op.fs*1e3,[-5 5],'r','linewidth',2);
hold off
xlim(xlim_curr);
ylim(yybnd);
xlabel('Time (ms)');
ylabel('Volt');

setappdata(0,'gui_call_op',gui_call_op);


% --- Executes when user attempts to close figure1.
function figure1_CloseRequestFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if isappdata(0,'gui_call_op')
    rmappdata(0,'gui_call_op');
end
if isappdata(0,'edit_call_gui_handles')
    rmappdata(0,'edit_call_gui_handles');
end

% Hint: delete(hObject) closes the figure
delete(hObject);
