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

% Last Modified by GUIDE v2.5 16-Nov-2015 21:23:57

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

plot_everything(handles);

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

% store data
data.proc.call_align_short_se_idx(gui_op.current_call_idx,gui_call_op.curr_ch,1) = gui_call_op.se_idx(gui_call_op.curr_ch,1);
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

% store data
data.proc.call_align_short_se_idx(gui_op.current_call_idx,gui_call_op.curr_ch,2) = gui_call_op.se_idx(gui_call_op.curr_ch,2);
setappdata(0,'gui_call_op',gui_call_op);
setappdata(0,'data',data);
reload_bp_timeseries(hh_main,gui_op.mic_config);



function reload_bp_timeseries(hh,mic_config)
plot_time_series_in_gui(hh);  % display time series of first call
if strcmp(mic_config,'rb_cross');
    plot_bp_cross(hh);  % display beampattern
else
    plot_bp_2d(hh);  % display beampattern
end


% --- Executes on button press in button_ch_next.
function button_ch_next_Callback(hObject, eventdata, handles)
% hObject    handle to button_ch_next (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
gui_call_op = getappdata(0,'gui_call_op');

tmp = mod(gui_call_op.curr_ch+1,gui_call_op.num_ch_in_file);
if tmp==0
    gui_call_op.curr_ch = gui_call_op.num_ch_in_file;
else
    gui_call_op.curr_ch = tmp;
end
set(handles.edit_ch,'String',num2str(gui_call_op.curr_ch));

setappdata(0,'gui_call_op',gui_call_op);

plot_everything(handles);  % update spectrogram and time series


% --- Executes on button press in butto_ch_previous.
function butto_ch_previous_Callback(hObject, eventdata, handles)
% hObject    handle to butto_ch_previous (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
gui_call_op = getappdata(0,'gui_call_op');

tmp = mod(gui_call_op.curr_ch-1,gui_call_op.num_ch_in_file);
if tmp==0
    gui_call_op.curr_ch = gui_call_op.num_ch_in_file;
else
    gui_call_op.curr_ch = tmp;
end
set(handles.edit_ch,'String',num2str(gui_call_op.curr_ch));

setappdata(0,'gui_call_op',gui_call_op);

plot_everything(handles);  % update spectrogram and time series


function edit_ch_Callback(hObject, eventdata, handles)
% hObject    handle to edit_ch (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
gui_call_op = getappdata(0,'gui_call_op');
tmp = mod(str2double(get(hObject,'String')),gui_call_op.num_ch_in_file);
if tmp==0
    gui_call_op.curr_ch = gui_call_op.num_ch_in_file;
else
    gui_call_op.curr_ch = tmp;
end
setappdata(0,'gui_call_op',gui_call_op);

plot_everything(handles);  % update spectrogram and time series


function plot_everything(handles)
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
ylabel('Frequency (kHz)');

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
