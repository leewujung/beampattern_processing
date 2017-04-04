function varargout = bp_check_gui(varargin)
% BP_CHECK_GUI MATLAB code for bp_check_gui.fig
%      BP_CHECK_GUI, by itself, creates a new BP_CHECK_GUI or raises the existing
%      singleton*.
%
%      H = BP_CHECK_GUI returns the handle to a new BP_CHECK_GUI or the handle to
%      the existing singleton*.
%
%      BP_CHECK_GUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in BP_CHECK_GUI.M with the given input arguments.
%
%      BP_CHECK_GUI('Property','Value',...) creates a new BP_CHECK_GUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before bp_check_gui_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to bp_check_gui_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help bp_check_gui

% Last Modified by GUIDE v2.5 31-Mar-2017 14:15:42

% 2015 10 13  -- feed bat head aim from data
%             -- use new format of mic sensitivity and beampattern
% 2015 10 19  -- use all 34 channels with corresponding sensitivity and
%                beampattern
%             -- make mic gain a separate sheet to be loaded for flexiblity
% 2015 10 20  -- get rid of mic_seq in the code, just order mic_loc according to
%                ch number
%             -- track index of extracted sections
%             -- kick out bad channel
%             -- checkbox for bad call/click
% 2015 10 21  -- plot for linear configuration  --> beampattern_gui_v6.m
% ------------------------------------------------------------------------
% 2015 10 21  -- move all calculation outside of GUI and use GUI only for
%                visualization and checking each call
% 2015 10 22  -- fixed everything to load processed data
% Things to fix: -- need to load mic_data
%                -- no call_align/call_no_align and need to plot time trace
%                -- need to initialize mic_config part
% TODO: -- enable changing call duration in GUI
%       -- make save button



% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @bp_check_gui_OpeningFcn, ...
                   'gui_OutputFcn',  @bp_check_gui_OutputFcn, ...
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


% --- Executes just before bp_check_gui is made visible.
function bp_check_gui_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to bp_check_gui (see VARARGIN)

% Choose default command line output for bp_check_gui
handles.output = hObject;

% Enable figure toolbar
set(hObject,'toolbar','figure');

% Update handles structure
guidata(hObject, handles);

% Set main GUI handle to appdata
setappdata(0,'bp_check_gui_handles',handles);

% UIWAIT makes bp_check_gui wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = bp_check_gui_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;



% --- Executes on button press in button_proc_file.
function button_proc_file_Callback(hObject, eventdata, handles)
% hObject    handle to button_proc_file (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
gui_op = getappdata(0,'gui_op');

% Load processed file
if ispref('bp_check_gui') && ispref('bp_check_gui','proc_file_path')
    pname = getpref('bp_check_gui','proc_file_path');
else
    pname = '';
end
[fname,pname] = uigetfile('*.mat','Pick processed file',pname);
if isequal(pname,0)
    return
end
setpref('bp_check_gui','proc_file_path',pname)
disp(['Loading processed file: ',fname]);
data = load(fullfile(pname,fname));  % load all processed data
set(handles.edit_proc_file,'String',fname);

gui_op.base_dir = fileparts(fileparts(pname));  % base_dir is one level up from the processed data folder
if ~exist(gui_op.base_dir,'dir')  % if base_dir doesn't exist
    gui_op.base_dir = pwd;
end
% cd(gui_op.base_dir)

gui_op.current_call_idx = 1;
data.path.proc_data = pname;
data.files.proc_data = fname;
disp('Processed results loaded');

% Load mic data
disp('Loading raw mic signals...')
if ~ispref('bp_check_gui','mic_data_path')
    md_path = fullfile(gui_op.base_dir,data.path.mic_data);
else
    md_path = getpref('bp_check_gui','mic_data_path');
end
if ~exist(fullfile(md_path,[data.files.mic_data,'.mat']),'file')  % if cannot find mic_data in path then select new one
    md_path = uigetdir(md_path,'Select the directory with mic_data files');
    if isequal(md_path,0)
        return;
    end
    if ~exist(fullfile(md_path,[data.files.mic_data,'.mat']),'file')
        disp('Specified mic data file not found in the path, try again');
        fprintf('File wanted: %s',[data.files.mic_data,'.mat']);
        return;
    end
end
setpref('bp_check_gui','mic_data_path',md_path);
gui_op.path_mic_data = md_path;
A = load(fullfile(gui_op.path_mic_data,data.files.mic_data));
data.mic_data.sig = A.sig;  % save mic signals into the current data structure


% Get current mic config display setting
gui_op.mic_config = handles.config_radio_grp.SelectedObject.Tag;  % configuration selection
gui_op.linlog = handles.loglin_radio_grp.SelectedObject.Tag;  % linear/log selection
gui_op.interp = handles.interp_radio_grp.SelectedObject.Tag;  % interpolation selection

% Get default caxis info
freq_wanted = str2double(get(handles.edit_bp_freq,'String'))*1e3;  % beampattern frequency [Hz]
call_dB = nan(data.mic_data.num_ch_in_file,1);
for iM=1:data.mic_data.num_ch_in_file
    freq = data.proc.call_freq_vec{gui_op.current_call_idx,iM};
    [~,fidx] = min(abs(freq-freq_wanted));
    call_dB(iM) = data.proc.call_psd_dB_comp_re20uPa_withbp{gui_op.current_call_idx,iM}(fidx);
end
gui_op.caxis_raw_default = [floor(min(min(call_dB))/5)*5, ceil(max(max(call_dB))/5)*5];
gui_op.caxis_norm_default = [floor(-range(call_dB)/5)*5, 0];
gui_op.caxis_raw_current = gui_op.caxis_raw_default;
gui_op.caxis_norm_current = gui_op.caxis_norm_default;

% Change folder and save stuff
setappdata(0,'data',data);
setappdata(0,'gui_op',gui_op);

% Open bat track window
view_track_gui;

% Update figure and display info for call#1
set(handles.text_current_call1,'String',num2str(gui_op.current_call_idx));
set(handles.text_current_call2,'String',['/',num2str(length(data.mic_data.call_idx_w_track))]);
update_good_call(handles);  % udpate good call checkbox
update_head_aim_mkr(handles);  % update head aim source
update_ch_ex(handles);     % update list of channel to be excluded
update_caxis(handles);     % update color axis for bp display
update_edit_call_gui;  % update edit_call_section GUI if exist

plot_bat_mic_vector;   % plot bat2mic vector
plot_time_series_in_gui(handles);  % display time series of first call
update_bp_plots(handles,gui_op);
disp('Loaded ')


function update_caxis(handles)
gui_op = getappdata(0,'gui_op');
if get(handles.checkbox_norm,'Value')==0  % if "normalized" not checked
    set(handles.edit_cmin,'String',num2str(gui_op.caxis_raw_current(1)));
    set(handles.edit_cmax,'String',num2str(gui_op.caxis_raw_current(2)));
else
    set(handles.edit_cmin,'String',num2str(gui_op.caxis_norm_current(1)));
    set(handles.edit_cmax,'String',num2str(gui_op.caxis_norm_current(2)));
end


% --- Executes on button press in button_previous_call.
function button_previous_call_Callback(hObject, eventdata, handles)
% hObject    handle to button_previous_call (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
data = getappdata(0,'data');
gui_op = getappdata(0,'gui_op');
if gui_op.current_call_idx~=1
    gui_op.current_call_idx = gui_op.current_call_idx-1;
end
setappdata(0,'data',data);
setappdata(0,'gui_op',gui_op);

% Update info
update_good_call(handles);
update_head_aim_mkr(handles);  % update head aim source
update_ch_ex(handles);
update_call_num(handles);  % updated displayed call number
update_caxis(handles);     % update color axis for bp display
move_call_circle_on_track;  % update current call location on track
update_edit_call_gui;  % update edit_call_section GUI if exist

% Update plot
plot_bat_mic_vector;   % plot bat2mic vector
plot_time_series_in_gui(handles);  % display time series of first call
update_bp_plots(handles,gui_op);


function update_bp_plots(handles,gui_op)
if strcmp(gui_op.mic_config,'rb_cross')
  plot_bp_cross(handles);  % display beampattern
elseif strcmp(gui_op.mic_config,'rb_multi_freq')
  plot_bp_multi_freq(handles);
else
  plot_bp_2d(handles);  % display beampattern
end


% --- Executes on button press in button_next_call.
function button_next_call_Callback(hObject, eventdata, handles)
% hObject    handle to button_next_call (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
data = getappdata(0,'data');
gui_op = getappdata(0,'gui_op');
if gui_op.current_call_idx~=length(data.mic_data.call_idx_w_track)
    gui_op.current_call_idx = gui_op.current_call_idx+1;
end
setappdata(0,'data',data);
setappdata(0,'gui_op',gui_op);

% Update info
update_good_call(handles);
update_head_aim_mkr(handles);  % update head aim source
update_ch_ex(handles);
update_call_num(handles);  % updated displayed call number
update_caxis(handles);     % update color axis for bp display
move_call_circle_on_track();  % update current call location on track
update_edit_call_gui;  % update edit_call_section GUI if exist

% Update plot
plot_bat_mic_vector;   % plot bat2mic vector
plot_time_series_in_gui(handles);  % display time series of first call
update_bp_plots(handles,gui_op);


function update_edit_call_gui()
if isappdata(0,'edit_call_gui_handles')
    edit_call_section;
end


function text_current_call1_Callback(hObject, eventdata, handles)
% hObject    handle to text_current_call1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
data = getappdata(0,'data');
gui_op = getappdata(0,'gui_op');
entered_call_num=str2double(get(hObject,'String'));
if entered_call_num > length(data.mic_data.call_idx_w_track)
    entered_call_num=length(data.mic_data.call_idx_w_track);
elseif entered_call_num < 1
    entered_call_num=1;
end
gui_op.current_call_idx = entered_call_num;
setappdata(0,'data',data);
setappdata(0,'gui_op',gui_op);

% Update info
update_good_call(handles);
update_head_aim_mkr(handles);  % update head aim source
update_ch_ex(handles);
update_call_num(handles);  % updated displayed call number
update_caxis(handles);     % update color axis for bp display
move_call_circle_on_track();  % update current call location on track

% Update plot
plot_bat_mic_vector;   % plot bat2mic vector
plot_time_series_in_gui(handles);  % display time series of first call
update_bp_plots(handles,gui_op);



function update_ch_ex(handles)
data = getappdata(0,'data');
gui_op = getappdata(0,'gui_op');
vv = data.proc.ch_ex{gui_op.current_call_idx};

if ~isempty(vv)
    ss = [];
    for iV=1:length(vv)
        ss = [ss,', ',num2str(vv(iV))];
    end
    ss = ss(3:end);
    set(handles.edit_ch_ex,'String',ss);
else
    set(handles.edit_ch_ex,'String','');
end


function update_call_num(handles)
gui_op = getappdata(0,'gui_op');
set(handles.text_current_call1,'String',num2str(gui_op.current_call_idx));
% Update call number in the track and mic GUI
if isappdata(0,'track_gui_handles');
    track_gui_handles = getappdata(0,'track_gui_handles');
    set(track_gui_handles.edit_call_idx1,'String',num2str(gui_op.current_call_idx));
end


function move_call_circle_on_track()
data = getappdata(0,'data');
gui_op = getappdata(0,'gui_op');
if isappdata(0,'track_gui_handles');
    track_gui_handles = getappdata(0,'track_gui_handles');
    axes(track_gui_handles.axes2);
    set(track_gui_handles.h_curr_call,...
        'Xdata',data.track.track_interp(data.track.call_loc_idx_on_track_interp(gui_op.current_call_idx),1),...
        'Ydata',data.track.track_interp(data.track.call_loc_idx_on_track_interp(gui_op.current_call_idx),2),...
        'Zdata',data.track.track_interp(data.track.call_loc_idx_on_track_interp(gui_op.current_call_idx),3));
end



% --- Executes on button press in checkbox_norm.
function checkbox_norm_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_norm (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
gui_op = getappdata(0,'gui_op');
data = getappdata(0,'data');
update_caxis(handles);     % update color axis for bp display
if isfield(data.proc,'call_psd_dB_comp_re20uPa_withbp')  % if data already loaded
  update_bp_plots(handles,gui_op);
end
% Hint: get(hObject,'Value') returns toggle state of checkbox_norm



% --- Executes when user attempts to close figure1.
function figure1_CloseRequestFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% close track figure
hh2 = getappdata(0,'track_gui_handles');
if isfield(hh2,'figure1')
    delete(hh2.figure1);
end
hh3 = getappdata(0,'edit_call_gui_handles');
if isfield(hh3,'figure1')
    delete(hh3.figure1);
end

% delete appdata
if isappdata(0,'data')
    rmappdata(0,'data');
end
if isappdata(0,'gui_op')
    rmappdata(0,'gui_op');
end
if isappdata(0,'gui_call_op')
    rmappdata(0,'gui_call_op');
end
if isappdata(0,'track_gui_handles')
    rmappdata(0,'track_gui_handles');
end
if isappdata(0,'edit_call_gui_handles')
    rmappdata(0,'edit_call_gui_handles');
end
if isappdata(0,'bp_check_gui_handles')
    rmappdata(0,'bp_check_gui_handles');
end

% Hint: delete(hObject) closes the figure
delete(hObject);



% --- Executes on button press in button_save.
function button_save_Callback(hObject, eventdata, handles)
% hObject    handle to button_save (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
data = getappdata(0,'data');
[save_fname,save_pname] = uiputfile('*.mat','Save detection results',...
  fullfile(data.path.proc_data,data.files.proc_data));
if isfield(data.mic_data,'sig')
    data.mic_data.sig = [];
end    
if isequal(save_fname,0)
    return
end
disp('Saving...')
save(fullfile(save_pname,save_fname),'-struct','data');
disp('Saved')




% --- Executes on button press in button_export_good_calls.
function button_export_good_calls_Callback(hObject, eventdata, handles)
% hObject    handle to button_export_good_calls (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
data = getappdata(0,'data');
[save_fname,save_pname] = uiputfile('*.mat','Export checked results',...
    fullfile(data.path.proc_data,[data.files.proc_data(1:end-4) '_checked.mat']));
if isequal(save_fname,0)
    return
end

export_data=struct();
export_data.proc.chk_good_call = data.proc.chk_good_call;
export_data.proc.ch_ex = data.proc.ch_ex;

save(fullfile(save_pname,save_fname),'-struct','export_data');
disp('Exported')

function edit_bp_freq_Callback(hObject, eventdata, handles)
% hObject    handle to edit_bp_freq (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
gui_op = getappdata(0,'gui_op');
data = getappdata(0,'data');
update_caxis(handles);     % update color axis for bp display
if isfield(data.proc,'call_psd_dB_comp_re20uPa_withbp')  % if data already loaded
  update_bp_plots(handles,gui_op);
end
if isappdata(0,'track_gui_handles')
  track_gui_handles = getappdata(0,'track_gui_handles'); 
  if get(track_gui_handles.checkbox_top_view_bp,'value')
    plot_bat_mic_vector;
  end
end

% Hints: get(hObject,'String') returns contents of edit_bp_freq as text
%        str2double(get(hObject,'String')) returns contents of edit_bp_freq as a double


% --- Executes during object creation, after setting all properties.
function edit_bp_freq_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_bp_freq (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in checkbox_good_call.
function checkbox_good_call_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_good_call (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
data = getappdata(0,'data');
gui_op = getappdata(0,'gui_op');
data.proc.chk_good_call(gui_op.current_call_idx) = get(hObject,'Value');
setappdata(0,'data',data);
setappdata(0,'gui_op',gui_op);


function update_good_call(handles)
data = getappdata(0,'data');
gui_op = getappdata(0,'gui_op');
if data.proc.chk_good_call(gui_op.current_call_idx)==1
    set(handles.checkbox_good_call,'Value',1);
else
    set(handles.checkbox_good_call,'Value',0);
end



function edit_ch_ex_Callback(hObject, eventdata, handles)
% hObject    handle to edit_ch_ex (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
s = get(hObject,'String');

data = getappdata(0,'data');
gui_op = getappdata(0,'gui_op');
data.proc.ch_ex{gui_op.current_call_idx} = sscanf(s,'%d,');
setappdata(0,'data',data);

gui_op = getappdata(0,'gui_op');
update_bp_plots(handles,gui_op);
if isappdata(0,'track_gui_handles')
  track_gui_handles = getappdata(0,'track_gui_handles'); 
  if get(track_gui_handles.checkbox_top_view_bp,'value')
    plot_bat_mic_vector;
  end
end

% Hints: get(hObject,'String') returns contents of edit_ch_ex as text
%        str2double(get(hObject,'String')) returns contents of edit_ch_ex as a double


% --- Executes during object creation, after setting all properties.
function edit_ch_ex_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_ch_ex (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes when selected object is changed in config_radio_grp.
function config_radio_grp_SelectionChangedFcn(hObject, eventdata, handles)
% hObject    handle to the selected object in config_radio_grp 
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
gui_op = getappdata(0,'gui_op');
gui_op.mic_config = get(eventdata.NewValue, 'Tag');
setappdata(0,'gui_op',gui_op);

cla(handles.axes_bp_contour,'reset');
cla(handles.axes_bp,'reset');

data = getappdata(0,'data');
if isfield(data.proc,'call_psd_dB_comp_re20uPa_withbp')
  update_bp_plots(handles,gui_op);
end
if ~get(handles.show_bp_plot,'Value') && ~strcmp(gui_op.mic_config,'rb_cross')
  set(handles.axes_bp,'visible','off')
end




function edit_proc_file_Callback(hObject, eventdata, handles)
% hObject    handle to edit_proc_file (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_proc_file as text
%        str2double(get(hObject,'String')) returns contents of edit_proc_file as a double


% --- Executes during object creation, after setting all properties.
function edit_proc_file_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_proc_file (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in rb_log.
function rb_log_Callback(hObject, eventdata, handles)
% hObject    handle to rb_log (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of rb_log


% --- Executes when selected object is changed in loglin_radio_grp.
function loglin_radio_grp_SelectionChangedFcn(hObject, eventdata, handles)
% hObject    handle to the selected object in loglin_radio_grp 
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
gui_op = getappdata(0,'gui_op');
gui_op.linlog = get(eventdata.NewValue, 'Tag');
setappdata(0,'gui_op',gui_op);

data = getappdata(0,'data');
if isfield(data.proc,'call_psd_dB_comp_re20uPa_withbp')  % if data already loaded
  update_bp_plots(handles,gui_op);
end
if isappdata(0,'track_gui_handles')
  track_gui_handles = getappdata(0,'track_gui_handles'); 
  if get(track_gui_handles.checkbox_top_view_bp,'value')
    plot_bat_mic_vector;
  end
end


% --- Executes when selected object is changed in interp_radio_grp.
function interp_radio_grp_SelectionChangedFcn(hObject, eventdata, handles)
% hObject    handle to the selected object in interp_radio_grp 
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
gui_op = getappdata(0,'gui_op');
gui_op.interp = get(eventdata.NewValue, 'Tag');
setappdata(0,'gui_op',gui_op);

data = getappdata(0,'data');
if isfield(data.proc,'call_psd_dB_comp_re20uPa_withbp')  % if data already loaded
  update_bp_plots(handles,gui_op);
end



function edit_cmin_Callback(hObject, eventdata, handles)
% hObject    handle to edit_cmin (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
gui_op = getappdata(0,'gui_op');
if get(handles.checkbox_norm,'Value')==0  % if "normalized" not checked
    gui_op.caxis_raw_current(1) = str2double(get(hObject,'String'));
else
    gui_op.caxis_norm_current(1) = str2double(get(hObject,'String'));
end
setappdata(0,'gui_op',gui_op);
update_bp_plots(handles,gui_op);
% Hints: get(hObject,'String') returns contents of edit_cmin as text
%        str2double(get(hObject,'String')) returns contents of edit_cmin as a double


% --- Executes during object creation, after setting all properties.
function edit_cmin_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_cmin (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function edit_cmax_Callback(hObject, eventdata, handles)
% hObject    handle to edit_cmax (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
gui_op = getappdata(0,'gui_op');
if get(handles.checkbox_norm,'Value')==0  % if "normalized" not checked
    gui_op.caxis_raw_current(2) = str2double(get(hObject,'String'));
else
    gui_op.caxis_norm_current(2) = str2double(get(hObject,'String'));
end
setappdata(0,'gui_op',gui_op);
update_bp_plots(handles,gui_op);
% Hints: get(hObject,'String') returns contents of edit_cmax as text
%        str2double(get(hObject,'String')) returns contents of edit_cmax as a double


% --- Executes during object creation, after setting all properties.
function edit_cmax_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_cmax (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_peak_db_Callback(hObject, eventdata, handles)
% hObject    handle to edit_peak_db (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_peak_db as text
%        str2double(get(hObject,'String')) returns contents of edit_peak_db as a double


% --- Executes during object creation, after setting all properties.
function edit_peak_db_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_peak_db (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on key press with focus on figure1 and none of its controls.
function figure1_KeyPressFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  structure with the following fields (see MATLAB.UI.FIGURE)
%	Key: name of the key that was pressed, in lower case
%	Character: character interpretation of the key(s) that was pressed
%	Modifier: name(s) of the modifier key(s) (i.e., control, shift) pressed
% handles    structure with handles and user data (see GUIDATA)
switch eventdata.Key
    case {'numpad4','leftarrow'}
        button_previous_call_Callback(handles.button_previous_call, eventdata, handles);
    case {'numpad6','rightarrow'}
        button_next_call_Callback(handles.button_next_call, eventdata, handles);
end



% --- Executes on button press in button_edit_call_sec.
function button_edit_call_sec_Callback(hObject, eventdata, handles)
% hObject    handle to button_edit_call_sec (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
edit_call_section;



% --- Executes on button press in checkbox_head_aim_mkr.
function checkbox_head_aim_mkr_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_head_aim_mkr (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_head_aim_mkr


function update_head_aim_mkr(handles)
data = getappdata(0,'data');
gui_op = getappdata(0,'gui_op');
if data.proc.source_head_aim(gui_op.current_call_idx)==1
    set(handles.checkbox_head_aim_mkr,'Value',1);
else
    set(handles.checkbox_head_aim_mkr,'Value',0);
end


% --- Executes on button press in bp_freq_anim.
function bp_freq_anim_Callback(hObject, eventdata, handles)
% hObject    handle to bp_freq_anim (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

data = getappdata(0,'data');
gui_op = getappdata(0,'gui_op');

%check for whether you are in multi_freq
if strcmp(gui_op.mic_config,'rb_multi_freq')
  disp('Can''t animate in multi freq mode.')
  return
end

start_freq = str2double(get(handles.bp_freq_anim_start,'String'))*1e3;
end_freq = str2double(get(handles.bp_freq_anim_end,'String'))*1e3;

current_call = gui_op.current_call_idx;
freq_vec = data.proc.call_freq_vec{current_call,1};

def_interv = diff(freq_vec(1:2));
new_interv = str2double(get(handles.anim_khz_size,'String'))*1e3;
skip_idx = round(new_interv/def_interv);
freq_vec = freq_vec(1:skip_idx:end);

freq_idx = find(freq_vec>=start_freq+new_interv,1):find(freq_vec>=end_freq,1);
for ff=1:length(freq_idx)
  freq_wanted = freq_vec(freq_idx(ff));
  if strcmp(gui_op.mic_config,'rb_cross')
    plot_bp_cross(handles,freq_wanted);  % display beampattern
  else
    plot_bp_2d(handles,freq_wanted);  % display beampattern
  end
  title(['Freq: ' num2str(round(freq_wanted/1e3))])
  drawnow;
end

%returning to original freq display:
freq_wanted = str2double(get(handles.edit_bp_freq,'String'))*1e3;
if strcmp(gui_op.mic_config,'rb_cross')
  plot_bp_cross(handles,freq_wanted);  % display beampattern
else
  plot_bp_2d(handles,freq_wanted);  % display beampattern
end

% fig_h = figure(3);
% for ff=1:length(freq_idx)
%   plot those data in a new window with (hopefully) minimal delay
% end


function bp_freq_anim_end_Callback(hObject, eventdata, handles)
% hObject    handle to bp_freq_anim_end (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of bp_freq_anim_end as text
%        str2double(get(hObject,'String')) returns contents of bp_freq_anim_end as a double


function bp_freq_anim_start_Callback(hObject, eventdata, handles)
% hObject    handle to bp_freq_anim_start (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of bp_freq_anim_start as text
%        str2double(get(hObject,'String')) returns contents of bp_freq_anim_start as a double


% --- Executes during object creation, after setting all properties.
function bp_freq_anim_start_CreateFcn(hObject, eventdata, handles)
% hObject    handle to bp_freq_anim_start (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes during object creation, after setting all properties.
function bp_freq_anim_end_CreateFcn(hObject, eventdata, handles)
% hObject    handle to bp_freq_anim_end (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function anim_khz_size_Callback(hObject, eventdata, handles)
% hObject    handle to anim_khz_size (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of anim_khz_size as text
%        str2double(get(hObject,'String')) returns contents of anim_khz_size as a double


% --- Executes during object creation, after setting all properties.
function anim_khz_size_CreateFcn(hObject, eventdata, handles)
% hObject    handle to anim_khz_size (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in show_bp_plot.
function show_bp_plot_Callback(hObject, eventdata, handles)
% hObject    handle to show_bp_plot (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of show_bp_plot
gui_op = getappdata(0,'gui_op');
if strcmp(gui_op.mic_config,'rb_2d')
  if get(hObject,'Value')
    plot_bp_2d(handles);
  else
    axes(handles.axes_bp);
    cla(handles.axes_bp,'reset');
    set(handles.axes_bp,'visible','off');
  end
end
