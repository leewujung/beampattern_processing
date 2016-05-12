function varargout = view_track_gui(varargin)
% VIEW_TRACK_GUI MATLAB code for view_track_gui.fig
%      VIEW_TRACK_GUI, by itself, creates a new VIEW_TRACK_GUI or raises the existing
%      singleton*.
%
%      H = VIEW_TRACK_GUI returns the handle to a new VIEW_TRACK_GUI or the handle to
%      the existing singleton*.
%
%      VIEW_TRACK_GUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in VIEW_TRACK_GUI.M with the given input arguments.
%
%      VIEW_TRACK_GUI('Property','Value',...) creates a new VIEW_TRACK_GUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before view_track_gui_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to view_track_gui_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help view_track_gui

% Last Modified by GUIDE v2.5 11-May-2016 12:05:23

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @view_track_gui_OpeningFcn, ...
                   'gui_OutputFcn',  @view_track_gui_OutputFcn, ...
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


% --- Executes just before view_track_gui is made visible.
function view_track_gui_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to view_track_gui (see VARARGIN)

% Enable figure toolbar
set(hObject,'toolbar','figure');
cla reset

% Get main GUI data
data = getappdata(0,'data');
gui_op = getappdata(0,'gui_op');
curr_call_idx = gui_op.current_call_idx;
call_loc_on_track_idx = data.track.call_loc_idx_on_track_interp;  % call emission location in terms of idx of interpolated track
bat_loc_at_call = data.track.track_interp(call_loc_on_track_idx(curr_call_idx),:);
 
% Update call and track idx information
set(handles.edit_call_idx1,'String',num2str(curr_call_idx));
set(handles.edit_call_idx2,'String',['/',num2str(length(data.mic_data.call_idx_w_track))]);

% Plot track and call locations
axes(handles.axes2);
corder = get(gca,'colororder');
% plot all tracks
h_all_track = plot3(data.track.track_interp(:,1),...
                    data.track.track_interp(:,2),...
                    data.track.track_interp(:,3),'color',corder(1,:));
hold on
% plot call location
h_all_call = plot3(data.track.track_interp(call_loc_on_track_idx,1),...
                   data.track.track_interp(call_loc_on_track_idx,2),...
                   data.track.track_interp(call_loc_on_track_idx,3),'r.','markersize',8);
% plot mic location
h_mic_loc = plot3(data.mic_loc(:,1),data.mic_loc(:,2),data.mic_loc(:,3),'ko');
h_mic_num = text(data.mic_loc(:,1),data.mic_loc(:,2),data.mic_loc(:,3),num2str([1:data.mic_data.num_ch_in_file]'));
% plot current call location
h_curr_call = plot3(bat_loc_at_call(1),bat_loc_at_call(2),bat_loc_at_call(3),'o','color',corder(1,:));

grid on
axis equal
xlabel('X (m)');
ylabel('Y (m)');
zlabel('Z (m)');

% Save plot handles
handles.h_all_track = h_all_track;
handles.h_all_call = h_all_call;
handles.h_curr_call = h_curr_call;
handles.h_mic_loc = h_mic_loc;
handles.h_mic_num = h_mic_num;

% Save this GUI handles
setappdata(0,'track_gui_handles',handles);

% Choose default command line output for view_track_gui
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes view_track_gui wait for user response (see UIRESUME)
% uiwait(handles.figure1);



% --- Outputs from this function are returned to the command line.
function varargout = view_track_gui_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;



function edit_call_idx1_Callback(hObject, eventdata, handles)
% hObject    handle to edit_call_idx1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_call_idx1 as text
%        str2double(get(hObject,'String')) returns contents of edit_call_idx1 as a double


% --- Executes during object creation, after setting all properties.
function edit_call_idx1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_call_idx1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_track_idx1_Callback(hObject, eventdata, handles)
% hObject    handle to edit_track_idx1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_track_idx1 as text
%        str2double(get(hObject,'String')) returns contents of edit_track_idx1 as a double


% --- Executes during object creation, after setting all properties.
function edit_track_idx1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_track_idx1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in checkbox_top_view_bp.
function checkbox_top_view_bp_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_top_view_bp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_top_view_bp
plot_bat_mic_vector;