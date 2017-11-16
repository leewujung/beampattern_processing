## Overview
This set of Matlab code is written for processing and reconstructing beampatterns emitted by animals. The code includes a processing component and a GUI that allow users to go through each of the calls for visual inspection and excluding/including specific channels if necessary.

## Data Processing

### Batch processing input file
The entry point of the code is `batch_bp_proc_example.m` in which an example is given to set up the various parameters and filenames for feeding into the main processing function `bp_proc`. The paths to folders with the various types of data (video, acoustic, calibration) are specified in the next section. The set up for path and some global parameters are:
```matlab
username = getenv('username');
pname = ['C:\Users\',username,'\Dropbox\0_ANALYSIS\bp_processing'];  % base path
fname = 'rousettus_20150825_file_match.xlsx';  % spreadsheet containing all matching files of different types
trial_to_proc = 5:28;   % row index of the trials to process in the spreadsheet above
chk_indiv_call = 0;     % whether to display the cut-out section for each call/channel
track_cut_idx = 1:800;  % index for frames with acoustic data
save_dir = ['C:\Users\',username,'\Dropbox\0_ANALYSIS\bp_processing\proc_output_rousettus_new'];  % path to save processing output
if ~exist(save_dir,'dir')
    mkdir(save_dir);
end
load(['C:\Users\',username,'\Dropbox\0_CODE\beampattern_processing\bpf30.mat']);  % filter use only when detecting Rousettus clicks
``` 

Below are other user-specified parameters needed for beampattern processing:
```matlab
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
data.track.fs = 200;   % video frame rate [Hz]
data.track.smooth_len = 10;  % number of points used to smooth tracks
data.track.head_aim_est_time_diff = 50; % [ms] time span between points used for estimating head aim from animal position
data.track.head_n_prescribe = [0,0,1];  % precribed head normal vector (only used in 1-marker case, ignored in 3-marker case)

data.param.tempF = 75;  % temperature [deg F]
data.param.humid = 50;  % humidity (relative in %)
data.param.extract_call_len = 5;        % [ms] length of recording section isolated from around a given call
data.param.tolernace = 2;               % tolerance for call misalignment, make it bigger when mic location uncertainty is large
data.param.tukeywin_proportion = 0.25;  % proportional of tukeywin for call tapering
data.param.dura_flag = 0;   % 1-use duration marked in the mic detect file (FM bats)
                            % 0-use default detection stuff (Rousettus)

% Below are only used when data.param.dura_flag=1, such as for Rousettus clicks
data.param.call_short_len = 0.5;        % [ms] desired length of extracted call
data.param.call_portion_front = 0.2;    % proportion of extracted call before the peak
data.param.click_th = 0.1;              % threshold for extracting click ranges
data.param.click_bpf = bpf30;           % bandpass filter for use in click range estimation

data.param.mic_floor_idx = [4,5,7,17,24,26,27];  % index of mics on floor, used to estimate head normal vector
% in missing marker scenarios for 3-marker case
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
```


### Folder structure
The various types of data should be placed into specific structures to be loaded by the code. Here's how the folder structure should look like:

![folder structure](/img/folder_structure.png "Folder structure")

There should be a spreadsheet (`SPECIES_DATE_file_match.xlsx` in the above image) containing a list of matching files that specifies which combination of audio, animal position, and mic-related files are from the same experimental trial and should be processed together. The different data types and corresponding folders are explained in the following section. Below is an example showing how things go in the file matching list:

![file matching list](/img/file_matching_list.png "List of matching files")


### Data types
Data required for beampattern processing are the audio data, animal location, mic calibration/receiving beampattern, and mic info files that contain the other misc info from the experiment. The example files used here are from an experiment with:
* Audio recording: 34 channels, 4 second long, sampled at 250 kHz
* Three-dimensional bat position recorded at 200 Hz also for 4 seconds
* The audio and bat position recordings are synchronized using a stop-trigger, meaning the data saved were the 4 seconds *before* the trigger.

#### Audio data (folder: mic_data and mic_detect)
Audio data are MAT files with fields `sig` and `fs`. `sig` contains microphone recordings with data from each channel stored in each column. `fs` contains the sampling frequency of the audio data in Hz. For example, the example mentioned above would look like this when loaded into Matlab:
```
  Name            Size               Bytes  Class     Attributes

  fs              1x1                     8  double              
  sig       1000002x34            264000528  double    
```

#### Animal position (folder: bat_pos)
Animal positions are MAT files containing [x,y,z] positions of markers mounted on the animal's head or just the animal's rough position. The position of each marker should be stored in different cells. For data set in which there are 3 markers (explained below), the structure should look like this when loaded into Matlab:
```
  Name            Size            Bytes  Class     Attributes

  bat_pos         1x3             59376  cell                
  markers         1x3               360  cell             
```
in which
```
bat_pos = 
    [800x3 double]    [800x3 double]    [800x3 double]
```
and 
```
markers = 
    'Tip'    'Left'    'Right'
```
This format was used in my experiment where three markers on a triangular frame was mounted on the animal's head for monitoring the head aim and head roll during recording. The code calculates the head aim and head normal vectors based on these marker positions (in subfunction `load_bat_pos`). For video frames with missing markers, the code uses the trajectory differential as the head aim vector in x-y plane (with z=0), and calculates a head normal that is orthogonal to all the mics on the floor (mics identified by `data.param.mic_floor_idx` in the batch processing input file). This orthogonal head normal vector can also be derived using all the mics placed as a horizontal 1D array on the wall.

For data set with only one position at each time point, the structure would then be:
```
  Name            Size            Bytes  Class     Attributes

  bat_pos         1x1             19792  cell                
```
The code calculates the head aim vector based on the trajectory differential. The user should supply the head normal vector by setting `data.track.head_n_prescribe` in the batch processing input file.

Note the trajectories are smoothed when loaded. The length of smoothing is set by `data.track.smooth_len` in the batch processing input file. The time span over which the trajectory differential is used for estimating the head aim is set by `data.track.head_aim_est_time_diff` in [sec].


#### Mic info (folder: mic_info)
This is a MAT file containing miscellaneous information about the mics during experiment. The content should look like this when loaded into Matlab:
```
  Name           Size            Bytes  Class     Attributes

  mic_gain      34x1               272  double              
  mic_loc       34x3               816  double              
  mic_vec       34x3               816  double              
  mic_vh        34x1                 0  double              
```
Here, `mic_gain` are the gains applied to the mic channels *in dB*, `mic_loc` are [x,y,z] locations of the mics, `mic_vec` are [x,y,z] vectors representing the base-to-tip direction of the mics. `mic_vh` is only used if the mics are arranged in a cross configuration. In this case the array entry is 1 if a particular mic belongs to the horizontal line and 0 if it belongs to the vertical line, and 0.5 to be considered in both (e.g. center mic). `mic_vh` is empty if the mics are arranged in a irregular, non-cross configuration. This variable will be used for plotting in the GUI.

#### Mic receiving beampattern (folder: mic_bp)
This is a MAT file containing information about the receiving beampattern of each of the mics. The content should look like this when loaded into Matlab:
```
 Name        Size                  Bytes  Class     Attributes

  bp         91x101x34            2499952  double              
  freq       91x1                     728  double              
  src         1x1                 1563678  struct              
  theta       1x101                   808  double     
```
Here, `bp` stores the receiving beampattern at `nf` frequencies listed in `freq` across `nt` angles (degrees) listed in `theta`, for `nm` mic channels (in the example `nm=34`). Therefore, the dimension of `bp` is `[nf x nt x nm]`. The values in `bp` are normalized to the axis of the mic at `theta=0` in dB scale. `src` contains information about from which files these parameters come from. This variable is not used in the processing. Note that the sequence of the `nm` channels here should be ordered according to the sequence of the mic locations.

#### Mic sensitivity (folder: mic_sens) ####
This is a MAT file containing information about the sensitivity of each of the mics. The content should look like this when loaded into Matlab:
```
  Name        Size            Bytes  Class     Attributes

  freq      129x1              1032  double              
  sens      129x34            35088  double              
  src         1x1             96210  struct  
```
Here, `sens` stores the sensitivity at `nf` frequencies listed in `freq` for `nm` mic channels. Therefore, the dimension of `sens` is `[nf x nm]`. The values in `sens` are in dB re 1V/Pa. `src` contains information about from which files these parameters come from. This variable is not used in the processing. Note that the sequence of the `nm` channels here should be ordered according to the sequence of the mic locations.

## Beampattern visualization GUI ##
The GUI displays the results of beampattern processing and allows the user to manually modify a few important parameters. The GUI is called by `bp_check_gui.m`. The GUI looks like:

![Beampattern GUI](/img/bp_gui_main.png "Beampattern visualization GUI main window")

The left plot axis shows the interpolated beampattern across the whole globe in filled contour form. The right plot axis shows the interpolated beampattern within a limited azimuth range in full colormap. In this axis the contour denotes the -3dB contour from the peak intensity. The rightmost panel displays the time series (with mic-animal distance compensated) from all channels. The uncompensated time series are shown in gray, the compensated time series are shown in blue, and the included section for beampattern analysis are shown in red.

The various buttons/check boxes are explained below:

- `Pick proc'ed file`  Select the processed file to view
- `Save :)`  Save all edited results through this GUI
- `<< Previous` and `Next >>` Navigate through the calls; the white space right next to **Call #** can also be used to directly jump to another call
- `Edit call section`  Brings up another window and allow user to manually change the included section of each channel to be used for beampattern calculation
- `Beampattern freq (kHz)`  Frequency at which the beampattern is to be visualized
- `Mic config`  Beampattern display according to mic configuration. It can be `2D` in which mics that are randomly located in the azimuth-elevation planes are used to do the interpolation. Alternatively, for `Cross` the mics are arranged in a vertical and horizontal configuration and the beampatterns are displayed as two polar plots in azimuth and elevation separately.
- `Log/Linear`  The scale at which the polar plots in `Cross` config will be displayed. In `2D` the beampattern is always displayed in the log scale.
- `Interpolation`  The interpolation scheme used for the displayed `2D` beampattern. `Natural neighbor` is default from MATLAB, and `Radial basis fcn` is a very useful interpolation method in computer graphics.
- `Head aim from marker`  If checked it means the *local* animal coordinate (head aim and head normal vectors) displayed are derived from >1 markers mounted on the animal's head. Otherwise these vectors are derived from the trajectory and user input.
- `Good call`  If this call is good to be included in follow-up analysis. This can be determined manually here in the GUI or programatically by specifying conditions to modify `data.proc.chk_good_call`.
- `Normalize bp plot`  If the displayed `2D` beampattern plot is normalized. I found it easier to visualize single beampattern characteristics when it's normalized, but for inter-call/click comparison it's easier without normalization.
- `Channels to be excluded`  Manually specify the mic channels to be excluded from beampattern analysis
- `cmin` `cmax`  The range of color axis in the `2D` beampattern display
- `Peak dB`  The maximum intensity measured in all mics, in dB scale

The animal position at which the call was emitted is plotted along with the entire trajectory in a separate window:

![Beampattern GUI track view](/img/bp_gui_track_view.png "Animal orientation/track and call locations")

The *local* animal coordinate is plotted as vectors:

- red: head aim
- orange: head normal
- purple: cross-product of the above two vectors, pointing to the right of the animal
The azimuth and elevation of the beampattern representation shown in the main GUI window are calculated based on this *local* coordinate.
