## Overview ##
This set of Matlab code is written for processing and reconstructing beampatterns emitted by animals. The code include a processing component and a GUI that allow users to go through each of the calls for visual inspection and excluding/including specific channels if necessary.

## Processing Code ##
#### Batch processing starting file ####
The entry point of the code is `batch_bp_proc_example.m` in which an example is given to set up the various parameters and filenames for feeding into the main processing function `bp_proc`. 
Paths to the folders in which the audio data, animal location, mic calibration, and mic receiving beampattern files are specified in the section below. The path set up for various types of data rare
```matlab
username = getenv('username');
pname = ['C:\Users\',username,'\Dropbox\0_ANALYSIS\bp_processing'];
fname = 'eptesicus_20150824_file_match.xlsx';
trial_to_proc = 26:34;
chk_indiv_call = 1;
track_cut_idx = 1:800;
save_dir = ['C:\Users\',username,'\Dropbox\0_ANALYSIS\bp_processing\proc_output_eptesicus'];
if ~exist(save_dir,'dir')
    mkdir(save_dir);
end
load(['C:\Users\',username,'\Dropbox\0_CODE\beampattern_processing\bpf30.mat']);
``` 
### Data types ###
Data required for beampattern processing are the audio data, animal location, mic calibration/receiving beampattern, and mic info files that contain the other misc info from the experiment. Each type of the data are explained below. The example files are from an experiment with:
* Audio recording: 33 channels, 4 second long, sampled at 250 kHz
* Three-dimensional bat position recorded at 200 Hz also for 4 seconds
* The audio and bat position recordings are synchronized using a stop-trigger, meaning the data saved were the 4 seconds *before* the trigger.

#### Audio data ####
Audio data are `MAT` files with fields `sig` and `fs`. `sig` contains microphone recordings with data from each channel stored in each column. `fs` contains the sampling frequency of the audio data in Hz. For example, the example mentioned above would look like this when loaded into Matlab:
```
  Name            Size               Bytes  Class     Attributes

  fs              1x1                     8  double              
  sig       1000002x33            264000528  double    
```

#### Animal position ####
Animal positions are `MAT` files containing information from the three markers mounted on the animal's head. It should look something this when loaded into Matlab:
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
This format is used because in my experiment three markers on a triangular frame was mounted on the bat's head for monitoring the head aim and head roll during the recording. For data set with only one position at each time point, make the content of all three cells in `bat_pos` identical. For data set with two positions (e.g., front and back markers), put the position of the front marker in the first cell and the position of the back marker in the second and third cell in `bat_pos`. The code will extract head aim and head roll based on these data.

* **NOTE**: add description of the processing procedure for extracting head aim/roll, including smoothing and interpolation here.

#### Mic info ####
This is a `MAT` file containing miscellaneous information about the mics during experiment. They should look like this when loaded into Matlab:
```
  Name           Size            Bytes  Class     Attributes

  mic_gain      33x1               272  double              
  mic_loc       33x3               816  double              
  mic_vec       33x3               816  double              
  mic_vh         0x0                 0  double              
```

## GUI display ##
