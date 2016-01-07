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

## Folder structure ##
The various types of data should be placed into specific structures to be loaded by the code. Here's how the folder structure should look like:

![folder structure](/img/folder_structure.png "Folder structure")

There should be a spreadsheet (SPECIES_DATE_file_match.xlsx in above image) containing a list of matching files that specifies which combination of audio, animal position, and mic-related files are from the same experimental trial and should be processed together. Below is an example file showing how things go:

![file matching list](/img/file_matching_list.png "List of matching files")

The different data types and corresponding folders they should go in are explained in the following section.

## Data types ##
Data required for beampattern processing are the audio data, animal location, mic calibration/receiving beampattern, and mic info files that contain the other misc info from the experiment. Each type of the data are explained below. The example files are from an experiment with:
* Audio recording: 34 channels, 4 second long, sampled at 250 kHz
* Three-dimensional bat position recorded at 200 Hz also for 4 seconds
* The audio and bat position recordings are synchronized using a stop-trigger, meaning the data saved were the 4 seconds *before* the trigger.

### Audio data ###
Audio data are MAT files with fields `sig` and `fs`. `sig` contains microphone recordings with data from each channel stored in each column. `fs` contains the sampling frequency of the audio data in Hz. For example, the example mentioned above would look like this when loaded into Matlab:
```
  Name            Size               Bytes  Class     Attributes

  fs              1x1                     8  double              
  sig       1000002x34            264000528  double    
```

### Animal position ###
Animal positions are MAT files containing [x,y,z] positions of markers mounted on the animal's head or just the bat's rough position. It should look something this when loaded into Matlab:
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
This format is used because in my experiment three markers on a triangular frame was mounted on the bat's head for monitoring the head aim and head roll during the recording. For data set with only one position at each time point, make the content of all three cells in `bat_pos` identical. For data set with two positions (e.g., front and back markers), put the positions of the front marker in the first cell and the position of the back marker in the second and third cell in `bat_pos`. The code will extract head aim and head roll based on these data.

* **NOTE**: add description of the processing procedure for extracting head aim/roll, including smoothing and interpolation here.
* **NOTE**: need to update code by adding a switch for dealing with different number of markers

### Mic info ###
This is a MAT file containing miscellaneous information about the mics during experiment. The content should look like this when loaded into Matlab:
```
  Name           Size            Bytes  Class     Attributes

  mic_gain      34x1               272  double              
  mic_loc       34x3               816  double              
  mic_vec       34x3               816  double              
  mic_vh        34x1                 0  double              
```
Here, `mic_gain` are the gains applied to the mic channels *in dB*, `mic_loc` are [x,y,z] locations of the mics, `mic_vec` are [x,y,z] vectors representing the base-to-tip direction of the mics. `mic_vh` is only used if the mics are arranged in a cross configuration. In this case the array entry is 1 if a particular mic belongs to the horizontal line and 0 if it belongs to the vertical line. `mic_vh` is empty if the mics are arranged in a irregular, non-cross configuration. This variable will be used for plotting in the GUI.

#### Mic receiving beampattern ####
This is a MAT file containing information about the receiving beampattern of each of the mics. The content should look like this when loaded into Matlab:
```
 Name        Size                  Bytes  Class     Attributes

  bp         91x101x34            2499952  double              
  freq       91x1                     728  double              
  src         1x1                 1563678  struct              
  theta       1x101                   808  double     
```
Here, `bp` stores the receiving beampattern at `nf` frequencies listed in `freq` across `nt` polar angles listed in `theta`, for `nm` mic channels (in the example `nm=34`). Therefore, the dimension of `bp` is `[nfx nt x nm]`. The values in `bp` are normalized to the axis of the mic at `theta=0` in dB scale. `src` contains information about from which files these parameters come from. This variable is not used int the processing.

#### Mic sensitivity ####
This is a MAT file containing information about the sensitivity of each of the mics. The content should look like this when loaded into Matlab:
```
  Name        Size            Bytes  Class     Attributes

  freq      129x1              1032  double              
  sens      129x34            35088  double              
  src         1x1             96210  struct  
```
Here, `sens` stores the sensitivity at `nf` frequencies listed in `freq` for `nm` mic channels. Therefore, the dimension of `sens` is `[nf x nm]`. The values in `sens` are in dB re 1V/Pa. `src` contains information about from which files these parameters come from. This variable is not used int the processing.

## GUI display ##
