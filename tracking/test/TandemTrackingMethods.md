# Diacamma Tandem Tracking Methdos
#### N. Mizumoto

### 1. Tracking using FastTrack  
- File -> Open -> File name  
- Image Options  
  - Background (average, 100)
  - ROI
  - Binarization + Operation + Detection (Need some experiments to determine the parameter values), so that tandem pair is tracked as if single individual.
- Tracking Options  
  - Tuning parameters  
  - Preview  
  - Track  
- Tracking Inspector (manual correct)  
  - Name focused tandem as “0” and
  - go through the whole video
  - The data will be autosaved. The results can be found in the directory “Tracking_Result_(Filename)”. The “Export” button will generate a video with annotations.  
  - Get tracking.txt in the result directory and rename it.
  
### 2. Convert tracking.txt into .csv file for python analysis  
- Place .txt files at [data/input/](./tracking/test/data/input) (multiple files can be processed at the same time)
- Run [TrackingConverter.R](./tracking/test/TrackingConverter.R)
- The csv file will be found in [data/output/](./tracking/test/data/output)

### 3. Make a video with focused tandem using Python
- Place original video files at [VideoPrepPython/videos](./tracking/test/VideoPrepPython/videos) (multiple files can be processed at the same time)
- Place csv files at [VideoPrepPython/coordinate](./tracking/test/VideoPrepPython/coordinate)
- Follow the python program in [VideoPrepPython/coordinate](./tracking/test/VideoPrepPython/coordinate)
- The results will be found as VideoPrepPython/videos/*_extract.mp4

### 4. Use UMATracker
- Filter generator
  - Create background 
  - Rectangle selection
  - BGRToGray, Threshold
- Tracking
  - Drag and drop video file and filter file
  - num of object = 2
- Tracking Corrector

