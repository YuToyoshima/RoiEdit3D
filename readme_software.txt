Readme for using roiedt3d_25_12.


Written by Yu Toyoshima (ytoyo@bs.s.u-tokyo.ac.jp), 2020/Jan/28.


1. The purpose of the programs

RoiEdit3d has following functionalities:
+ automatic detection of cell nuclei in  2D/3D/4D image(s) 
+ automatic annotation for the cell nuclei (when the atlas was provided)
+ automatic tracking the cell nuclei
+ visulaizes the results and provides an efficient way for correction


2. Environmental requirement 

The programs were developped under Windows 7 Pro SP1 64bit.
The executables for Windows 64bit environment will be located in binary directory (bin/).

For other environments, run the code in matlab.
The instructions for compiling the program from source will be bundled (compilation.txt).


3. Installation 

This section is for using bundled executables on Windows.


3-1. Install Matlab Compiler Runtime (R2017a)

Download the Matlab Compiler Runtime (R2017a) for Windows 64bit from Mathworks web site
(http://jp.mathworks.com/products/compiler/mcr/) and install.


3-2. Locate required library

Run RoiEdit3D as an Administrator.
If everything goes well, you will see the instuction for restarting the program.

If somthing goes wrong, you will see messages like as "Copy failed", and
you should locate guava (java library) manually.
Please see the section 5 (Locate guava) of the instructions for building the program from source.




///////////////////////////////////////////////////////////////////

Please note that the following descriptions are for older version.
Descriptions for current version is under construction.

Please see the tutorial materials in figshare:
http://dx.doi.org/10.6084/m9.figshare.8341088

///////////////////////////////////////////////////////////////////







4. How to use peak_detection

Drug and Drop *.tif file(s) to the executable, the program will detect cell nulcei and 
output thr results for each image file.
Peak_detection requires a parameter file that reflects the character of the input image.
Peak_detection_12.exe searches peak_detction_12_param.txt in the same directory for the image.

Peak_detection will make following four files as an output;
+ (image_name).mat contains the information of positions and shapes of the detected nuclei. 
+ (image_name)_aligned.tif is the input image that the parallel shift was cancelled.
+ (image_name)_blurred.tif is the image for the nuclei detection and gaussian mixture optimization.
+ (image_name)_cum2tmp.mat contains the information for parallel shift.
The former 3 files are required for RoiEdit3D.


5. How to use RoiEdit3D

Drug & Drop an (image_name).mat to roiedit3d_25_12.exe. 
Then the program will be started with following windows;
+ Command prompt
+ FIJI/ImageJ main window
+ RoiEdit3D main window (titled ÅeRoiEdit3D)
+ Aligned image (with orthogonal views for z-sliced image)
+ Customized ROI Manager


5-1. RoiEdit3D main window

+ Main panel
  + Load / Save: Loading and saving the result of detection and correction of nuclei.
  + Visualize radius: Scaling the size of ellipsoidal ROIs.
  + Dist. threshold (um): Distance of the threshold between centers of ROIs.
    The unit is micro meter.
  + Logarithm intensity: Show the image in logarithmic scale.
  + MaxProjZ: Show max projected image for Z-slice.

+ ROI_disp panel
  + ROI: Toggle display ROIs.
  + Name: Toggle display ROI names.
  + (Upper pulldown menu): Specify unselected ROI color.
  + (Lower pulldown menu): Specify selected ROI color

+ Filtering panel
  Filter the ROIs displayed in Customized ROI Manager window.
  + < distance threshold: Display ROIs that have neighbor within Dist. threshold.
  + OutOfArea: Display ROIs whose positions are out of area of original blob.
    Valid only for the direct output of PeakDetection. 
  + Recovered: Display ROIs that were split based on curvature.
    Valid only for the direct output of PeakDetection. 
  + Unchecked: Display ROIs that were not ÅecheckedÅf in the leftmost column of Customized ROI Manager.

+ Optimize panel
  + Current: Optimize parameters of ROIs in the image of displayed time frame. 
             Please use this function after manual correction (add/remove) of ROIs.
  + Åe>Åf  : Copy ROIs in the image of current time frame to that of the next time frame 
             and optimize parameters. (i.e. Tracking ROIs by one frame after).
  + Åe>>Åf : Repeat Åe>Åf automatically for all time frames (i.e. Tracking ROIs in ascending manner).
  + Åe<Åf  : Copy ROIs in the image of current time frame to that of the previous time frame 
             and optimize parameters. (i.e. Tracking ROIs by one frame before).
  + Åe<<Åf : Repeat Åe<Åf automatically for all time frames (i.e. Tracking ROIs in descending manner).

+ Rotation panel
  Rotate the images and ROIs by angle of 180 degree with the X, Y, or Z axis as a center.


5-2. Customize ROI Manager window

Each line corresponds to each ROI in the image of one time frame.
ROI(s) can be selected by left click (modification by shift or ctrl is valid).
Selected ROIs can be removed by pressing delete.
Modify the names and parameters of ROIs by editing the content in the cell of the table.
ROIs can be sorted and/or filtered by clicking the column headers.


5-3. Aligned image window

The aligned image of the input is shown.
For Z-sliced image, XZ section and YZ section of the image will be shown as trihedral figures 
using Orthogonal Views plugin. 
The Orthogonal Views plugin used in RoiEdit3D 
(Fiji/ImageJ main window -> Plugins -> MyOrthogonal Views) 
is the modified version of standard Orthogonal Views plugin in Fiji/ImageJ 
(Fiji/ImageJ main window -> Image -> Stack -> Orthogonal Views).

+ Left click: Update the XZ and YZ section (only for Orthogonal Views)
+ Left click with Shift: Search and select a ROI near the specified position.
+ Left click with Ctrl: Put a new ROI at the specified position.
+ Left click with Ctrl & Shift: Search and remove a ROI near the specified position.
+ Space key: Toggle display ROIs.
+ Space key with Shift: Toggle display ROI names
+ Space key with Ctrl: Optimize ROIs in the image of current time frame.


6. Note

+ RoiEdit3D may freeze when the operations are too fast. Please stop operating 
  when the buttons and checkboxes in RoiEdit3D main panel are displayed in grayout state. 

+ RoiEdit3D save the results automatically when the operations are specified.
  When the program stopped abnormally, the results might be recovered by 
  loading (working_mat_file_name)_tmp.mat.

+ For Z-sliced image, all the ROIs in the frame may be displayed in all the slices.
  If so, please check "Associate Show All ROIs with slices" through "More >> optionsÅc" 
  in the ROI Manager of the Fiji/ImageJ 
  (Fiji/ImageJ main window -> Analyze -> Tools -> ROI ManagerÅc).




7. Copyright notice

This software utilizes many functionalities made by other people such as
Matlab, Fiji/ImageJ, guava, MorphoLibJ, dftregistration, MinMaxFilter.
The copyright and related rights of the functionalities have been owned by 
the copyright holder. The rights of the other part of the software have 
been owned by Yu Toyoshima. 


This project is licensed under the terms of the MIT license.

Copyright (c) 2016 Yu Toyoshima

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
THE SOFTWARE.
