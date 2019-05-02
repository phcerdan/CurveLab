****************************************
* CurveletUtils                        
* version: 0.1 (Dec 17, 2007)          
* 	first version
* version: 0.2 (Feb 11, 2008)
*	major changes to EdgeDetect GUI with many options/thresholds to be set
*	extraction of sheets and edges, and measurements of length and area
* 	crvlt_getmagnitude mex file extracts total magnitudes
* version: 0.3 (Jul 2, 2008)
* 	added batch mode to EdgeDetect GUI and CViewer GUI
****************************************

The CurveletUtils package contains MATLAB utilities to use curvelets
and display results in GUIs. 
Developed on Mac OS X with MATLAB (R2008a). Other platforms have
not been tested, but it should work.

The software is free for non-commercial purposes, but see section 5 on
citation information. For commercial use, please contact the address below.

Copyright 2008 Tobias Gebaeck, ETH Zurich
For info, suggestions and bug reports contact: tobias.gebaeck@inf.ethz.ch

1. Contents

The CurveletUtils contain three MATLAB GUI:s for visualising and experimenting 
with curvelets:
* cviewer: general tool for thresholding and visualisation of curvelet 
  transformed images, predefined shapes and single curvelets.
* cmoviemaker: a tool to make movies out of several images, applying a curvelet 
  or wavelet transform to each frame.
* edgedetect: a tool for experimenting with the crvlt_extractdirs utility, 
  finding edges in images using the curvelet transform. 

In addition, there are a few useful tools in the utils/ directory, and some 
example programs named crvlt_*.m.


2. Installation

* Download and Install CurveLab (www.curvelet.org) (latest tested
  version is 2.1.2) according to the instructions (including installing
  FFTW ver 2 from www.fftw.org), including the matlab mex-functions
  (write "make matlab" in the CurveLab directory).
* Edit the makefile in the 'src' subdirectory of your CurveletUtils
  directory and set the directories to the proper CurveLab and FFTW
  directories, as well as the path to MATLAB's 'mex' utility.
* Run "make"
* Edit "setpath.m" and set the proper paths to CurveLab and CurveletUtils
* Start MATLAB, change directory to your CurveletUtils directory, and
  type "cviewer" at the prompt to start the CViewer GUI, or
  "edgedetect" to start the EdgeDetect GUI, or start using 
  the other tools.


3. Quick start

3.1 CViewer

To view the reconstruction of an image using only levels 3 and 4 of
the curvelet transform, do the following:

* Start CViewer by writing "cviewer" at the MATLAB command prompt.
* Ensure that the "Display mode" radio button is set to "Image". 
* Click the "Show an image" drop-down list and select "Load
   image...". In the dialog box that appears, browse to the image you
   want to display and click "Open". It is advisable to use an image
   with resolution less than 1000 x 1000 pixels to ensure reasonable
   computation times when trying out the GUI.
* In the "Domain" drop-down next to the bottom image, select
  "Reconstr. - physical space".
* In "Levels to show", select levels 3 and 4 (using Shift to select
  multiple levels).
* Click the "GO!" button. After a while, the reconstructed image will
   appear as the bottom image.


3.2 EdgeDetect

To detect edges in an image using the curvelet-based scheme, do the
following: 

* Start EdgeDetect by writing "edgedetect" at the MATLAB command prompt.
* Select "Open file..." in the "File" menu and select an image
  file. It is advisable to use an image with less than 1000 x 1000
  pixels to ensure reasonable computation times.
* Select for example level 4 in the "Levels" list box.
* Check the "Edges" check box at the bottom.
* The detected edges on the specified level should now show up in
  green. If no edges are shown, try lowering the "Canny thresholds"
  on the left (the left value must be lower than the right value),
  until some edges are seen.


4. Features

4.1 CViewer

CViewer uses the fast discrete curvelet transform with default number
of levels and directions, and with wavelets on the finest level.


4.1.1 Basic features

The CViewer GUI has the following features:

* View single curvelets: Select "single curvelet" in the "Display
  mode" box, select desired levels and directions on the right and click
  "GO!"

* View predefined images and their curvelet transforms: Select a shape
  from the "Show an image" menu, select levels and directions on the
  right, click "GO!", select desired view in "Domain" and "Values"	
  menus, pan and zoom using the tools in the toolbar

* Load images from files: Select "Load image..." from the "Show an
  image" menu, then proceed as above.

* Threshold curvelet coefficients: On the right, select Threshold and
  enter a value in the box to use a fixed threshold, select Number and
  enter the number of the largest coefficients to keep, or select Percentage and
  enter a fraction between 0 and 1 to keep the largest fraction of the nonzero
  coefficients (applied after the level and direction selection). Then
  click "GO!" to apply the extraction.

* Link the axes when zooming: Click the up/down arrow in the toolbar
  to link/unlink the top and bottom axes.

* Show positions of nonzero curvelet coefficients: Select an item in
  "Show curvelet positions" menu (only useful if thresholding is
  applied). Different levels get different colors.

* View directions of curvelets: Click the rightmost image in the
  toolbar, then click on the top axes (when an image is loaded) to see
  the directions and magnitudes of the curvelet coefficients at that	
  point. Select levels on the right. Choose polar/rectangular plot,
  running average, summing across levels at will.

* Export data: click this button to export the current data to the
  main window (in the variable 'cviewdata')

* Show axes in new window: Click "Show in separate window" to 
  open a new figure and show the current plot in that figure. From here
  you can edit, print and export the plot.

* Make movies: Select the menu item "Tools->Make movie..." to open the
  CMovieMaker GUI, which enables you to make a movie out of a series
  of images, showing curvelet reconstructions, errors, wavelet
  reconstructions etc with similar features as above. Choose layout of
  the figure on bottom right, click Preview to see what it looks like,
  and click "Make movie!" to create a movie and save as an AVI-file.


4.1.2 Description of each button/selection in CViewer

4.1.2.1 Buttons / drop-down menus

Display mode: Image/Single curvelet. The "Image" mode lets the user
display an image or shape using the "Show an image" drop down, and
view reconstructions using selected curvelet levels and directions,
and threshold the curvelet coefficients. The "Single curvelet" mode
lets the user select curvelet levels and directions on the right and
view the resulting curvelet image (obtained by setting one
curvelet coefficient in the center of the image on each selected level
and direction to 1, and computing the inverse curvelet transform).

Show an image: This drop-down box lets the user select a pre-defined
shape for display and for viewing its curvelet decomposition. When
selecting "Load image...", a dialog box appears for opening an image
file. Any image format recognized by MATLAB can be opened, but only
gray-scale images and images with one non-zero color channel are
treated correctly (since the curvelet transform only works for one
color channel at a time).

Smooth shapes: Checking this checkbox applies some smoothing to the
shapes selected in the "Show an image" drop-down.

Colormap: provides a set of colormaps for displaying the images.

Invert colormap: when checked, inverts the colormap,
i.e. white-to-black instead of black-to-white for the gray colormap

Show curvelet positions: Selecting one of the options here displays
the position of the active curvelet positions, i.e. the curvelets
selected using the "Levels to show" and "Directions to show" list
boxes, and that pass any applied thresholding. Selecting "Large dots"
or "Small dots" shows every curvelet as a dot, with different color
for different detail levels, while selecting "Scaled arrows" or
"Same-size arrows" shows a small arrow instead, indicating the
direction of the curvelet. The arrow points in the direction of the
short axis of the curvelet in physical space, i.e. in the oscillating
direction. No distinction is made between directions separated by 180
degrees, and only one of these directions is shown.

Levels to show: In this list box one or several detail levels of the
curvelet transform may be selected for inclusion in the reconstructed
image. Use Shift and/or Ctrl/Command to select multiple
levels. Changes appear when the "GO!" button is clicked. Only the
levels that exist in the curvelet transform of the current image are
displayed. 

Directions to show: In this list box, the curvelet directions included
in the reconstruction may be selected, just as for the "Levels to
show" list box. The directions available on the coarsest selected
curvelet level are shown, and curvelets on finer levels with more
directions are included if they belong to a subdivision of the
selected direction. Thus, if levels 2 to 4 are selected, and 16
directions are shown in the "Directions to show" list box, but level 4
has 32 directions, then selecting direction 1 in the list box will
include directions 1 and 2 on level 4.

Symmetric directions: Checking this box, forces the selection of
directions to be symmetric, i.e. when one direction is selected, the
direction opposite to it (shifted by 180 degrees) is also included. It
is often wise to enable this option, because it ensures that
reconstructed images will be real.

Domain: This drop-down box determines what will be shown in the image
axes next to it. The choices are "Original image", "Reconstructed
image", "Error" (reconstructed minus original) and "Coefficients" (the
value of the curvelet coefficients on the coarsest selected level and
first selected direction for all space points - this is not equal to
the reconstruction). Original, reconstructed and error images can be
viewed either in physical space or in Fourier space. And the
coefficients can be viewed either for the original image or for the
reconstructed (and thresholded) image.

Values: This determines which values are shown in the image axes next
to it - real part, imaginary part, absolute value or logarithm of
absolute values.

Thresholding: In this box, the thresholding of curvelet coefficients
for reconstruction can be specified. Selecting "Threshold" in the
drop-down keeps only coefficients with absolute value larger than the
value specified in the edit box. Selecting "Number" keeps the largest
coefficients, but only the number of coefficients specified in the
edit box. The "Percentage" option keeps the largest fraction specified
in the edit box, i.e. writing 0.1 keeps the 10% largest coefficients.

Show in separate window: Clicking these buttons opens a new figure
window and displays the image currently displayed, but in a new window
- useful for comparisons, and for printing and saving images.

Export data: Exports current internal program data to a variable
"cviewdata" in the matlab workspace.

l^2 error and PSNR: Show the error of the reconstruction as l^2-error,
||im - rec||_2, and as
 PSNR (= 20 * log10((max(im) - min(im)) / ||im - rec||_2) 


4.1.2.2 The Toolbar

The toolbar contains the usual MATLAB Print, Zoom and Pan tools, but
also the following tools:

Link axes (the up-down arrow): Links the zoom and pan operations in
the two axes (useful for comparing original and reconstructed images)

Plot directions (the red star-like shape): Turns the bottom axes into
a polar plot, when the user clicks in the top image (displaying the
original or reconstructed image). The polar plot shows the magnitude
of curvelet coefficients in different directions for the selected
detail levels at the position clicked in the image. By the check-boxes
that appear to the left, the user can turn the polar plot into a
cartesian plot, specify a sliding average over a specified number of
neighboring curvelets before plotting, and choose to sum coefficients
on all the selected levels.

Norm of directions (the multi-colored disk): Displays the l^2-norm of
the curvelet coefficients at all levels and directions, computed
separately at each level and direction.


4.1.2.3 Menu items

The Tools menu contains two menu items: "Make movie" and "Batch
analysis...". The "Make movie" item opens the CMovieMaker window that
enables the user to create movies of a series of images, where
curvelet and/or wavelet transforms and thresholding can be applied to
each image (with similar options like in the CViewer GUI).
The "Batch analysis" item opens the "Batch window", where the user can
select a directory, and all images in this directory are loaded and
the current level, direction and thresholding options selected in the
CViewer window are applied to each image. For each image, a mat-file
is created, containing the same data as when clicking "Export data",
and the reconstructed image is saved in a subdirectory named Output,
and with a filename appening "_analyzed" to the original filename.



4.2 The EdgeDetect GUI

4.2.1 Basic features

The edgedetect GUI implements the algorithm presented in the paper
"T. Geback, P. Koumoutsakos; Edge detection in microscopy images using
curvelets; 2008"
for edge detection using the fast discrete curvelet transform.

The algorithm extracts a directional field with direction and
magnitude at each point given by the largest curvelet coefficient on
the chosen level. The resolution is given by the finest chosen level of curvelet
coefficients. The Canny edge detector is applied to trace along ridges
of the field magnitude to mark the edges.

4.2.2 All buttons and boxes in EdgeDetect GUI

4.2.2.1 Buttons and drop-downs

The "Directional field extraction" group contains options for
selecting which curvelet levels to use, and how to include them.
The "Levels" list box lets the user select which detail levels to
use. Checking "All curvelets" uses curvelets on the finest level
instead of wavelets. "Circular average size" specifies the size of the
averaging over directions before extracting curvelet magnitudes (1
means going 1 step left and 1 right, meaning a moving average of size
3). "Curvelet extent" specifies in the same way how many neighboring
positions should be averaged over (1 here means 1 left, 1 right, 1 up,
1 down - total 5). "# directional fields" means how many directional
fields should be extracted and used for the edge detection, the first
one being the largest coefficients, the second field the second
largest coefficients, if separated from the first by the number of
directions given in "Direction spacing".
Finally, "Filter type" specifies if we use the magnitude of the
curvelet coefficient, or the real part (corresponding to an even
filter), imaginary part (an odd filter), or sum of the real and
imaginary parts.

The "Edge detection" group lets the user specify thresholds for the edge
detection. More information about the Canny thresholds may be found by
writing "help edge" in the Matlab window. "Apply tophat transform"
applies a tophat-transform (original - opening of image) using a disk
of the specified radius as the structural element.
The "Sheet threshold" is used for extracting "sheets", which are areas
where the total curvelet magnitude (summed over all directions) is
large. Finally, "Exclude sheets from edges" prevents the marking of
edges in any areas that are marked as sheets.

The "Postprocessing" group provides some opportunities for
manipulating the extracted edges. An "Extension length" larger than 0
extends the found edges by at most the specified number of steps,
moving in the direction of the field. If another edge is encountered,
the travelled path between is included as an edge, connecting the two
edges. The global extension threshold gives the fraction of the
maximum magnitude that a pixel must exceed to be considered for
inclusion, while the local threshold is the fraction of the current
pixel value (that we are moving from).
"Thin edges" applies the matlab command "bwmorph('thin')" to the
edges. "Dilate and thin" applies the "bwmorph('dilate')" and the
"bwmorph('thin')" commands. "Remove spurs" applies "bwmorph('spur')".
"Skeletonize" uses bwmorph('thin', Inf) on the edges mapped to the
original image before displaying them.

Details about the edge detection and post-processing can be found in
curvecanny_multi.m.

In the "Display" group, the user can choose to display the image or
the curvelet magnitude image as background to the edges, choose to
display edges and/or sheets and/or the directional field (with small
arrows indicating direction and magnitude). Finally, the user can
select which of the extracted fields to show (if more than one).

By clicking "Export data", the user exports the current data
(including field, edges, sheets, etc) to the Matlab command window in
the variable edgedata.


4.2.2.2 The toolbar

The toolbar contains the standard Matlab tools for zooming and
panning.


4.2.2.3 The menu

The "Open file..." item opens a dialog box where the user can select
an image file to open. All image formats supported by Matlab's imread
function are supported, plus the .dcm DICOM format.

The "Batch analysis..." item opens the batch window, where the user
can select a directory from which all image files are read and
analysed using the current settings in the main window. For each
image, a Matlab mat-file is written in a subdirectory name Output,
containing similar data as exported by the "Export data" button. The
analyzed image is also written, corresponding to the options for
display specified in the main window. Some statistics is also
displayed in the Matlab command window for each analyzed file.



4.3 Other files and functions

For information on the other files and functions included in
CurveletUtils, see each of the separate files. The utils/ subdirectory
contains some simpler utilities that makes working with the Curvelet
data structure a bit simples. The src/ subdirectory contains source
for the mex-files used by the GUIs.




5 Citation information

If you use the EdgeDetect GUI or the underlying algorithm for any
sientific publication, please cite the following paper:
T. Geback, P. Koumoutsakos; "Edge detection in microscopy images using
curvelets", 2008





