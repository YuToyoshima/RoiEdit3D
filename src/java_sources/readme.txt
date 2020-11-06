HOW_TO_USE

Compile all *.java files with Fiji (ImageJ), Life-Line version, 2013 July 15.
Compress all *.class files to a zip file and rename it as 'roiedit3d_subclass.jar', then copy it to (FIJI_INSTALL_PATH)\plugins\.


Readme for roiedit3d_subclass

Written by Yu Toyoshima (ytoyo@bs.s.u-tokyo.ac.jp), 2020/Jan/28.
Yu Toyoshima @ Univ. of Tokyo



This library was written for handling 5-dimensional (xyzct) images in a GUI roiedit3d.
This library also works as a plugin in ImageJ/Fiji and can be used for reading image files including *.tif, *.dcv, *.sif, and *.mif.

*.dcv: A file format produced by High Speed Recording software (Hamamatsu Photonics) for sCMOS camera ORCA-Flash 4.0 series 
*.sif: A file format produced by iQ software (Andor) for EMCCD camera iXon Ultra series
*.tif: A file format produced by ImageJ/Fiji. Please note that a tiff file format allows several different types of encoding methods. ImageJ/Fiji and this library support only a part of the whole encoding methods.
*. mif: A intermediate file format. Human readable text. One can integrate multiple *.dcv, *.sif, and *.tif images into a single image file.

How to use this library as a ImageJ/Fiji plugin:
1. Put guava-18.0.jar and roiedit3d_subclass.jar to (ImageJ/Fiji install folder)/plugins.
2. (Re)start ImageJ/Fiji, then choose Plugns -> OpenMif in ImageJ main menu.

We confirmed that the plugin is successfully working with Fiji Life-line version ( 20130715, ImageJ 1.47v). The plugin may not work correctly using other version of Fiji/ImageJ.

About *.mif file format
The *.mif file is human readable text. 
Any text editor can be used for editing the contents. 
The file extension should be *.mif. 
Specify a property of image(s) in each line as matlab/java style:
(property1) = value; % for number values
(property2) = "value"; % for string values
(property2) = true; % for binary values
There are two types of properties: properties for image and properties for each channel.
Properties for image are:
	width: image width
	height: image height
	slices: number of images in a Z-stack
	frames: number of volumes (time lapse)
	subBG: whether subtract background (optional; default is false)
Properties for each channel are:
	path*: specify location (file path) for *.dcv, *.sif, and *.tif file. 
		Both relative path and absolute path can be used.
		If a time lapse image set was recorded as a multiple files (split), 
		you can specify multiple files separating with comma.
	flipX*: Whether flip the image along with the X axis (optional, default is false)
	flipY*: Whether flip the image along with the Y axis (optional, default is false)
	rotXY*: The angle for rotating the image in the XY plane (optional, default is 0)
	transX*: The length (pixels) of moving the image with the X axis (optional, default is 0)
	transY*: The length (pixels) of moving the image with the X axis (optional, default is 0)
	numC*: The number of color channel of the image file (optional, default is 0)
	useC*: The index of color channel used from the image file (optional, default is 0)
	padZ*: The number of padded slices for z-stack of the image file.
		The pixel values of slices from 1st to the specified will be set as 0, 
		in order to hide motion artifacts by Z piezo. (optional, default is 0)


