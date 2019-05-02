# DHM Reconstruction Plugin for FIJI

### How to install:

Follow the instructions at the start of https://imagej.net/Installing_3rd_party_plugins by adding the following update sites:
https://sites.imagej.net/DHMReconstruction/ 

### How to use

Once installed, the plugin will appear near the bottom of the plugins dropdown. The settings small constraint and large constraint define how big the zero,real,imaginary orders can be during thresholding process. Small constraint removes all blobs where the area is smaller than (width/small_constraint)*(height/small_constraint). Similarly, large constraint will ensure all blobs are smaller than the same definition. 

### To Do List:

* Autofocus
* Add errors for invalid values
* Automatically convert RGB images to 16Bit images
* Improve runtime