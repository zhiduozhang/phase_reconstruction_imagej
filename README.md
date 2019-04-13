# DHM Reconstruction Plugin for FIJI

How to install:

1. Follow the instructions at the start of https://imagej.net/Installing_3rd_party_plugins by adding the following update sites:
https://sites.imagej.net/DHMReconstruction/
https://sites.imagej.net/IJ-OpenCV/

NOTE: IJ-OpenCV must appear after DHMReconstruction on the update list, otherwise some libraries will not download correctly. 

2. Open {fiji directory}\Fiji.app\jars and delete ejml-0.24.jar


To Do List:

* Save settings
* Test for image stacks
* Autofocus
* Deleting ejml-0.24 seems to break Fiji updater. Need to resolve dependency conflict.
* Redesign GUI
