fileName = "Con_002.tif";

intensity_TimeLapse = "C:/Users/rptheart/Downloads/MEL_Quick_Test/" + fileName;
eventImage_timelapse = "C:/Users/rptheart/Downloads/MEL_Quick_Test/Events_" + fileName + "/";
output_path = "C:/Users/rptheart/Downloads/MEL_Quick_Test/Output" + fileName + "/";

if(!File.exists(output_path))
	File.makeDirectory(output_path);

//eventList = getFileList(eventImage_timelapse);

close("*");

	
run("Bio-Formats", "open=[" + intensity_TimeLapse + "] color_mode=Default rois_import=[ROI manager] split_timepoints view=Hyperstack stack_order=XYCZT");
numFrames = nImages;
for(i = 0; i < numFrames - 1; i++)
{
	frameText = "" + (i+1);
	if(i < 10)
		frameText = "0" + frameText;
	print("Processing: " + frameText);


	open(eventImage_timelapse + "Events" + frameText + ".tif");
	run("Gaussian Blur 3D...", "x=0.2 y=0.2 z=1");
	rename("Events");
	
	selectImage(intensity_TimeLapse + " - T=" + i);
	rename("Intensity");

	numSlices = nSlices;
	for (s = 1; s <= numSlices; s++)
	{
		selectImage("Events");
		setSlice(s);
		selectImage("Intensity");
		setSlice(s);
		run("Add Image...", "image=Events x=0 y=0 opacity=50 zero");
	}
	selectImage("Intensity");
	run("Flatten", "stack");
	
	rename(frameText);
	run("RGB Color");

	wait(10);
		
//	selectImage("Intensity");
//	close();	
	selectImage("Events");
	close();	
}
selectImage(intensity_TimeLapse + " - T=" + i);
close();	

run("Concatenate...", "all_open open");
saveAs("Tiff", output_path + fileName.substring(0, fileName.lastIndexOf(".")) + ".tif");
run("Enhance Contrast...", "saturated=0.35 process_all use");
saveAs("Tiff", output_path + fileName.substring(0, fileName.lastIndexOf(".")) + "_enhanced.tif");
close("*");
