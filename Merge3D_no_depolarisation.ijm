// ONLY MODIFY THESE TWO LINES
fileName = "Con_003.2.tif";
rootPath = "/Users/rensu/Documents/Temp/Sholto Remove Depolarisation/"


intensity_TimeLapse = rootPath + fileName;
eventImage_timelapse = rootPath + "Events_" + fileName + "/";
output_path = rootPath + "Output" + fileName + "/";


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
	
	// THIS SECTION REMOVES THE DEPOLARISATION EVENTS
	run("Split Channels");
	run("Merge Channels...", "c1=[Events (red)] c2=[Events (green)]");
	rename("Events");
	//
	
	if(isOpen(intensity_TimeLapse + " - T=" + i))
	{
		print("Selecting: " + intensity_TimeLapse + " - T=" + i);
		selectImage(intensity_TimeLapse + " - T=" + i);
	}
	else(isOpen(fileName + " - T=" + i))
	{
		print("Selecting: " + fileName + " - T=" + i);
		selectImage(fileName + " - T=" + i);
	}
		
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
if(isOpen(intensity_TimeLapse + " - T=" + i))
{
	close(intensity_TimeLapse + " - T=" + i);
}
else(isOpen(fileName + " - T=" + i))
{
	close(fileName + " - T=" + i);
}
		

run("Concatenate...", "all_open open");
saveAs("Tiff", output_path + fileName.substring(0, fileName.lastIndexOf(".")) + ".tif");
run("Enhance Contrast...", "saturated=0.35 process_all use");
saveAs("Tiff", output_path + fileName.substring(0, fileName.lastIndexOf(".")) + "_enhanced.tif");
close("*");
