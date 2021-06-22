var stack = true;

//var extenstion = ".tif";
//open("C:/Users/rptheart/Dropbox/Research/MitoMorph/MEL2/Temp/Test/Frame1_Thresholded" + extenstion);
//rename("Frame1");
//open("C:/Users/rptheart/Dropbox/Research/MitoMorph/MEL2/Temp/Test/Frame2_Thresholded" + extenstion);
//rename("Frame2");

open("C:/RESEARCH/MEL/Thresholded/Con001_otsu.tif")
rename("Timelapse");

// extract fames from timelapse
run("Make Substack...", "slices=1-4 frames=1");
rename("Frame1");

selectImage("Timelapse");
run("Make Substack...", "slices=1-4 frames=2");
rename("Frame2");

run("Console");
run("MEL Process", "frame_1_title=Frame1 frame_2_title=Frame2  min_structure_volume=5 min_overlap_percentage=0.1 skeleton_distance_threshold=20 depolarisation_range_threshold=50 depolarisation_volume_similarity_threshold=0.2 debug_output=false"); 





selectWindow("Labels Frame1");
run("mpl-viridis");

selectWindow("Labels Frame2");
run("mpl-viridis");

//selectWindow("Matched graph");
//run("mpl-viridis");

selectWindow("Labeled Skeleton F1 - small structures might be missing");
run("mpl-viridis");

selectWindow("Labeled Skeleton F2 - small structures might be missing");
run("mpl-viridis");

selectWindow("Matched graph Fusion");
run("mpl-viridis");

selectWindow("Matched graph Fission");
run("mpl-viridis");


selectWindow("Fusion events");

if(stack)
{
	//run("Duplicate...", "duplicate");
	run("8-bit");
	run("Convert to Mask", "method=Otsu background=Dark calculate black");
	run("Dilate", "stack");
	run("Dilate", "stack");
	run("Smooth", "stack");
	run("Smooth", "stack");
	// run("Z Project...", "projection=[Max Intensity]");
	
	selectWindow("Fission events");
	//run("Duplicate...", "duplicate");
	run("8-bit");
	run("Convert to Mask", "method=Otsu background=Dark calculate black");
	run("Dilate", "stack");
	run("Dilate", "stack");
	run("Smooth", "stack");
	run("Smooth", "stack");
	// run("Z Project...", "projection=[Max Intensity]");
	
	selectWindow("Depolarisation events");
	//run("Duplicate...", "duplicate");
	run("8-bit");
	run("Convert to Mask", "method=Otsu background=Dark calculate black");
	run("Dilate", "stack");
	run("Dilate", "stack");
	run("Smooth", "stack");
	run("Smooth", "stack");
	// run("Z Project...", "projection=[Max Intensity]");
	
	//run("Merge Channels...", "c1=[MAX_Fission events-1] c2=[MAX_Fusion events-1] c3=[MAX_Depolarisation events-1] keep");
	//run("Merge Channels...", "c1=[MAX_Fission events] c2=[MAX_Fusion events] c3=[MAX_Depolarisation events]");
	run("Merge Channels...", "c1=[Fission events] c2=[Fusion events] c3=[Depolarisation events]");

	
	run("Calculator Plus", "i1=Frame1 i2=Frame1 operation=[Scale: i2 = i1 x k1 + k2] k1=0.5 k2=0 create");
	run("RGB Color");
	run("Calculator Plus", "i1=Result i2=RGB operation=[Add: i2 = (i1+i2) x k1 + k2] k1=1 k2=0 create");
	rename("Output Time " + 0)
	close("Result")
}
else {
	
	//run("Duplicate...", "duplicate");
	run("8-bit");
	run("Convert to Mask", "method=Otsu background=Dark calculate black");
	run("Dilate", "stack");
	run("Dilate", "stack");
	run("Smooth", "stack");
	run("Smooth", "stack");
	
	selectWindow("Fission events");
	//run("Duplicate...", "duplicate");
	run("8-bit");
	run("Convert to Mask", "method=Otsu background=Dark calculate black");
	run("Dilate", "stack");
	run("Dilate", "stack");
	run("Smooth", "stack");
	run("Smooth", "stack");
	
	selectWindow("Depolarisation events");
	//run("Duplicate...", "duplicate");
	run("8-bit");
	run("Convert to Mask", "method=Otsu background=Dark calculate black");
	run("Dilate", "stack");
	run("Dilate", "stack");
	run("Smooth", "stack");
	run("Smooth", "stack");
	
	//run("Merge Channels...", "c1=[MAX_Fission events-1] c2=[MAX_Fusion events-1] c3=[MAX_Depolarisation events-1] keep");
	run("Merge Channels...", "c1=[Fission events] c2=[Fusion events] c3=[Depolarisation events]");

	// REMOVE THIS LINE
//selectWindow("Frame1_Thresholded.tif");
//run("Add Image...", "image=RGB x=0 y=0 opacity=50");

selectWindow("Frame1");
run("Duplicate...", " ");
run("Add Image...", "image=RGB x=0 y=0 opacity=50");
}




run("Tile");
run("Synchronize Windows");


