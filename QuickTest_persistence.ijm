inputPathToTimelapse = "C:/RESEARCH/Sholto_Conf/4.Thresholded/100-10/LME.tif";
outputPath = "C:/RESEARCH/Sholto_Conf/5.MEL_Output/100-10/";

inputPathToTimelapse = "/Volumes/Extreme SSD/RESEARCH/Vanderbilt Samples/2023/Marina/Output/Untreated/Thresh/MH15.2_hNPC_notrt_001.tif"
outputPath = "/Volumes/Extreme SSD/RESEARCH/Vanderbilt Samples/2023/Marina/Output/Untreated/MEL/"

if(!File.exists(outputPath))
	File.makeDirectory(outputPath);

SF = 1;
min_structure_volume = 5;

var stack = true;
run("Console");

close("*");
open(inputPathToTimelapse);
rename("Timelapse");
getDimensions(width, height, channels, slices, frames);
if(channels > 1){
	print("MEL only works with single channel images that were thresholded");
	return 0;
}
	
print(frames + " frames detected");
print(slices + " slices detected");

//f = 1;
// Note Fiji starts counting at 1
for(f = 1; f <= (frames-2); f++)
{		
	fText = "" + f + ".tif";
	if(f < 10)
		fText = "0" + fText;

	// extract fames from timelapse
	run("Make Substack...", "slices=1-" + slices +" frames=" + f);
	rename("Frame1");
	
	selectWindow("Timelapse");
	run("Make Substack...", "slices=1-" + slices +"  frames=" + (f+1));
	rename("Frame2");
	
	selectWindow("Timelapse");
	run("Make Substack...", "slices=1-" + slices +"  frames=" + (f+2));
	rename("Frame3");

	
	
	run("MEL Process", "frame_1_title=Frame1 frame_2_title=Frame2 frame_3_title=Frame3 event_persistence=true "+
	"min_structure_volume="+(min_structure_volume*SF)+" min_overlap_percentage=0.5 skeleton_distance_threshold="+(20*SF)+" "+
	"depolarisation_range_threshold="+(50*SF)+" depolarisation_structure_similarity_threshold=2.0"+
	" remove_duplicates=true duplicate_range="+(10*SF)+" debug_output=false " +
	"save_event_location_and_stats=true path_to_event_csv="+outputPath+"EventLocations"+f+".csv " + 
	"path_to_f1_stats_csv="+outputPath+"f1_stats"+f+".csv path_to_f2_stats_csv="+outputPath+"f2_stats"+f+".csv");
	
		
	
	selectWindow("Labels Frame1");
	save(outputPath + fText+"_L1");
	run("mpl-viridis");
	
	selectWindow("Labels Frame2");
	save(outputPath + fText+"_L2");
	run("mpl-viridis");
	
	//selectWindow("Matched graph");
	//run("mpl-viridis");
	
//	selectWindow("Labeled Skeleton F1 - small structures might be missing");
//	run("mpl-viridis");
//	
//	selectWindow("Labeled Skeleton F2 - small structures might be missing");
//	run("mpl-viridis");
	
//	selectWindow("Matched graph Fusion");
//	run("mpl-viridis");
//	
//	selectWindow("Matched graph Fission");
//	run("mpl-viridis");
	
	
	selectWindow("Fusion events");
	
	operation = "Erode"; // "Erode" or "Dilate" - this differs based on your default fiji settings
	if(stack)
	{	
		//run("Duplicate...", "duplicate");
		run("8-bit");
		run("Convert to Mask", "method=Otsu background=Dark calculate black");
		
		run(operation, "stack");
		run(operation, "stack");
		run(operation, "stack");	
		run("Smooth", "stack");
		run("Smooth", "stack");
		run("Smooth", "stack");
		run("Smooth", "stack");

		// run("Z Project...", "projection=[Max Intensity]");
		
		selectWindow("Fission events");
		//run("Duplicate...", "duplicate");
		run("8-bit");
		run("Convert to Mask", "method=Otsu background=Dark calculate black");
		
		run(operation, "stack");
		run(operation, "stack");
		run(operation, "stack");	
		run("Smooth", "stack");
		run("Smooth", "stack");
		run("Smooth", "stack");
		run("Smooth", "stack");

		// run("Z Project...", "projection=[Max Intensity]");
		
		selectWindow("Depolarisation events");
		//run("Duplicate...", "duplicate");
		run("8-bit");
		run("Convert to Mask", "method=Otsu background=Dark calculate black");
		
		run(operation, "stack");
		run(operation, "stack");
		run(operation, "stack");	
		run("Smooth", "stack");
		run("Smooth", "stack");
		run("Smooth", "stack");
		run("Smooth", "stack");

		// run("Z Project...", "projection=[Max Intensity]");
		
		//run("Merge Channels...", "c1=[MAX_Fission events-1] c2=[MAX_Fusion events-1] c3=[MAX_Depolarisation events-1] keep");
		//run("Merge Channels...", "c1=[MAX_Fission events] c2=[MAX_Fusion events] c3=[MAX_Depolarisation events]");
		run("Merge Channels...", "c1=[Fission events] c2=[Fusion events] c3=[Depolarisation events]");

		// save Events RGB image
		save(outputPath + "Events" + fText);
		run("Calculator Plus", "i1=Frame1 i2=Frame1 operation=[Scale: i2 = i1 x k1 + k2] k1=0.5 k2=0 create");
		run("RGB Color");
		run("Calculator Plus", "i1=Result i2=RGB operation=[Add: i2 = (i1+i2) x k1 + k2] k1=1 k2=0 create");
		rename("Output Time " + fText);
		
		// save Events overlaid on Thresholded Frame 1
		save(outputPath + fText);
		
		close("Result");
	}
	else {
		
		//run("Duplicate...", "duplicate");
		run("8-bit");
		run("Convert to Mask", "method=Otsu background=Dark calculate black");
		
		run(operation, "stack");
		run(operation, "stack");
		run(operation, "stack");
		run("Smooth", "stack");
		run("Smooth", "stack");
		
		selectWindow("Fission events");
		//run("Duplicate...", "duplicate");
		run("8-bit");
		run("Convert to Mask", "method=Otsu background=Dark calculate black");

		run(operation, "stack");
		run(operation, "stack");
		run(operation, "stack");
		run("Smooth", "stack");
		run("Smooth", "stack");
		
		selectWindow("Depolarisation events");
		//run("Duplicate...", "duplicate");
		run("8-bit");
		run("Convert to Mask", "method=Otsu background=Dark calculate black");

		run(operation, "stack");
		run(operation, "stack");
		run(operation, "stack");
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
	
	// run("Tile");
	// run("Synchronize Windows");

	selectWindow("Timelapse");
	close("\\Others");
}

selectWindow("MEL Results"); 
saveAs("Results", outputPath + "MEL_Results.csv");





