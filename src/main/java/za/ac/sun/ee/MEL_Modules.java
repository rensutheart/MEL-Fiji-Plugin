/*
 * To the extent possible under law, the ImageJ developers have waived
 * all copyright and related or neighboring rights to this tutorial code.
 *
 * See the CC0 1.0 Universal license for details:
 *     http://creativecommons.org/publicdomain/zero/1.0/
 */

package za.ac.sun.ee;

import net.imagej.Dataset;
import net.imagej.ImageJ;
import net.imagej.ops.OpService;
import net.imglib2.type.numeric.RealType;

import org.jgrapht.Graph;
import org.jgrapht.graph.DefaultEdge;
import org.jgrapht.graph.DefaultUndirectedGraph;
import org.jgrapht.traverse.DepthFirstIterator;
import org.scijava.ItemIO;
import org.scijava.command.Command;
import org.scijava.plugin.Attr;
import org.scijava.plugin.Parameter;
import org.scijava.plugin.Plugin;
//import org.scijava.table.DefaultGenericTable;
//import org.scijava.table.DoubleColumn;
//import org.scijava.table.GenericTable;
//import org.scijava.table.IntColumn;
//import org.scijava.table.Table;
import org.scijava.ui.UIService;

import ij.IJ;
import ij.ImagePlus;
import ij.WindowManager;
import ij.gui.NewImage;
import ij.measure.ResultsTable;
import mcib3d.geom.Object3DVoxels;
import mcib3d.geom.Point3D;
import mcib3d.geom.Vector3D;
import mcib3d.geom.Voxel3D;
import mcib3d.image3d.ImageByte;
import mcib3d.image3d.ImageHandler;
import mcib3d.image3d.ImageInt;
import mcib3d.image3d.ImageLabeller;

import java.io.File;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.Set;

//, attrs = {@Attr(name = "no-legacy") }
@Plugin(type = Command.class, name = "MEL", description = "Automatically calculate the mitochondrial fission, fusion and depolarisation event locations ", menuPath = "Plugins>MEL Process", headless = false)
public class MEL_Modules<T extends RealType<T>> implements Command {
//	@Parameter
//	private Dataset currentData;
//
//	@Parameter
//	private UIService uiService;
//
//	@Parameter
//	private OpService opService;

	// the "required = false" part is so that I can run it automatically through a
	// macro using the run() function without seeing a dialog box
	@Parameter(label = "The window title of Frame 1", required = false)
	private String frame_1_title = "";

	@Parameter(label = "The window title of Frame 2", required = false)
	private String frame_2_title = "";

	@Parameter(label = "The minimum volume of the thresholded structures (voxels)", required = false)
	private int min_structure_volume = 5; // TODO: This can be dependent on Voxel size (see Mitochondrial Analyzer macro)

	@Parameter(label = "The minimum percentage of volume overlap to be considered a match (as 0-1)", required = false)
	private float min_overlap_percentage = 0.1f;

	@Parameter(label = "The distance threholds between the skeletons in two frames to be considered a match (pixels)", required = false)
	private float skeleton_distance_threshold = 20;

	@Parameter(label = "The distance a structure can move before being considered as depolarised (pixels)", required = false)
	private float depolarisation_range_threshold = 50;

	@Parameter(label = "For a structure that moved, what volume similarity must it have to be considered a match (as 0-1)", required = false)
	private float depolarisation_volume_similarity_threshold = 0.2f; // this is a percentage: 0 means that they must be exactly the same, 0.2 means
																		// that the other structure may be 20% larger or smaller
	@Parameter(label = "Display full debug output in the console", persist = false, required = false)
	private boolean debug_output;

	@Override
	public void run() {
		long startTime = System.currentTimeMillis();
		// I'm assuming the input images are Pre-processed and thresholded.

		/*
		 * LOAD THRESHOLDED FRAMES
		 */
//		String[] imageTitles = WindowManager.getImageTitles();
//		for (String imgTitle : imageTitles) {
//			System.out.println(imgTitle);
//		}
//		String title_Frame1 = imageTitles[0];// "Frame1_Thresholded.tif";
//		String title_Frame2 = imageTitles[1];// "Frame2_Thresholded.tif";
//
//		ImagePlus imagePlus_Frame1 = WindowManager.getImage(title_Frame1);
//		ImagePlus imagePlus_Frame2 = WindowManager.getImage(title_Frame2);

		ImagePlus imagePlus_Frame1 = WindowManager.getImage(frame_1_title);
		ImagePlus imagePlus_Frame2 = WindowManager.getImage(frame_2_title);

		// convert from whatever bit depth and format to 8-bit for further processing
		imagePlus_Frame1.setProcessor(imagePlus_Frame1.getProcessor().convertToByteProcessor());
		imagePlus_Frame2.setProcessor(imagePlus_Frame2.getProcessor().convertToByteProcessor());

		/*
		 * LABEL MITOCHONDRIAL STRUCTURES
		 */
		// Get the segmented labeled mitochondrial structures
		ImageInt labels_F1 = getLabeledImage(imagePlus_Frame1, min_structure_volume);
		ImageInt labels_F2 = getLabeledImage(imagePlus_Frame2, min_structure_volume);

		// this 1 removes the background and only retains mitochondrial structures
		Object3DVoxels[] labelVoxels_F1 = labelImageTo3DVoxelArray(labels_F1, 1);
		Object3DVoxels[] labelVoxels_F2 = labelImageTo3DVoxelArray(labels_F2, 1);

		// Display results
		labels_F1.show("Labels Frame1");
		labels_F2.show("Labels Frame2");

		System.out.println("Number of Labels F1: " + labelVoxels_F1.length + " Labels Max: " + labels_F1.getMax());
		System.out.println("Number of Labels F2: " + labelVoxels_F2.length + " Labels Max: " + labels_F2.getMax());

		/*
		 * CALCULATE MEL PARAMTERS
		 */
		// Get the overlapping volumes - IMPORTANT - None of these include the
		// background, hence the index
		// is off by 1 compared to the label image label number
		int[][] overlappingVolumes = getOverlappingVolumes(labelVoxels_F1, labelVoxels_F2);

		Point3D[] centerOfStructures_F1 = getCenterOfStructures(labelVoxels_F1);
		Point3D[] centerOfStructures_F2 = getCenterOfStructures(labelVoxels_F2);

		int[] numVoxelsInStructures_F1 = getNumVoxelsInStructures(labelVoxels_F1);
		int[] numVoxelsInStructures_F2 = getNumVoxelsInStructures(labelVoxels_F2);

		List<List<Integer>> associatedLabelsBetweenFrames_F1toF2 = getAssociatedLabelsBetweenFrames(overlappingVolumes, centerOfStructures_F1, numVoxelsInStructures_F1, centerOfStructures_F2,
				numVoxelsInStructures_F2, depolarisation_range_threshold, depolarisation_volume_similarity_threshold);
		List<List<Integer>> associatedLabelsBetweenFrames_F2toF1 = getAssociatedLabelsBetweenFrames(transposeMatrix(overlappingVolumes), centerOfStructures_F2, numVoxelsInStructures_F2,
				centerOfStructures_F1, numVoxelsInStructures_F1, depolarisation_range_threshold, depolarisation_volume_similarity_threshold);

		// I don't seem to need this anymore for the new approach
//		List<List<Integer>> associatedLabelsWithinFrame_F1 = getAssociatedLabelsWithinFrame(associatedLabelsBetweenFrames_F1toF2, associatedLabelsBetweenFrames_F2toF1);
//		List<List<Integer>> associatedLabelsWithinFrame_F2 = getAssociatedLabelsWithinFrame(associatedLabelsBetweenFrames_F2toF1, associatedLabelsBetweenFrames_F1toF2);

		List<List<Float>> percentageOverlap_F1toF2 = getRelativePercentageOverlap(overlappingVolumes, associatedLabelsBetweenFrames_F1toF2);
		List<List<Float>> percentageOverlap_F2toF1 = getRelativePercentageOverlap(transposeMatrix(overlappingVolumes), associatedLabelsBetweenFrames_F2toF1);

		List<List<Integer>> reducedAssociatedLabelsBetweenFrames_F1toF2 = reduceAssociateLabels(associatedLabelsBetweenFrames_F1toF2, percentageOverlap_F1toF2, min_overlap_percentage);
		List<List<Integer>> reducedAssociatedLabelsBetweenFrames_F2toF1 = reduceAssociateLabels(associatedLabelsBetweenFrames_F2toF1, percentageOverlap_F2toF1, min_overlap_percentage);

		// Display Results
		if (debug_output)
			System.out.println("associatedLabelsBetweenFrames_F1toF2.size " + associatedLabelsBetweenFrames_F1toF2.size());
		if (debug_output)
			System.out.println("associatedLabelsBetweenFrames_F2toF1.size " + associatedLabelsBetweenFrames_F2toF1.size());

//		if(debugOutput) System.out.println("associatedLabelsWithinFrame_F1.size " + associatedLabelsWithinFrame_F1.size());
//		if(debugOutput) System.out.println("associatedLabelsWithinFrame_F2.size " + associatedLabelsWithinFrame_F2.size());

		if (debug_output)
			System.out.println("percentageOverlap_F1toF2.size " + percentageOverlap_F1toF2.size());
		if (debug_output)
			System.out.println("percentageOverlap_F2toF1.size " + percentageOverlap_F2toF1.size());

		if (debug_output)
			System.out.println("reducedAssociatedLabelsBetweenFrames_F1toF2.size " + reducedAssociatedLabelsBetweenFrames_F1toF2.size());
		if (debug_output)
			System.out.println("reducedAssociatedLabelsBetweenFrames_F2toF1.size " + reducedAssociatedLabelsBetweenFrames_F2toF1.size());

		// Viewer3D_Utils can perform marching cubes and .obj export

		/*
		 * SKELETONIZE
		 */
		// These skeletons does not include the background as an index, hence
		// does the label image, but the same compared to the MEL parameter Lists
		List<Object3DVoxels> labelsSkeletons_F1 = labelsToSkeleton(labels_F1, "Labeled Skeleton F1");
		List<Object3DVoxels> labelsSkeletons_F2 = labelsToSkeleton(labels_F2, "Labeled Skeleton F2");

		List<Graph<Vector3D, DefaultEdge>> labelsSkeletonGraphs_F1 = allSkeletonsToGraphs(labelsSkeletons_F1);
		List<Graph<Vector3D, DefaultEdge>> labelsSkeletonGraphs_F2 = allSkeletonsToGraphs(labelsSkeletons_F2);

		// Display Results
		if (debug_output)
			System.out.println("labelsSkeletons_F1.size " + labelsSkeletons_F1.size());
		if (debug_output)
			System.out.println("labelsSkeletons_F2.size " + labelsSkeletons_F2.size());

		if (debug_output)
			System.out.println("labelsSkeletonGraphs_F1.size " + labelsSkeletonGraphs_F1.size());
		if (debug_output)
			System.out.println("labelsSkeletonGraphs_F2.size " + labelsSkeletonGraphs_F2.size());

		/*
		 * DETECT FUSION
		 */
		// Find the matched skeleton F1 to F2
		if (debug_output)
			System.out.println("\nFIND FUSION");
		Graph<GraphNode, DefaultEdge> matchedGraphs_onF2 = matchGraphNodesBetweenFrames(labelsSkeletonGraphs_F1, labelsSkeletonGraphs_F2, reducedAssociatedLabelsBetweenFrames_F1toF2,
				skeleton_distance_threshold);
		showGraphNodesAsImage(matchedGraphs_onF2, labels_F1.sizeX, labels_F1.sizeY, labels_F1.sizeZ, true, "Matched graph Fusion");
//		showGraphNodesAsImage(matchedGraphs_onF2, labels_F1.sizeX, labels_F1.sizeY, labels_F1.sizeZ, false, "Unmatched graph");

		// Find all the events F1 to F2 (fusion)
		List<Vector3D> fusionEventLocations = findEvents(matchedGraphs_onF2, labelsSkeletonGraphs_F1, true, 10);
		ImageInt fusionEventsImage = eventsToImage(fusionEventLocations, labels_F1.sizeX, labels_F1.sizeY, labels_F1.sizeZ, "Fusion events");

		/*
		 * DETECT FISSION
		 */
		// Find the matched skeleton F2 to F1
		if (debug_output)
			System.out.println("\nFIND FISSION");
		Graph<GraphNode, DefaultEdge> matchedGraphs_onF1 = matchGraphNodesBetweenFrames(labelsSkeletonGraphs_F2, labelsSkeletonGraphs_F1, reducedAssociatedLabelsBetweenFrames_F2toF1,
				skeleton_distance_threshold);
		showGraphNodesAsImage(matchedGraphs_onF1, labels_F1.sizeX, labels_F1.sizeY, labels_F1.sizeZ, true, "Matched graph Fission");
//		showGraphNodesAsImage(matchedGraphs_onF1, labels_F1.sizeX, labels_F1.sizeY, labels_F1.sizeZ, false, "Unmatched graph");

		// Find all the events F2 to F1 (fission)
		List<Vector3D> fissionEventLocations = findEvents(matchedGraphs_onF1, labelsSkeletonGraphs_F2, true, 10);
		ImageInt fissionEventsImage = eventsToImage(fissionEventLocations, labels_F1.sizeX, labels_F1.sizeY, labels_F1.sizeZ, "Fission events");

		/*
		 * DETECT DEPOLARISATION
		 */
		if (debug_output)
			System.out.println("\nFIND DEPOLARISATION");
		// NOTE: I'm NOT using the reducedAssociatedLabelsBetweenFrames_F1toF2 list
		// here, since for depolarisation I would err on the side of caution, and if
		// there is a slight possibility that the even did join to another structure or
		// moved, then I don't want to mark it
		List<Vector3D> depolarisationEventLocations = findDepolarisationEvents(associatedLabelsBetweenFrames_F1toF2, labelVoxels_F1);
		ImageInt depolarisationEventsImage = eventsToImage(depolarisationEventLocations, labels_F1.sizeX, labels_F1.sizeY, labels_F1.sizeZ, "Depolarisation events");

		int fusionEventCount = fusionEventLocations.size();
		int fissionEventCount = fissionEventLocations.size();
		int depolarisationEventCount = depolarisationEventLocations.size();
		System.out.println(String.format("Number of events:\n\tFusion = %d\n\tFission = %d\n\tDepolarisation = %d", fusionEventCount, fissionEventCount, depolarisationEventCount));
		System.out.println("Fission:Fusion ratio = " + ((float) fissionEventCount / (float) fusionEventCount));

		ResultsTable table = ResultsTable.getResultsTable();
		table.incrementCounter();

		table.addValue("Fusion events", fusionEventCount);
		table.addValue("Fission events", fissionEventCount);
		table.addValue("Depolarisation events", depolarisationEventCount);
		// table.addColumns();

		table.show("MEL Results");

		/*
		 * // TEST CODE: Note the - 1, that is since background is removed, now label 1
		 * is // at index 0. int label_F1_from = 40; int label_F1_to =
		 * associatedLabelsWithinFrame_F1.get(label_F1_from - 1).get(0) + 1;
		 * 
		 * if(debugOutput) System.out.println("Structures associated within with label "
		 * + label_F1_from); for (int label_to :
		 * associatedLabelsWithinFrame_F1.get(label_F1_from - 1)) { if(debugOutput)
		 * System.out.println(label_to + 1); }
		 * 
		 * if(debugOutput) System.out.println("Finding 5 closest points between label "
		 * + label_F1_from + " and " + label_F1_to); // IMPORTANT NOTE: labelsSkeletons
		 * matches the Label image and therefore // includes the background, hence no -1
		 * for label number List<VectorPair> closestVectorPair =
		 * getClosestPointsBetweenStructures(labelsSkeletons_F1.get(label_F1_from),
		 * labelsSkeletons_F1.get(label_F1_to), 5);
		 */

//		Img<UnsignedIntType> skeletonImg_Frame1 = ImageJFunctions.wrap(WindowManager.getImage(title_Frame1));
//		Img<UnsignedIntType> skeletonImg_Frame2 = ImageJFunctions.wrap(WindowManager.getImage(title_Frame2));

//		IJ.run(skeleton_F1, "Analyze Skeleton (2D/3D)", "prune=[none] calculate show");
		// Capture results

		// TODO: Find a way to match the skeleton to the labeled structure:
		// Option 1: Calculate the skeletonization many times (once for each structure)
		// Option 2: Sample the value of a single skeleton point to determine the label
		// value
		// Option 3: multiply binarized skeleton with label image, and then extract
		// labels of skeleton?

//		WindowManager.closeAllWindows();

		long endTime = System.currentTimeMillis();
		System.out.println("MEL - Total execution time: " + (endTime - startTime) + "ms");
	}

	public ImageInt getLabeledImage(ImagePlus imgPlus, int minStructureVolume) {
		// Create Labeler
		ImageLabeller labeler = new ImageLabeller();
		labeler.setMinSize(minStructureVolume);
		// labeler.setMaxsize(1048576);

		// Load binarized image in the ImageHandler format
		ImageHandler img = ImageHandler.wrap(imgPlus);
		ImageInt bin = img.thresholdAboveInclusive(128); // not necessary for thresholded images, but for safety

		// Get and return labled image
		return labeler.getLabels(bin);
	}

	public int[][] getOverlappingVolumes(Object3DVoxels[] labelVoxels_F1, Object3DVoxels[] labelVoxels_F2) {
		long startTime = System.currentTimeMillis();
		System.out.println("Started: getOverlappingVolumes()");

		int numLabels_F1 = labelVoxels_F1.length;
		int numLabels_F2 = labelVoxels_F2.length;

		int[][] overlappingVolumes = new int[numLabels_F1][numLabels_F2];

		// NOTES: This excludes background if labelImageTo3DVoxelArray() was called with
		// a start label of 1
		for (int label_F1 = 0; label_F1 < numLabels_F1; ++label_F1) {
			for (int label_F2 = 0; label_F2 < numLabels_F2; ++label_F2) {
				// TODO: This line of code is very slow - consider multiplying a single label
				// structure with the entire image in the other frame
				Object3DVoxels intersection = labelVoxels_F1[label_F1].getIntersectionObject(labelVoxels_F2[label_F2]);

				int volumeOverlap = 0;
				if (intersection != null) // structures are not disjoint
					volumeOverlap = intersection.getVoxels().size();

				overlappingVolumes[label_F1][label_F2] = volumeOverlap;

				// This assumes a startLabel of 1
				if (volumeOverlap != 0) {
					if (debug_output)
						System.out.println("Label F1 " + (label_F1 + 1) + "  Label F2 " + (label_F2 + 1) + "  Volume Overlap " + volumeOverlap);
				}
			}
		}

		long endTime = System.currentTimeMillis();
		System.out.println("getOverlappingVolumes() - Total execution time: " + (endTime - startTime) + "ms");

		return overlappingVolumes;
	}

	public Point3D[] getCenterOfStructures(Object3DVoxels[] labelVoxels) {
		Point3D[] centerOfStructures = new Point3D[labelVoxels.length];

		for (int i = 0; i < centerOfStructures.length; i++) {
			centerOfStructures[i] = labelVoxels[i].getCenterAsPoint();
		}

		return centerOfStructures;
	}

	public int[] getNumVoxelsInStructures(Object3DVoxels[] labelVoxels) {
		int[] numVoxelsInStructure = new int[labelVoxels.length];

		for (int i = 0; i < numVoxelsInStructure.length; i++) {
			numVoxelsInStructure[i] = labelVoxels[i].getVoxels().size();
		}

		return numVoxelsInStructure;
	}

	// The startLabel is primarily used for removing background (therefore
	// startLabel = 1)
	public Object3DVoxels[] labelImageTo3DVoxelArray(ImageInt labeledImage, int startLabel) {
		int maxLabel = (int) labeledImage.getMax();
		Object3DVoxels[] outputArray = new Object3DVoxels[maxLabel - startLabel + 1];

		for (int i = startLabel; i <= maxLabel; ++i) {
			outputArray[i - startLabel] = new Object3DVoxels(labeledImage, i);
		}

		return outputArray;
	}

	public Object3DVoxels[] labelImageTo3DVoxelArray(ImageInt labeledImage) {
		return labelImageTo3DVoxelArray(labeledImage, 0);
	}

	public List<List<Integer>> getAssociatedLabelsBetweenFrames(int[][] overlappingVolumes, Point3D[] centerOfStructures_F1, int[] numVoxelsInStructures_F1, Point3D[] centerOfStructures_F2,
			int[] numVoxelsInStructures_F2, float depolarisationRange, float depolarisationVolumeSimilarity) {
		long startTime = System.currentTimeMillis();

		List<List<Integer>> associatedLabelsBetweenFrames = new ArrayList<List<Integer>>(overlappingVolumes.length);

		// NOTES: This excludes background if labelImageTo3DVoxelArray() was called with
		// a start label of 1
		for (int i = 0; i < overlappingVolumes.length; ++i) {
			List<Integer> tempList = new ArrayList<Integer>();
			for (int j = 0; j < overlappingVolumes[0].length; ++j) {
				if (overlappingVolumes[i][j] > 0) {
					tempList.add(j);
					// This assumes a startLabel of 1
					if (debug_output)
						System.out.println("Associated between " + (i + 1) + " and " + (j + 1));
				}
			}

			/*
			 * This section tries to find "Close" structures of similar volume, before
			 * allowing something to be considered as depolarisation
			 */
			// there were no associated labels (depolarized? from Frame 1 to Frame 2, or
			// appeared from Frame 1 to Frame 2 when this function is called with the
			// transpose)
			if (tempList.size() == 0) {
				Point3D centerCurrent = centerOfStructures_F1[i];
				int numVoxelsCurrent = numVoxelsInStructures_F1[i];

				double minDistance = Float.POSITIVE_INFINITY;
				double minVoxelRange = Float.POSITIVE_INFINITY;
				int F2_index = Integer.MAX_VALUE;

				// loop through the structures in the other frame and try to find the closest
				// structure close to the current center that is also similar volume
				for (int j = 0; j < centerOfStructures_F2.length; j++) {
					double newDistance = centerCurrent.distance(centerOfStructures_F2[j]);
					if (newDistance < depolarisationRange) {
						double newVoxelRange = Math.abs(numVoxelsCurrent / numVoxelsInStructures_F2[j] - 1);
						// check if they are of similar volume, could be larger or smaller, looking for
						// percentage difference
						if (newVoxelRange < depolarisationVolumeSimilarity) {
							// Prioritise the distance, if closer, and still within the voxel num range,
							// then set as new one.
							if (newDistance <= minDistance) { // && newVoxelRange < minVoxelRange
								if (debug_output)
									System.out.println("NOT DEPOLARISATION. Link between " + (i + 1) + " and " + (j + 1) + " with distance " + newDistance + " and voxel similarity " + newVoxelRange);
								minDistance = newDistance;
								minVoxelRange = newVoxelRange;
								F2_index = j;
							}
						}
					}
				}

				// Only add, if a match could be found
				if (F2_index != Integer.MAX_VALUE) {
					tempList.add(F2_index);

					if (debug_output)
						System.out.println("... SAVED Associated between " + (i + 1) + " and " + (F2_index + 1));
				}
			}

			associatedLabelsBetweenFrames.add(tempList);
		}

		long endTime = System.currentTimeMillis();
		System.out.println("getAssociatedLabelsBetweenFrames() - Total execution time: " + (endTime - startTime) + "ms");

		return associatedLabelsBetweenFrames;
	}

	// https://stackoverflow.com/questions/26197466/transposing-a-matrix-from-a-2d-array
	public int[][] transposeMatrix(int[][] matrix) {
		long startTime = System.currentTimeMillis();

		int m = matrix.length;
		int n = matrix[0].length;

		int[][] transposedMatrix = new int[n][m];

		for (int x = 0; x < n; x++) {
			for (int y = 0; y < m; y++) {
				transposedMatrix[x][y] = matrix[y][x];
			}
		}

		long endTime = System.currentTimeMillis();
		System.out.println("transposeMatrix() - Total execution time: " + (endTime - startTime) + "ms");

		return transposedMatrix;
	}

	// This is also known as Back and Forth Structure Matching
	public List<List<Integer>> getAssociatedLabelsWithinFrame(List<List<Integer>> associatedLabelsBetweenFrames_F1, List<List<Integer>> associatedLabelsBetweenFrames_F2) {
		long startTime = System.currentTimeMillis();

		List<List<Integer>> associatedLabelsWithinFrame = new ArrayList<List<Integer>>(associatedLabelsBetweenFrames_F1.size());

		// Reference counter
		int labelIndex_F1 = 0;
		for (List<Integer> list_F1 : associatedLabelsBetweenFrames_F1) {
			// This assumes a startLabel of 1
			if (debug_output)
				System.out.println("Frame 1 Label: " + (labelIndex_F1 + 1));
			Set<Integer> tempSet = new HashSet<Integer>();

			if (list_F1.size() > 0) {
				for (int labelNum_F2 : list_F1) {
					if (debug_output)
						System.out.println("\t Frame 2 Label: " + (labelNum_F2 + 1));
					for (int labelNum_F1 : associatedLabelsBetweenFrames_F2.get(labelNum_F2)) {
						if (labelIndex_F1 != labelNum_F1) // ensure the label is not associated to itself
						{
							tempSet.add(labelNum_F1);
							// This assumes a startLabel of 1
							if (debug_output)
								System.out.println("\t\t Label " + (labelIndex_F1 + 1) + " is associated with label " + (labelNum_F1 + 1));
						}
					}
				}
			} else {
				if (debug_output)
					System.out.println("\t Label " + (labelIndex_F1 + 1) + " has no associated structures in Frame 2");
			}

			List<Integer> tempList = new ArrayList<Integer>();
			for (Integer label : tempSet) {
				tempList.add(label);
			}
			associatedLabelsWithinFrame.add(tempList);
			labelIndex_F1++;
		}

		long endTime = System.currentTimeMillis();
		System.out.println("getAssociatedLabelsWithinFrame() - Total execution time: " + (endTime - startTime) + "ms");

		return associatedLabelsWithinFrame;
	}

	public List<List<Float>> getRelativePercentageOverlap(int[][] overlappingVolumes, List<List<Integer>> associatedLabelsBetweenFrames_F1) {
		long startTime = System.currentTimeMillis();

		List<List<Float>> percentageOverlap = new ArrayList<List<Float>>(associatedLabelsBetweenFrames_F1.size());

		for (int labelNum_F1 = 0; labelNum_F1 < associatedLabelsBetweenFrames_F1.size(); ++labelNum_F1) {
			List<Float> labelVolumeOverlaps_F1 = new ArrayList<Float>();
			for (int labelNum_F2 : associatedLabelsBetweenFrames_F1.get(labelNum_F1)) {
				labelVolumeOverlaps_F1.add((float) overlappingVolumes[labelNum_F1][labelNum_F2]);
			}

			float sum = 0;
			for (float d : labelVolumeOverlaps_F1)
				sum += d;

			for (int i = 0; i < labelVolumeOverlaps_F1.size(); ++i) {
				labelVolumeOverlaps_F1.set(i, labelVolumeOverlaps_F1.get(i) / sum);

				// System.out.println("Percentage overlap: " + labelVolumeOverlaps_F1.get(i) + "
				// (" + labelNum_F1 + ", " + i + ")");
			}

			percentageOverlap.add(labelVolumeOverlaps_F1);
		}

		long endTime = System.currentTimeMillis();
		System.out.println("getRelativePercentageOverlap() - Total execution time: " + (endTime - startTime) + "ms");

		return percentageOverlap;
	}

	private List<List<Integer>> reduceAssociateLabels(List<List<Integer>> associatedLabelsBetweenFrames_F1toF2, List<List<Float>> percentageOverlap_F1toF2, float cutoffPercentage) {
		long startTime = System.currentTimeMillis();

		List<List<Integer>> reducedAssociatedLabelsBetweenFrames = new ArrayList<List<Integer>>();

		// NOTES: This excludes background if labelImageTo3DVoxelArray() was called with
		// a start label of 1
		for (int i = 0; i < associatedLabelsBetweenFrames_F1toF2.size(); ++i) {
			List<Integer> tempList = new ArrayList<Integer>();
			for (int j = 0; j < associatedLabelsBetweenFrames_F1toF2.get(i).size(); ++j) {
				if (percentageOverlap_F1toF2.get(i).get(j) > cutoffPercentage)
					tempList.add(associatedLabelsBetweenFrames_F1toF2.get(i).get(j));
				else if (debug_output)
					System.out.println("Removed match between label " + (i + 1) + " and label " + (j + 1) + " in the other frame with percentage " + percentageOverlap_F1toF2.get(i).get(j));
			}
			reducedAssociatedLabelsBetweenFrames.add(tempList);
		}

		long endTime = System.currentTimeMillis();
		System.out.println("reduceAssociateLabels() - Total execution time: " + (endTime - startTime) + "ms");

		return reducedAssociatedLabelsBetweenFrames;
	}

	// Returns the skeleton voxels for each label structure
	public List<Object3DVoxels> labelsToSkeleton(ImageInt labeledImage, String title) {
		long startTime = System.currentTimeMillis();

//		Object3DVoxels[] labelVoxels = labelImageTo3DVoxelArray(labeledImage);

		List<Object3DVoxels> labelsSkeletons;// = new ArrayList<Object3DVoxels>((int) labeledImage.getMax());
		// initialize ArrayList
//		for (int i = 0; i < labelVoxels.length; ++i) {
//			labelsSkeletons.add(new Object3DVoxels());
//		}

		// Use the labeled image to ensure small structures are removed already in
		// previous step (minStructureSize)
		ImageByte binarizedLabel = labeledImage.thresholdAboveInclusive(1);
		ImagePlus skeleton = binarizedLabel.getImagePlus();

//		// upscaled times 2 since normal skeletonization tends to entirely remove small structures
//		IJ.run(skeleton, "Size...", "width=" + skeleton.getWidth() * 2 + " height=" + skeleton.getHeight() * 2 + " depth=" + skeleton.getNSlices() + " constrain average interpolation=Bicubic");

		// Note that these macros change the image itself so if using the original
		// images, use .duplicate()
		IJ.run(skeleton, "Skeletonize (2D/3D)", null);

		// IMPLEMENTATION 1
		// NOTE PROBLEM: I don't know how to take the result of Calculator plus and use
		// it int the plugin
//		ImageHandler skeletonBinary = ImageHandler.wrap(skeleton);		
//		skeletonBinary.show("Skeleton_Binary");		
//		labeledImage.show("Labeled_Image");
//		ImagePlus labeledSkeleton = null;
//		IJ.run("Calculator Plus", "i1=Skeleton_Binary i2=Labeled_Image operation=[Multiply: i2 = (i1*i2) x k1 + k2] k1="+skeletonBinary.getMax()+" k2=0 create");
//		IJ.selectWindow("Result");

		// IMPLEMENTATION 2 (reasonably fast)
		ImageInt labeledSkeletonImage = ImageInt.wrap(skeleton);
		ImagePlus labeledSkeletonImagePlus = labeledSkeletonImage.multiplyImage(labeledImage, (float) (1.0f / labeledSkeletonImage.getMax())).getImagePlus();
//		labeledSkeletonImagePlus.show(); 
		// this is a 32 bit image, and as soon as I wrap it to int it becomes 16-bit,
		// but then all the numbers change (rescales)... There might be a more elegant
		// fix...(TODO)
		labeledSkeletonImage = ImageInt.wrap(labeledSkeletonImagePlus);
		labeledSkeletonImage.multiplyByValue((float) (labeledImage.getMax() + 1) / 65536.0f);
		// labeledSkeletonImage.multiplyByValue((float)
		// (labeledImage.getMax()/labeledSkeletonImage.getMax()));
		labeledSkeletonImage.show(title + " - small structures might be missing");

		labelsSkeletons = Arrays.asList(labelImageTo3DVoxelArray(labeledSkeletonImage, 1));

		if (debug_output)
			System.out.println("Number of skeletons: " + labelsSkeletons.size());

		// ensure that no structures were accidentally completely removed during
		// skeletonization (which sometimes happens for small structures)
		for (int i = 0; i < (int) labeledImage.getMax(); i++) { // labelsSkeletons.size()
			if (labelsSkeletons.get(i).getVoxels().size() == 0) {
				Point3D centerPoint = new Object3DVoxels(labeledImage, (i + 1)).getCenterAsPoint();
				centerPoint.x = Math.round(centerPoint.x);
				centerPoint.y = Math.round(centerPoint.y);
				centerPoint.z = Math.round(centerPoint.z);
				Voxel3D centerVoxel = new Voxel3D(centerPoint, (i + 1));
				LinkedList<Voxel3D> tempList = new LinkedList<Voxel3D>();
				tempList.add(centerVoxel);
				labelsSkeletons.set(i, new Object3DVoxels(tempList));

				if (debug_output)
					System.out.println("SKELETON, Added for label " + (i + 1) + " at " + centerVoxel);
			}

			if (debug_output)
				System.out.println(
						"Label " + (i + 1) + " skeleton size " + labelsSkeletons.get(i).getVoxels().size() + " label volume " + (new Object3DVoxels(labeledImage, (i + 1)).getVoxels().size()));
		}

		// IMPLEMENTATION 3 (quite slow)
//		Object3DVoxels skeletonVoxels = new Object3DVoxels(ImageHandler.wrap(skeleton));
//		for (int i = 0; i < labelVoxels.length; ++i) {
//			labelsSkeletons.set(i, skeletonVoxels.getIntersectionObject(labelVoxels[i]));
//		}

		long endTime = System.currentTimeMillis();
		System.out.println("labelsToSkeleton() - Total execution time: " + (endTime - startTime) + "ms");

		return labelsSkeletons;
	}

	public Graph<Vector3D, DefaultEdge> skeletonToGraph(Object3DVoxels skeleton) {
		Graph<Vector3D, DefaultEdge> skeletonGraph = new DefaultUndirectedGraph<>(DefaultEdge.class);

		// add all vertices to graph
		List<Voxel3D> skeletonVoxels = skeleton.getVoxels();
		for (Voxel3D voxel : skeletonVoxels) {
			skeletonGraph.addVertex(voxel.getVector3D());
		}

		// connect graph vertices
		for (Voxel3D voxel : skeletonVoxels) {
			Vector3D vTemp = voxel.getVector3D();

			// look in 3x3x3 grid around current voxel
			for (int x = -1; x <= 1; ++x) {
				for (int y = -1; y <= 1; ++y) {
					for (int z = -1; z <= 1; ++z) {
						if (x == 0 && y == 0 && z == 0) // since this clearly exist already
						{
							Vector3D testVector = vTemp.add(new Vector3D(x, y, z));
							if (skeletonGraph.containsVertex(testVector)) {
								skeletonGraph.addEdge(vTemp, testVector);
							}
						}
					}
				}
			}

		}

		return skeletonGraph;
	}

	public List<Graph<Vector3D, DefaultEdge>> allSkeletonsToGraphs(List<Object3DVoxels> labelsSkeletons) {
		long startTime = System.currentTimeMillis();

		List<Graph<Vector3D, DefaultEdge>> labelsSkeletonGraphs = new ArrayList<Graph<Vector3D, DefaultEdge>>(labelsSkeletons.size());

		for (Object3DVoxels skeletonVoxels : labelsSkeletons) {
			labelsSkeletonGraphs.add(skeletonToGraph(skeletonVoxels));
		}

		if (debug_output)
			System.out.println("Number of graphs: " + labelsSkeletonGraphs.size());

		long endTime = System.currentTimeMillis();
		System.out.println("allSkeletonsToGraphs() - Total execution time: " + (endTime - startTime) + "ms");

		return labelsSkeletonGraphs;
	}

	class VectorPair implements Comparable<VectorPair> {
		public Vector3D vectorA;
		public Vector3D vectorB;
		public float distance;

		public VectorPair(Vector3D vA, Vector3D vB, float dist) {
			vectorA = vA;
			vectorB = vB;
			distance = dist;
		}

		public String toString() {
			return "A: " + this.vectorA + " B: " + this.vectorB + " Distance " + this.distance;
		}

		@Override
		public int compareTo(VectorPair otherVectorPair) {
			return Float.compare(distance, otherVectorPair.distance);
		}

//	    @Override
//	    public int compare(VectorPair vectorPairA, VectorPair vectorPairB) {
//	       return Integer.compare(vectorPairA.distance, vectorPairB.distance);
//	    }
	};

	public List<VectorPair> getClosestPointsBetweenStructures(Object3DVoxels structureA, Object3DVoxels structureB, int nClosest) {
		long startTime = System.currentTimeMillis();

		List<VectorPair> allVectorPairs = new ArrayList<VectorPair>(structureA.getVoxels().size() * structureB.getVoxels().size());

		for (Voxel3D voxelA : structureA.getVoxels()) {
			for (Voxel3D voxelB : structureB.getVoxels()) {
				Vector3D vA = voxelA.getVector3D();
				Vector3D vB = voxelB.getVector3D();
				allVectorPairs.add(new VectorPair(vA, vB, (float) vA.distance(vB)));
			}
		}
		Collections.sort(allVectorPairs); // , Collections.reverseOrder());

		long endTime = System.currentTimeMillis();
		System.out.println("getClosestPointsBetweenStructures() - Total execution time: " + (endTime - startTime) + "ms");

		for (VectorPair pair : allVectorPairs.subList(0, Math.min(nClosest, allVectorPairs.size()))) {
			if (debug_output)
				System.out.println(pair);
		}

		return allVectorPairs.subList(0, Math.min(nClosest, allVectorPairs.size()));
	}

	class GraphNode {
		public Vector3D location;
		public int relatedLabelInOtherFrame; // -1 means no match
		public int relatedLabelInThisFrame;
		public float distanceToRelated;

		public GraphNode(Vector3D location, int relatedLabel, int thisLabel, float distanceToRelated) {
			this.location = location;
			relatedLabelInOtherFrame = relatedLabel;
			relatedLabelInThisFrame = thisLabel;
			this.distanceToRelated = distanceToRelated;
		}

		public GraphNode(Vector3D location, int thisLabel) {
			this.location = location;
			relatedLabelInOtherFrame = -1; // means no match
			relatedLabelInThisFrame = thisLabel;
			distanceToRelated = Float.POSITIVE_INFINITY;
		}

		public String toString() {
			// +1 since background not included
			return "Node Location: " + this.location + " Related Label: " + (this.relatedLabelInOtherFrame + 1) + " This Label: " + (this.relatedLabelInThisFrame + 1) + " Distance to related: "
					+ this.distanceToRelated;
		}

		public boolean sameLocation(Vector3D other) {
			// return location == other;
			return Math.round(location.x) == Math.round(other.x) && Math.round(location.y) == Math.round(other.y) && Math.round(location.z) == Math.round(other.z);
		}

		public boolean sameLocation(GraphNode other) {
			return sameLocation(other.location);
		}
	};

	public GraphNode findDuplicateGraphNode(Graph<GraphNode, DefaultEdge> graph, Vector3D nodeLocation) {
		for (GraphNode graphNode : graph.vertexSet()) {
			if (graphNode.sameLocation(nodeLocation))
				return graphNode; // I'm assuming there would not be duplicate vertices that must be removed
		}

		return null;
	}

	public GraphNode findDuplicateGraphNode(Graph<GraphNode, DefaultEdge> graph, GraphNode node) {
		return findDuplicateGraphNode(graph, node.location);
	}

	// The index of the graph in the provided list matches the label number
	private Graph<GraphNode, DefaultEdge> matchGraphNodesBetweenFrames(List<Graph<Vector3D, DefaultEdge>> graphs_F1, List<Graph<Vector3D, DefaultEdge>> graphs_F2,
			List<List<Integer>> associatedLabelsBetweenFrames_F1toF2, float allowedDistance) {
		long startTime = System.currentTimeMillis();

		Graph<GraphNode, DefaultEdge> matchedGraphs_onF2 = new DefaultUndirectedGraph<>(DefaultEdge.class);

		// Loop through structures in Frame 1
		for (int labelNum_F1 = 0; labelNum_F1 < graphs_F1.size(); ++labelNum_F1) {
			Graph<Vector3D, DefaultEdge> labelGraph_F1 = graphs_F1.get(labelNum_F1);

			Graph<Voxel3D, DefaultEdge> tempAssociatedCompositeGraph = new DefaultUndirectedGraph<>(DefaultEdge.class);

			int associatedNodePairsMatched = 0;

			// only continue if there are in fact overlapping structures
			if (associatedLabelsBetweenFrames_F1toF2.get(labelNum_F1).size() == 0) {
				// TODO: WRITE THE END CODE
				// This label might be depolarising
				if (debug_output)
					System.out.println("CONTINUED for labelNum_F1 " + (labelNum_F1 + 1));
				continue;
			}

			// Loop through associated labeled structures in Frame 2 and create a graph
			// containing all the voxels of the Frame 2 structure
			// that is associated with the current label
			// input: associatedLabelsBetweenFrames_F1toF2, labelNum_F1, graphs_F2
			// output: tempAssociatedCompositeGraph - graph containing Frame 2 voxels
			// associated to Label in Frame 1
			for (int labelNum_F2 : associatedLabelsBetweenFrames_F1toF2.get(labelNum_F1)) {
				Graph<Vector3D, DefaultEdge> labelGraph_F2 = graphs_F2.get(labelNum_F2);

				for (Vector3D vertex_F2 : labelGraph_F2.vertexSet()) {
					// This code apparently already checks for duplicates
					Voxel3D voxel_F2 = new Voxel3D(vertex_F2.getPoint3f(), labelNum_F2);
					tempAssociatedCompositeGraph.addVertex(voxel_F2);
				}
			}

			//// DEBUG
			if (debug_output)
				System.out.println("Label F1 " + (labelNum_F1 + 1) + " has " + tempAssociatedCompositeGraph.vertexSet().size() + " voxels in F2");
			if (tempAssociatedCompositeGraph.vertexSet().size() == 0) {
				if (debug_output)
					System.out.println("associatedLabelsBetweenFrames_F1toF2 for label " + (labelNum_F1 + 1));
				for (int labelNum_F2 : associatedLabelsBetweenFrames_F1toF2.get(labelNum_F1)) {
					System.out.print((labelNum_F2 + 1) + ", ");
				}
				if (debug_output)
					System.out.println("");
			}
			////

			// These are only used in the situation where a node has no matches, then at the
			// end overwrite a non-optimal match but
			// at least ensure that there is a match, instead of missing events.
			GraphNode backupExistingNode = null;
			GraphNode backupNewGraphNode = null;

			// These two are related to the SAME LOCATION, but DIFFERENT Frame 1 labels,
			// allows to "swop out" options
			GraphNode backupExistingNodeShouldOverwrite = null;
			GraphNode backupNewGraphNodeShouldOverwrite = null;

			// this map stores the shortest distance for each associated Frame 1 label
			Map<Integer, Float> existingMatchedMinDistance = new HashMap<Integer, Float>();
			// this map stores the number of graph nodes for each associated Frame 1 label
			Map<Integer, Integer> numNodesPerF1Label = new HashMap<Integer, Integer>();

			// add matched node if close enough
			for (Voxel3D nodeVoxel_F2 : tempAssociatedCompositeGraph.vertexSet()) {
				// determine if this Frame 2 node has already got a "matched" node, and store
				// that existing node (from another Frame 1 Label)
				backupExistingNode = findDuplicateGraphNode(matchedGraphs_onF2, nodeVoxel_F2.getVector3D());

				// Find the closest Frame 1 voxel (for current label) to the currently
				// considered Frame 2 voxel (given that distance is within the allowedDistance
				// range)
				float distanceBetweenNodes = -1; // invalid value to test if a valid one was found
				for (Vector3D nodeVector_F1 : labelGraph_F1.vertexSet()) {
					float testDistanceBetweenNodes = (float) nodeVector_F1.distance(nodeVoxel_F2.getVector3D());
					if (testDistanceBetweenNodes <= allowedDistance) {
						if (distanceBetweenNodes == -1 || testDistanceBetweenNodes < distanceBetweenNodes) {
							distanceBetweenNodes = testDistanceBetweenNodes; // new minimum distance

							if (backupNewGraphNode == null || distanceBetweenNodes < backupNewGraphNode.distanceToRelated) {
								// the current "closest" Frame 1 voxel to the current Frame 2 voxel
								backupNewGraphNode = new GraphNode(nodeVoxel_F2.getVector3D(), labelNum_F1, (int) nodeVoxel_F2.value, distanceBetweenNodes);
							}
						}
					}
				}
				// if(debugOutput) System.out.println("Calculated closest distanceBetweenNodes =
				// " +
				// distanceBetweenNodes);

				// If it should be added, this is the new graph node that will be added to the
				// Matched Graph
				GraphNode newFrame2GraphNode = new GraphNode(nodeVoxel_F2.getVector3D(), labelNum_F1, (int) nodeVoxel_F2.value, distanceBetweenNodes);

				boolean shouldAdd = false; // don't add by default
				// check if it already exists in the graph for another Frame 1 label
				GraphNode existingNode = backupExistingNode;
				// if node already exist, but has a different label associated to it, but with a
				// greater distance, then choose to use smaller distance, hence remove
				if (existingNode != null) {
					// if(debugOutput) System.out.println("EXISTING " + existingNode);
					shouldAdd = false;
					// if there was a closest node && (new distance is less than existing distance
					// || existingNode in Frame 2 has no related structure in Frame 1)
					if (distanceBetweenNodes != -1 && (distanceBetweenNodes < existingNode.distanceToRelated || existingNode.distanceToRelated == -1)) {
						// This if statement prevents certain existing labels to be completely removed
						// since some other label is always closer if, however, it has been detected
						// before in some other location and that detected distance was smaller, now you
						// may happily replace it in this location
						float existingMinDistance = existingMatchedMinDistance.containsKey(existingNode.relatedLabelInOtherFrame)
								? existingMatchedMinDistance.get(existingNode.relatedLabelInOtherFrame)
								: Float.POSITIVE_INFINITY;
						if (existingMatchedMinDistance.containsKey(existingNode.relatedLabelInOtherFrame) && existingMinDistance < existingNode.distanceToRelated) {
							// remove existing node in same location
							if (debug_output)
								System.out.println("REPLACED: " + existingNode);
							matchedGraphs_onF2.removeVertex(existingNode);
							shouldAdd = true;
						}
						// if, however, it has already been detected before, and now the new distance is
						// smaller, then this is the one you should keep,
						// replace the previous one, and also save this case, to for future
						else if (existingMatchedMinDistance.containsKey(existingNode.relatedLabelInOtherFrame) && existingMinDistance >= existingNode.distanceToRelated) {
							if (debug_output)
								System.out.println("EXISTING " + existingNode);

							// remove old location and label from graph (this is not the same location as
							// the currently processed newFrame2GraphNode)
							if (debug_output)
								System.out.println("  was REPLACED: " + backupExistingNodeShouldOverwrite);
							matchedGraphs_onF2.removeVertex(backupExistingNodeShouldOverwrite);

							// replace that old location vertex with a new label
							matchedGraphs_onF2.addVertex(backupNewGraphNodeShouldOverwrite);
							if (debug_output)
								System.out.println("    ...WITH: " + backupNewGraphNodeShouldOverwrite);
							associatedNodePairsMatched++;

							// replace to store the new minimum distance
							existingMatchedMinDistance.replace(existingNode.relatedLabelInOtherFrame, existingNode.distanceToRelated);

							// store the currently considered graph nodes, in case a similar replacement
							// must be done in future
							backupExistingNodeShouldOverwrite = existingNode;
							backupNewGraphNodeShouldOverwrite = newFrame2GraphNode;
						}
						// in the default case, where existingMatchedMinDistance doesn't have the Frame
						// 1 label yet, never replace the first time if an existing
						// label is detected, since that might be the only place that label is stored.
						else {
							existingMatchedMinDistance.put(existingNode.relatedLabelInOtherFrame, existingNode.distanceToRelated);

							// store backup in case they need to be replaced in future
							backupExistingNodeShouldOverwrite = existingNode;
							backupNewGraphNodeShouldOverwrite = newFrame2GraphNode;

						}
					}
//					else {
//						 if(debugOutput) System.out.println("WON'T ADD label F1 " + (labelNum_F1 + 1) + " with distance " 
//								 + distanceBetweenNodes + " since existing is closer " +existingNode);
//					}
				} else if (newFrame2GraphNode.distanceToRelated != -1) { // No existing node, therefore NEW node, always
																			// add
					shouldAdd = true;
				}

				if (shouldAdd) {
					matchedGraphs_onF2.addVertex(newFrame2GraphNode);
					// if(debugOutput) System.out.println("ADDED: " + newFrame2GraphNode);
					associatedNodePairsMatched++;

					// if label already exists, then increment count
					if (numNodesPerF1Label.containsKey(newFrame2GraphNode.relatedLabelInOtherFrame)) {
						numNodesPerF1Label.replace(newFrame2GraphNode.relatedLabelInOtherFrame, numNodesPerF1Label.get(newFrame2GraphNode.relatedLabelInOtherFrame) + 1);
					} else // else create a new one and set count to 1
					{
						numNodesPerF1Label.put(newFrame2GraphNode.relatedLabelInOtherFrame, 1);
					}
				}

			}

			// if the newGraphNode was never added due to existing waiting to check for
			// duplicate, then add it now
			if (backupNewGraphNodeShouldOverwrite != null && findDuplicateGraphNode(matchedGraphs_onF2, backupNewGraphNodeShouldOverwrite) != null) {
				matchedGraphs_onF2.addVertex(backupNewGraphNodeShouldOverwrite);
				if (debug_output)
					System.out.println("\tADDED last New Graph node, since never added: " + backupNewGraphNode);
			}

			// if there were associatedLabels but now there are none due to distance
			// threshold removing them then mark this in some way to respond appropriately
			// to "early remove" potential events.
			// TODO: Check this "backup add" code if it still makes sense
			if (associatedNodePairsMatched == 0) {
				// TODO: This label has no associated fusion/fission event,
				// and therefore might depolarise or nothing (but remember, this is only due to
				// distance threshold)
				if (debug_output)
					System.out.println("\tERROR NO MATCHING for Label F1 " + (labelNum_F1 + 1));

				if (backupExistingNode != null)
					matchedGraphs_onF2.removeVertex(backupExistingNode);

				if (backupNewGraphNode != null) {
					matchedGraphs_onF2.addVertex(backupNewGraphNode);
					if (debug_output)
						System.out.println("\t...ADDED BACKUP: " + backupNewGraphNode);
					associatedNodePairsMatched = 1;
				} else {
					if (debug_output)
						System.out.println("\t... THERE WAS NO BACKUP TO ADD");
				}
			}
		}

		//// DEBUG: Find all potential duplicates which were added (which shouldn't be
		//// the case if the above code did its job)
		List<GraphNode> nodesToRemove = new ArrayList<GraphNode>();
		for (GraphNode nodeA : matchedGraphs_onF2.vertexSet()) {
			for (GraphNode nodeB : matchedGraphs_onF2.vertexSet()) {
				if (nodeA.sameLocation(nodeB) && (nodeA.distanceToRelated == -1 || nodeB.distanceToRelated == -1) && (nodeA.distanceToRelated != nodeB.distanceToRelated)) {
					if (debug_output)
						System.out.println("MATCHED " + nodeA + " and " + nodeB);
					if (nodeA.distanceToRelated == -1) {
						if (debug_output)
							System.out.println("   ...therefore REMOVED " + nodeA);
						nodesToRemove.add(nodeA);
					}

					if (nodeB.distanceToRelated == -1) {
						if (debug_output)
							System.out.println("   ...therefore REMOVED " + nodeB);
						nodesToRemove.add(nodeB);
					}
				}
			}
		}
		// remove all the ones that were picked up to be removed
		for (GraphNode removeNode : nodesToRemove) {
			matchedGraphs_onF2.removeVertex(removeNode);
		}
		////
//		
//		for (GraphNode node : matchedGraphs_onF2.vertexSet()) {
//			GraphNode existingNode = findDuplicateGraphNode(matchedGraphs_onF2, node);
//			if (existingNode != null && node.distanceToRelated == existingNode.distanceToRelated) {
//				if(debugOutput) System.out.println("FOUND DUPLICATE " + node);
//				if(existingNode.distanceToRelated == -1)
//					matchedGraphs_onF2.removeVertex(existingNode);
//			}
//		}

		long endTime = System.currentTimeMillis();
		System.out.println("matchGraphNodesBetweenFrames() - Total execution time: " + (endTime - startTime) + "ms");

		return automaticallyConnectGraphNodes(matchedGraphs_onF2);
	}

	// This function takes a non-connected graph (with only vertices) and
	// automatically connect the nodes to adjacent
	// nodes based on the vertex location.
	public Graph<GraphNode, DefaultEdge> automaticallyConnectGraphNodes(Graph<GraphNode, DefaultEdge> inputGraph) {
		long startTime = System.currentTimeMillis();

		for (GraphNode node1 : inputGraph.vertexSet()) {
			for (GraphNode node2 : inputGraph.vertexSet()) {
				Vector3D diff = new Vector3D(Math.abs(node1.location.x - node2.location.x), Math.abs(node1.location.y - node2.location.y), Math.abs(node1.location.z - node2.location.z));

				if (0 <= diff.x && diff.x <= 1 && 0 <= diff.y && diff.y <= 1 && 0 <= diff.z && diff.z <= 1) {
					if (!inputGraph.containsEdge(node1, node2)) {
						inputGraph.addEdge(node1, node2);
						// if(debugOutput) System.out.println("ADDED EDGE: from " + node1 + " to " +
						// node2);
					}

				}
			}
		}

		long endTime = System.currentTimeMillis();
		System.out.println("automaticallyConnectGraphNodes() - Total execution time: " + (endTime - startTime) + "ms");

		return inputGraph;
	}

	// This renders the graph as a an Image (mainly for debugging purposes
	public void showGraphNodesAsImage(Graph<GraphNode, DefaultEdge> inputGraph, int width, int height, int nSlices, boolean showMatched, String title) {
		ImagePlus graphImage = NewImage.createImage("Graph Image", width, height, nSlices, 16, NewImage.FILL_BLACK);

		ImageHandler ih = ImageHandler.wrap(graphImage);

		for (GraphNode node : inputGraph.vertexSet()) {
			Vector3D v = node.location;

			if (showMatched && node.distanceToRelated != -1)
				ih.setPixel((int) v.x, (int) v.y, (int) v.z, node.relatedLabelInOtherFrame + 1); // +1 since the
																									// skeleton does not
																									// include
																									// background
			else if (!showMatched && node.distanceToRelated == -1)
				ih.setPixel((int) v.x, (int) v.y, (int) v.z, Integer.MAX_VALUE); // +1 since the skeleton does not
																					// include background

		}

		ih.show(title);

//		
//		ImageProcessor graphImageProc = graphImage.getProcessor();
//		
//		for(GraphNode node: inputGraph.vertexSet())
//		{
//			Vector3D v = node.location;
//			graphImageProc.setSliceNumber((int) v.z);
//			graphImageProc.putPixel((int)v.x, (int)v.y, node.relatedLabelInOtherFrame);
//		}
//		
//		graphImage.setProcessor(graphImageProc);
//		graphImage.show();

	}

	// run along the graph until I find a transition between two labels and that is
	// presumably the event location.
	public List<Vector3D> findEvents(Graph<GraphNode, DefaultEdge> inputGraph, List<Graph<Vector3D, DefaultEdge>> labelsSkeletonGraphs, boolean removeDuplicates, float duplicateDistance) {
		long startTime = System.currentTimeMillis();

		List<Vector3D> eventList = new ArrayList<Vector3D>();
		Map<Integer, Vector3D> alreadyMatched = new HashMap<Integer, Vector3D>();

		System.out.println("Input graph size: " + inputGraph.vertexSet().size());

		// loop through each node/voxel combination
		for (GraphNode node1 : inputGraph.vertexSet()) {
			for (GraphNode node2 : inputGraph.vertexSet()) {
				// calculate the vector between the two nodes
				Vector3D diff = new Vector3D(Math.abs(node1.location.x - node2.location.x), Math.abs(node1.location.y - node2.location.y), Math.abs(node1.location.z - node2.location.z));
				// the two voxels are neighbouring
				if (0 <= diff.x && diff.x <= 1 && 0 <= diff.y && diff.y <= 1 && 0 <= diff.z && diff.z <= 1) {
					// .. and not the same voxel, AND is at a transition point
					if (node1.relatedLabelInOtherFrame != node2.relatedLabelInOtherFrame) {
						System.out.println("EVENT");

						Vector3D eventLocation = getHalfwayPoint(node1.location, node2.location, true);

						if (!eventList.contains(eventLocation)) {
							eventList.add(eventLocation);
							// if (debug_output)
							System.out.println("EVENT LOCATION " + eventLocation + " label 1 " + (node1.relatedLabelInOtherFrame + 1) + " label 2 " + (node2.relatedLabelInOtherFrame + 1)
									+ " with diff " + diff.toString());
						}
					}
				}
			}
		}

		// TODO: Currently no duplicates are removed
		if (removeDuplicates) {
			List<Vector3D> duplicateRemovedEventList = new ArrayList<Vector3D>();
			for (Vector3D vec : eventsToGraph(eventList, removeDuplicates, duplicateDistance).vertexSet()) {
				duplicateRemovedEventList.add(vec);
			}
			eventList = duplicateRemovedEventList;
//			List<List<Vector3D>> nearbyEventsList = new ArrayList<List<Vector3D>>(eventList.size());
//			for (Vector3D eventLocation : eventList) {
//				List<Vector3D> nearbyEvents = new ArrayList<Vector3D>();
//				for (Vector3D otherEventLocaiton : eventList) {
//					double dist = eventLocation.distance(otherEventLocaiton);
//					if (0 < dist && dist <= duplicateDistance) {
//						nearbyEvents.add(otherEventLocaiton);
//					}
//				}
//				nearbyEvents.add(eventLocation); // include same event as well (but only once)
//				nearbyEventsList.add(nearbyEvents);
//			}
//
//			for (int i = 0; i < nearbyEventsList.size(); i++) {
//				if(debugOutput) System.out.println(i);
//				List<Vector3D> nearbyEvents = nearbyEventsList.get(i);
//				if (nearbyEvents.size() > 1) { // since normal event also included
//					Vector3D averageEventLocation = nearbyEvents.get(0);
//					eventList = removeEventLocation(eventList, nearbyEvents.get(0));
//					for (int j = 1; j < nearbyEvents.size(); j++) {
//						averageEventLocation = averageEventLocation.add(nearbyEvents.get(j));
//						eventList = removeEventLocation(eventList, nearbyEvents.get(j));
//					}
//
//					averageEventLocation = averageEventLocation.multiply(1.0 / nearbyEvents.size());
//					averageEventLocation.x = Math.round(averageEventLocation.x);
//					averageEventLocation.y = Math.round(averageEventLocation.y);
//					averageEventLocation.z = Math.round(averageEventLocation.z);
//
//					eventList.add(averageEventLocation);
//
//					if(debugOutput) System.out.println("ADDED AVERAGE " + averageEventLocation);
//				}
//			}
		}

		long endTime = System.currentTimeMillis();
		System.out.println("findEvents() - Total execution time: " + (endTime - startTime) + "ms");

		return eventList;
	}

	// This function receives a 3D point from one labeled structure and receives the
	// entire labeled structure of the other label
	// it then finds the closest point in the other structures to the first point
	// and calculates the halfway point from those two.
	public Vector3D getMinHalfwayPoint(Vector3D point, Graph<Vector3D, DefaultEdge> labeledGraph, boolean makeIntegerOutput) {
		Vector3D halfwayPoint = null;
		float distance = Float.POSITIVE_INFINITY;

		for (Vector3D vertex : labeledGraph.vertexSet()) {
			if (point.distance(vertex) < distance) {
				distance = (float) point.distance(vertex);
				halfwayPoint = point.add(vertex).multiply(0.5);
			}
		}

		if (makeIntegerOutput) {
			halfwayPoint.x = Math.round(halfwayPoint.x);
			halfwayPoint.y = Math.round(halfwayPoint.y);
			halfwayPoint.z = Math.round(halfwayPoint.z);
		}
		return halfwayPoint;
	}

	public Vector3D getHalfwayPoint(Vector3D point_L1, Vector3D point_L2, boolean makeIntegerOutput) {
		Vector3D halfwayPoint = point_L1.add(point_L2).multiply(0.5);
		;

		if (makeIntegerOutput) {
			halfwayPoint.x = Math.round(halfwayPoint.x);
			halfwayPoint.y = Math.round(halfwayPoint.y);
			halfwayPoint.z = Math.round(halfwayPoint.z);
		}
		return halfwayPoint;
	}

	public List<Vector3D> findDepolarisationEvents(List<List<Integer>> associatedLabelsBetweenFrames, Object3DVoxels[] labeledVoxels) {
		long startTime = System.currentTimeMillis();

		List<Vector3D> eventList = new ArrayList<Vector3D>();

		for (int i = 0; i < associatedLabelsBetweenFrames.size(); ++i) {
			if (associatedLabelsBetweenFrames.get(i).size() == 0) {
				if (debug_output)
					System.out.println("Depolaristaion detected at label " + (i + 1));

				Point3D centerPoint = labeledVoxels[i].getCenterAsPoint();
				centerPoint.x = Math.round(centerPoint.x);
				centerPoint.y = Math.round(centerPoint.y);
				centerPoint.z = Math.round(centerPoint.z);
				Vector3D centerVector = new Vector3D(centerPoint);
//				Voxel3D centerVoxel = new Voxel3D(centerPoint, (i + 1));
				eventList.add(centerVector);
			}
		}

		long endTime = System.currentTimeMillis();
		System.out.println("findEvents() - Total execution time: " + (endTime - startTime) + "ms");

		return eventList;
	}

	public Graph<Vector3D, DefaultEdge> eventsToGraph(List<Vector3D> eventsList, boolean removeDuplicates, float duplicateDistance) {
		Graph<Vector3D, DefaultEdge> eventsGraph = new DefaultUndirectedGraph<>(DefaultEdge.class);

		// add all vertices to graph
		for (Vector3D eventLoc : eventsList) {
			eventsGraph.addVertex(eventLoc);
		}

		// connect graph vertices
		for (Vector3D eventLocA : eventsGraph.vertexSet()) {
			for (Vector3D eventLocB : eventsGraph.vertexSet()) {
				if (eventLocA.distance(eventLocB) <= duplicateDistance) {
					eventsGraph.addEdge(eventLocA, eventLocB);
				}
			}
		}

		if (removeDuplicates) {
			boolean didRemove = false;
			do {
				didRemove = false;
				for (Vector3D eventLocA : eventsGraph.vertexSet()) {
					if (eventsGraph.edgesOf(eventLocA).size() > 0) {
						didRemove = true;

						List<Vector3D> toRemove = new ArrayList<Vector3D>();
						Iterator<Vector3D> iterator = new DepthFirstIterator<>(eventsGraph, eventLocA);
						while (iterator.hasNext()) {
							toRemove.add(iterator.next());
						}
						Vector3D averageEventLocation = toRemove.get(0);
						eventsGraph.removeVertex(toRemove.get(0));
						if (debug_output)
							System.out.println("REMOVED duplicate " + toRemove.get(0));
						int count;
						for (count = 1; count < toRemove.size(); count++) {
							Vector3D removeVertex = toRemove.get(count);
							averageEventLocation = averageEventLocation.add(removeVertex);

							eventsGraph.removeVertex(removeVertex);
							if (debug_output)
								System.out.println("REMOVED " + removeVertex);
						}

						averageEventLocation = averageEventLocation.multiply(1.0 / count);
						averageEventLocation.x = Math.round(averageEventLocation.x);
						averageEventLocation.y = Math.round(averageEventLocation.y);
						averageEventLocation.z = Math.round(averageEventLocation.z);

						eventsGraph.addVertex(averageEventLocation);
						if (debug_output)
							System.out.println("ADDED AVERAGE of duplicates " + averageEventLocation);
						break; // since my loop is no longer valid with vertices added and removed
					}
				}
			} while (didRemove);
		}

		return eventsGraph;
	}

	public List<Vector3D> removeEventLocation(List<Vector3D> eventList, Vector3D location) {
		int initialSize = eventList.size();
		for (int i = 0; i < eventList.size(); ++i) {
			if (eventList.get(i).distance(location) == 0) {
				eventList.remove(i);
				if (debug_output)
					System.out.println("Succesfully removed " + location);
				break;
			}
		}

		if (initialSize == eventList.size()) {
			if (debug_output)
				System.out.println("NO DUPLICATES DEMOVED for " + location);
		}

		return eventList;
	}

	public ImageInt eventsToImage(List<Vector3D> eventList, int width, int height, int nSlices, String title) {
		ImagePlus graphImage = NewImage.createImage("Graph Image", width, height, nSlices, 16, NewImage.FILL_BLACK);

		ImageInt image = ImageInt.wrap(graphImage);

		for (Vector3D eL : eventList) {
			image.setPixel((int) eL.x, (int) eL.y, (int) eL.z, 1);
		}

		image.show(title);

		return image;
	}

	// This is primarily for debugging purposes
	public void printAllWindowTitles() {
		for (String title : WindowManager.getImageTitles()) {
			System.out.println(title);
		}
	}

	/**
	 * This main function serves for development purposes. It allows you to run the
	 * plugin immediately out of your integrated development environment (IDE).
	 *
	 * @param args whatever, it's ignored
	 * @throws Exception
	 */
	public static void main(final String... args) throws Exception {
		// create the ImageJ application context with all available services
		final ImageJ ij = new ImageJ();
//        ij.ui().showUI();

		// ask the user for a file to open
		final File file = ij.ui().chooseFile(null, "open");

//		open("C:/Users/rptheart/Dropbox/Research/MitoMorph/MEL_fiji/Frame1_Thresholded.tif");
//		open("C:/Users/rptheart/Dropbox/Research/MitoMorph/MEL_fiji/Frame2_Thresholded.tif");

		if (file != null) {
			// load the dataset
			final Dataset dataset = ij.scifio().datasetIO().open(file.getPath());

			// show the image
			ij.ui().show(dataset);

			// invoke the plugin
			ij.command().run(MEL_Modules.class, true);

		}
	}

}
