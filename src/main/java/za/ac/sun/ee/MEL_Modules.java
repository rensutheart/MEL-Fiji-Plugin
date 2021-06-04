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
import net.imagej.legacy.convert.ImagePlusToDatasetConverter;
import net.imagej.ops.OpService;
import net.imglib2.Cursor;
import net.imglib2.IterableInterval;
import net.imglib2.RandomAccess;
import net.imglib2.RandomAccessibleInterval;
import net.imglib2.img.Img;
import net.imglib2.img.display.imagej.ImageJFunctions;
import net.imglib2.ops.parse.token.Int;
import net.imglib2.type.Type;
import net.imglib2.type.numeric.IntegerType;
import net.imglib2.type.numeric.RealType;
import net.imglib2.type.numeric.integer.UnsignedIntType;

import org.jgrapht.Graph;
import org.jgrapht.graph.DefaultEdge;
import org.jgrapht.graph.DefaultUndirectedGraph;
import org.jgrapht.traverse.DepthFirstIterator;
import org.scijava.command.Command;
import org.scijava.plugin.Parameter;
import org.scijava.plugin.Plugin;
import org.scijava.ui.UIService;

import ij.IJ;
import ij.ImagePlus;
import ij.WindowManager;
import ij.gui.NewImage;
import ij.process.ImageProcessor;
import mcib3d.geom.Object3D;
import mcib3d.geom.Object3DLabel;
import mcib3d.geom.Object3DVoxels;
import mcib3d.geom.Object3D_IJUtils;
import mcib3d.geom.Point3D;
import mcib3d.geom.Vector3D;
import mcib3d.geom.Voxel3D;
import mcib3d.image3d.ImageByte;
import mcib3d.image3d.ImageHandler;
import mcib3d.image3d.ImageInt;
import mcib3d.image3d.ImageLabeller;
import mcib3d.image3d.ImagePlus_Utils;

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

/**
 * This example illustrates how to create an ImageJ {@link Command} plugin.
 * <p>
 * The code here is a simple Gaussian blur using ImageJ Ops.
 * </p>
 * <p>
 * You should replace the parameter fields with your own inputs and outputs, and
 * replace the {@link run} method implementation with your own logic.
 * </p>
 */
@Plugin(type = Command.class, menuPath = "Plugins>MEL Process")
public class MEL_Modules<T extends RealType<T>> implements Command {
	//
	// Feel free to add more parameters here...
	//

	@Parameter
	private Dataset currentData;

	@Parameter
	private UIService uiService;

	@Parameter
	private OpService opService;

	@Parameter
	private int minStructureVolume = 20; // TODO: This can be dependent on Voxel size (see Mitochondrial Analyzer macro)
	
	@Parameter
	private float minOverlapPercentage = 0.1f;
		
	@Parameter
	private float skeletonDistanceThreshold = 20; 

	@Override
	public void run() {
		long startTime = System.currentTimeMillis();
		// I'm assuming the input images are Pre-processed and thresholded.

		/*
		 * LOAD THRESHOLDED FRAMES
		 */
		String title_Frame1 = "Frame1_Thresholded.tif";
		String title_Frame2 = "Frame2_Thresholded.tif";

		ImagePlus imagePlus_Frame1 = WindowManager.getImage(title_Frame1);
		ImagePlus imagePlus_Frame2 = WindowManager.getImage(title_Frame2);

		/*
		 * LABEL MITOCHONDRIAL STRUCTURES
		 */
		// Get the segmented labeled mitochondrial structures
		ImageInt labels_F1 = getLabeledImage(imagePlus_Frame1, minStructureVolume);
		ImageInt labels_F2 = getLabeledImage(imagePlus_Frame2, minStructureVolume);

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
		List<List<Integer>> associatedLabelsBetweenFrames_F1toF2 = getAssociatedLabelsBetweenFrames(overlappingVolumes);
		List<List<Integer>> associatedLabelsBetweenFrames_F2toF1 = getAssociatedLabelsBetweenFrames(transposeMatrix(overlappingVolumes));

		// I don't seem to need this anymore for the new approach
//		List<List<Integer>> associatedLabelsWithinFrame_F1 = getAssociatedLabelsWithinFrame(associatedLabelsBetweenFrames_F1toF2, associatedLabelsBetweenFrames_F2toF1);
//		List<List<Integer>> associatedLabelsWithinFrame_F2 = getAssociatedLabelsWithinFrame(associatedLabelsBetweenFrames_F2toF1, associatedLabelsBetweenFrames_F1toF2);

		List<List<Float>> percentageOverlap_F1toF2 = getRelativePercentageOverlap(overlappingVolumes, associatedLabelsBetweenFrames_F1toF2);
		List<List<Float>> percentageOverlap_F2toF1 = getRelativePercentageOverlap(transposeMatrix(overlappingVolumes), associatedLabelsBetweenFrames_F2toF1);
		
		List<List<Integer>> reducedAssociatedLabelsBetweenFrames_F1toF2 = reduceAssociateLabels(associatedLabelsBetweenFrames_F1toF2, percentageOverlap_F1toF2, minOverlapPercentage);
		List<List<Integer>> reducedAssociatedLabelsBetweenFrames_F2toF1 = reduceAssociateLabels(associatedLabelsBetweenFrames_F2toF1, percentageOverlap_F2toF1, minOverlapPercentage);

		// Display Results
		System.out.println("associatedLabelsBetweenFrames_F1toF2.size " + associatedLabelsBetweenFrames_F1toF2.size());
		System.out.println("associatedLabelsBetweenFrames_F2toF1.size " + associatedLabelsBetweenFrames_F2toF1.size());

//		System.out.println("associatedLabelsWithinFrame_F1.size " + associatedLabelsWithinFrame_F1.size());
//		System.out.println("associatedLabelsWithinFrame_F2.size " + associatedLabelsWithinFrame_F2.size());

		System.out.println("percentageOverlap_F1toF2.size " + percentageOverlap_F1toF2.size());
		System.out.println("percentageOverlap_F2toF1.size " + percentageOverlap_F2toF1.size());
		
		System.out.println("reducedAssociatedLabelsBetweenFrames_F1toF2.size " + reducedAssociatedLabelsBetweenFrames_F1toF2.size());
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
		System.out.println("labelsSkeletons_F1.size " + labelsSkeletons_F1.size());
		System.out.println("labelsSkeletons_F2.size " + labelsSkeletons_F2.size());

		System.out.println("labelsSkeletonGraphs_F1.size " + labelsSkeletonGraphs_F1.size());
		System.out.println("labelsSkeletonGraphs_F2.size " + labelsSkeletonGraphs_F2.size());

		/*
		 * DETECT FUSION
		 */
		// Find the matched skeleton F1 to F2
		System.out.println("\nFIND FUSION");
		Graph<GraphNode, DefaultEdge> matchedGraphs_onF2 = matchGraphNodesBetweenFrames(labelsSkeletonGraphs_F1, labelsSkeletonGraphs_F2, reducedAssociatedLabelsBetweenFrames_F1toF2, skeletonDistanceThreshold);
		showGraphNodesAsImage(matchedGraphs_onF2, labels_F1.sizeX, labels_F1.sizeY, labels_F1.sizeZ, true, "Matched graph Fusion");
//		showGraphNodesAsImage(matchedGraphs_onF2, labels_F1.sizeX, labels_F1.sizeY, labels_F1.sizeZ, false, "Unmatched graph");

		// Find all the events F1 to F2 (fusion)
		List<Vector3D> fusionEventLocations = findEvents(matchedGraphs_onF2, labelsSkeletonGraphs_F1, true, 10);
		ImageInt fusionEventsImage = eventsToImage(fusionEventLocations, labels_F1.sizeX, labels_F1.sizeY, labels_F1.sizeZ, "Fusion events");

		/*
		 * DETECT FISSION
		 */
		// Find the matched skeleton F2 to F1
		System.out.println("\nFIND FISSION");
		Graph<GraphNode, DefaultEdge> matchedGraphs_onF1 = matchGraphNodesBetweenFrames(labelsSkeletonGraphs_F2, labelsSkeletonGraphs_F1, reducedAssociatedLabelsBetweenFrames_F2toF1, skeletonDistanceThreshold);
		showGraphNodesAsImage(matchedGraphs_onF1, labels_F1.sizeX, labels_F1.sizeY, labels_F1.sizeZ, true, "Matched graph Fission");
//		showGraphNodesAsImage(matchedGraphs_onF1, labels_F1.sizeX, labels_F1.sizeY, labels_F1.sizeZ, false, "Unmatched graph");

		// Find all the events F2 to F1 (fission)
		List<Vector3D> fissionEventLocations = findEvents(matchedGraphs_onF1, labelsSkeletonGraphs_F2, true, 10);
		ImageInt fissionEventsImage = eventsToImage(fissionEventLocations, labels_F1.sizeX, labels_F1.sizeY, labels_F1.sizeZ, "Fission events");

		
		/*
		 * DETECT DEPOLARISATION
		 */
		System.out.println("\nFIND DEPOLARISATION");
		List<Vector3D> depolarisationEventLocations = findDepolarisationEvents(reducedAssociatedLabelsBetweenFrames_F1toF2, labels_F1);
		ImageInt depolarisationEventsImage = eventsToImage(depolarisationEventLocations, labels_F1.sizeX, labels_F1.sizeY, labels_F1.sizeZ, "Depolarisation events");
		
		
		/*
		 * // TEST CODE: Note the - 1, that is since background is removed, now label 1
		 * is // at index 0. int label_F1_from = 40; int label_F1_to =
		 * associatedLabelsWithinFrame_F1.get(label_F1_from - 1).get(0) + 1;
		 * 
		 * System.out.println("Structures associated within with label " +
		 * label_F1_from); for (int label_to :
		 * associatedLabelsWithinFrame_F1.get(label_F1_from - 1)) {
		 * System.out.println(label_to + 1); }
		 * 
		 * System.out.println("Finding 5 closest points between label " + label_F1_from
		 * + " and " + label_F1_to); // IMPORTANT NOTE: labelsSkeletons matches the
		 * Label image and therefore // includes the background, hence no -1 for label
		 * number List<VectorPair> closestVectorPair =
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
				if (volumeOverlap != 0)
					System.out.println("Label F1 " + (label_F1 + 1) + "  Label F2 " + (label_F2 + 1) + "  Volume Overlap " + volumeOverlap);
			}
		}

		long endTime = System.currentTimeMillis();
		System.out.println("getOverlappingVolumes() - Total execution time: " + (endTime - startTime) + "ms");

		return overlappingVolumes;
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

	public List<List<Integer>> getAssociatedLabelsBetweenFrames(int[][] overlappingVolumes) {
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
					System.out.println("Associated between " + (i + 1) + " and " + (j + 1));
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
			System.out.println("Frame 1 Label: " + (labelIndex_F1 + 1));
			Set<Integer> tempSet = new HashSet<Integer>();

			if (list_F1.size() > 0) {
				for (int labelNum_F2 : list_F1) {
					System.out.println("\t Frame 2 Label: " + (labelNum_F2 + 1));
					for (int labelNum_F1 : associatedLabelsBetweenFrames_F2.get(labelNum_F2)) {
						if (labelIndex_F1 != labelNum_F1) // ensure the label is not associated to itself
						{
							tempSet.add(labelNum_F1);
							// This assumes a startLabel of 1
							System.out.println("\t\t Label " + (labelIndex_F1 + 1) + " is associated with label " + (labelNum_F1 + 1));
						}
					}
				}
			} else {
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
				
				//System.out.println("Percentage overlap: " + labelVolumeOverlaps_F1.get(i) + " (" + labelNum_F1 + ", " + i + ")");
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
				if(percentageOverlap_F1toF2.get(i).get(j) > cutoffPercentage)
					tempList.add(associatedLabelsBetweenFrames_F1toF2.get(i).get(j));
				else
					System.out.println("Removed match between label " +(i+1) + " and label " + (j+1) + " in the other frame with percentage " + percentageOverlap_F1toF2.get(i).get(j));
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

		// Use the labeled image to ensure small structures are removed already
		// (minStructureSize)
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

		// IMPLEMENTATION 2 (resonably fast)
		ImageInt labeledSkeletonImage = ImageInt.wrap(skeleton);
		ImagePlus labeledSkeletonImagePlus = labeledSkeletonImage.multiplyImage(labeledImage, (float) (1.0f / labeledSkeletonImage.getMax())).getImagePlus();
//		labeledSkeletonImagePlus.show(); 
		// this is a 32 bit image, and as soon as I wrap it to int it becomes 16-bit,
		// but then all the numbers change (rescales)...
		labeledSkeletonImage = ImageInt.wrap(labeledSkeletonImagePlus);
		labeledSkeletonImage.multiplyByValue((float) (labeledImage.getMax() + 1) / 65536.0f);
		// labeledSkeletonImage.multiplyByValue((float)
		// (labeledImage.getMax()/labeledSkeletonImage.getMax()));
		labeledSkeletonImage.show(title + " - small structures might be missing");

		labelsSkeletons = Arrays.asList(labelImageTo3DVoxelArray(labeledSkeletonImage, 1));

		System.out.println("Number of skeletons: " + labelsSkeletons.size());

		// ensure that no structures were accidentally completely removed during
		// skeletonization (which sometimes happens for small structures)
		for (int i = 0; i < labelsSkeletons.size(); i++) {
			if (labelsSkeletons.get(i).getVoxels().size() == 0) {
				Point3D centerPoint = new Object3DVoxels(labeledImage, (i + 1)).getCenterAsPoint();
				centerPoint.x = Math.round(centerPoint.x);
				centerPoint.y = Math.round(centerPoint.y);
				centerPoint.z = Math.round(centerPoint.z);
				Voxel3D centerVoxel = new Voxel3D(centerPoint, (i + 1));
				LinkedList<Voxel3D> tempList = new LinkedList<Voxel3D>();
				tempList.add(centerVoxel);
				labelsSkeletons.set(i, new Object3DVoxels(tempList));

				System.out.println("SKELETON, Added for label " + (i + 1) + " at " + centerVoxel);
			}

			System.out.println("Label " + (i + 1) + " skeleton size " + labelsSkeletons.get(i).getVoxels().size() + " label volume " + (new Object3DVoxels(labeledImage, (i + 1)).getVoxels().size()));
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
			return "Node Location: " + this.location + " Related Label: " + (this.relatedLabelInOtherFrame + 1) + 
					" This Label: " + (this.relatedLabelInThisFrame + 1) + " Distance to related: " + this.distanceToRelated;
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
				System.out.println("CONTINUED for labelNum_F1 " + (labelNum_F1 + 1));
				continue;
			}

			// Loop through associated labeled structures in Frame 2
			for (int labelNum_F2 : associatedLabelsBetweenFrames_F1toF2.get(labelNum_F1)) {
				Graph<Vector3D, DefaultEdge> labelGraph_F2 = graphs_F2.get(labelNum_F2);

				for (Vector3D vertex_F2 : labelGraph_F2.vertexSet()) {
					// This code apparently already checks for duplicates
					Voxel3D voxel_F2 = new Voxel3D(vertex_F2.getPoint3f(), labelNum_F2);
					tempAssociatedCompositeGraph.addVertex(voxel_F2);
				}
			}

			//// DEBUG
			System.out.println("Label F1 " + (labelNum_F1 + 1) + " has " + tempAssociatedCompositeGraph.vertexSet().size() + " voxels in F2");
			if (tempAssociatedCompositeGraph.vertexSet().size() == 0) {
				System.out.println("associatedLabelsBetweenFrames_F1toF2 for label " + (labelNum_F1 + 1));
				for (int labelNum_F2 : associatedLabelsBetweenFrames_F1toF2.get(labelNum_F1)) {
					System.out.print((labelNum_F2 + 1) + ", ");
				}
				System.out.println("");
			}
			////

			// These are only used in the situation where a node has no matches, then at the end overwrite a non-optimal match but
			// at least ensure that there is a match, instead of missing events.
			GraphNode backupExistingNode = null;
			GraphNode backupNewGraphNode = null;
			
			GraphNode backupExistingNodeShouldOverwrite = null;
			GraphNode backupNewGraphNodeShouldOverwrite = null;
			Map<Integer, Float> existingMatchedMinDistance = new HashMap<Integer, Float>();
			
			// add matched node if close enough
			for (Voxel3D nodeVector_F2 : tempAssociatedCompositeGraph.vertexSet()) {
				float distanceBetweenNodes = -1; // invalid value to test if a valid one was found
				for (Vector3D nodeVector_F1 : labelGraph_F1.vertexSet()) {
					float testDistanceBetweenNodes = (float) nodeVector_F1.distance(nodeVector_F2.getVector3D());
					if (testDistanceBetweenNodes <= allowedDistance) {
						if (distanceBetweenNodes == -1 || testDistanceBetweenNodes < distanceBetweenNodes) {
							distanceBetweenNodes = testDistanceBetweenNodes;
							
							if(backupNewGraphNode == null || distanceBetweenNodes < backupNewGraphNode.distanceToRelated)
							{
								backupNewGraphNode = new GraphNode(nodeVector_F2.getVector3D(), labelNum_F1, (int) nodeVector_F2.value, distanceBetweenNodes);
								backupExistingNode = findDuplicateGraphNode(matchedGraphs_onF2, nodeVector_F2.getVector3D());
							}
						}
					}
				}
				
				 System.out.println("Calculated distanceBetweenNodes = " + distanceBetweenNodes);

				boolean shouldAdd = true; // add by default
				// check if it already exists in the graph for another label
				// if node already exist, but has a different label associated to it, but with a
				// greater distance, then choose to use smaller distance, hence remove
				GraphNode existingNode = findDuplicateGraphNode(matchedGraphs_onF2, nodeVector_F2.getVector3D());
				if (existingNode != null) {
					System.out.println("EXISTING " + existingNode);
					shouldAdd = false;
					if (distanceBetweenNodes != -1 && (distanceBetweenNodes < existingNode.distanceToRelated || existingNode.distanceToRelated == -1)) {
						// This if statement prevents certain existing labels to be completely removed since some other label is always closer
						// if, however, it has been detected before, and the detected distance was smaller, now you may hapily replace it
						if(existingMatchedMinDistance.containsKey(existingNode.relatedLabelInOtherFrame) && 
								existingMatchedMinDistance.get(existingNode.relatedLabelInOtherFrame) < existingNode.distanceToRelated)
						{
							System.out.println("REPLACED: " + existingNode);
							matchedGraphs_onF2.removeVertex(existingNode);
							shouldAdd = true;
						}
						// if it has already been detected before, and now the new distance is smaller, then this is the one you should keep, 
						// replace the previous one, and also save this case, to for future
						else if(existingMatchedMinDistance.containsKey(existingNode.relatedLabelInOtherFrame) && 
								existingMatchedMinDistance.get(existingNode.relatedLabelInOtherFrame) >= existingNode.distanceToRelated)
						{
							System.out.println("REPLACED: " + backupExistingNodeShouldOverwrite);
							matchedGraphs_onF2.removeVertex(backupExistingNodeShouldOverwrite);
							
							matchedGraphs_onF2.addVertex(backupNewGraphNodeShouldOverwrite);
							System.out.println("ADDED: " + backupNewGraphNodeShouldOverwrite);
							associatedNodePairsMatched++;
							
							backupExistingNodeShouldOverwrite = existingNode;
							backupNewGraphNodeShouldOverwrite = new GraphNode(nodeVector_F2.getVector3D(), labelNum_F1, (int) nodeVector_F2.value, distanceBetweenNodes);
						}
						// in the default case never replace the first time if an existing label is detected
						else
						{
							existingMatchedMinDistance.put(existingNode.relatedLabelInOtherFrame,  existingNode.distanceToRelated);
							
							backupExistingNodeShouldOverwrite = existingNode;
							backupNewGraphNodeShouldOverwrite = new GraphNode(nodeVector_F2.getVector3D(), labelNum_F1, (int) nodeVector_F2.value, distanceBetweenNodes);
							
						}
					} 
//					else {
//						 System.out.println("WON'T ADD label F1 " + (labelNum_F1 + 1) + " with distance " 
//								 + distanceBetweenNodes + " since existing is closer " +existingNode);
//					}
				}

				if (shouldAdd) {
					GraphNode newGraphNode = new GraphNode(nodeVector_F2.getVector3D(), labelNum_F1, (int) nodeVector_F2.value, distanceBetweenNodes);
					matchedGraphs_onF2.addVertex(newGraphNode);
					System.out.println("ADDED: " + newGraphNode);
					associatedNodePairsMatched++;
				}

			}

			// if there were associatedLabels but now there are none due to distance
			// threshold removing them then mark this in some way to respond appropriately
			// to "early remove" potential events.
			if (associatedNodePairsMatched == 0) {
				// TODO: This label has no associated fusion/fission event,
				// and therefore might depolarise or nothing
				System.out.println("\tERROR NO MATCHING for Label F1 " + (labelNum_F1 + 1));
				
				matchedGraphs_onF2.removeVertex(backupExistingNode);
				matchedGraphs_onF2.addVertex(backupNewGraphNode);
				System.out.println("\t\tADDED BACKUP: " + backupNewGraphNode);
				associatedNodePairsMatched = 1;
				
			}
		}

		//// DEBUG: Find all potential duplicates which were added (which shouldn't be
		//// the case if the above code did its job)
		for (GraphNode nodeA : matchedGraphs_onF2.vertexSet()) {
			for (GraphNode nodeB : matchedGraphs_onF2.vertexSet()) {
				if (nodeA.sameLocation(nodeB) && (nodeA.distanceToRelated == -1 || nodeB.distanceToRelated == -1) && nodeA.distanceToRelated != nodeB.distanceToRelated) {
					System.out.println("MATCHED " + nodeA + " and " + nodeB);
				}
			}
		}
		////
//		
//		for (GraphNode node : matchedGraphs_onF2.vertexSet()) {
//			GraphNode existingNode = findDuplicateGraphNode(matchedGraphs_onF2, node);
//			if (existingNode != null && node.distanceToRelated == existingNode.distanceToRelated) {
//				System.out.println("FOUND DUPLICATE " + node);
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

				if (diff.x <= 1 && diff.y <= 1 && diff.z <= 1) {
					if (!inputGraph.containsEdge(node1, node2)) {
						inputGraph.addEdge(node1, node2);
						// System.out.println("ADDED EDGE: from " + node1 + " to " + node2);
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
				ih.setPixel((int) v.x, (int) v.y, (int) v.z, node.relatedLabelInOtherFrame + 1); // +1 since the skeleton does not include background
			else if (!showMatched && node.distanceToRelated == -1)
				ih.setPixel((int) v.x, (int) v.y, (int) v.z, Integer.MAX_VALUE); // +1 since the skeleton does not include background

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

	// run along the graph until I find a transition between two labels and that is presumably the event location.
	public List<Vector3D> findEvents(Graph<GraphNode, DefaultEdge> inputGraph, List<Graph<Vector3D, DefaultEdge>> labelsSkeletonGraphs, boolean removeDuplicates, float duplicateDistance) {
		long startTime = System.currentTimeMillis();

		List<Vector3D> eventList = new ArrayList<Vector3D>();
		Map<Integer, Vector3D> alreadyMatched = new HashMap<Integer, Vector3D>();

		for (GraphNode node1 : inputGraph.vertexSet()) {
			for (GraphNode node2 : inputGraph.vertexSet()) {
				Vector3D diff = new Vector3D(Math.abs(node1.location.x - node2.location.x), Math.abs(node1.location.y - node2.location.y), Math.abs(node1.location.z - node2.location.z));

				// the two voxels are neighbouring
				if (diff.x <= 1 && diff.y <= 1 && diff.z <= 1) {
					//.. and not the same voxel, and is at a transition point
					if (node1.relatedLabelInOtherFrame != node2.relatedLabelInOtherFrame) {
						Vector3D eventLocation = getMinHalfwayPoint(node1.location, labelsSkeletonGraphs.get(node1.relatedLabelInOtherFrame)); // (node1.location.add(node2.location).multiply(0.5));
						eventLocation.x = Math.round(eventLocation.x);
						eventLocation.y = Math.round(eventLocation.y);
						eventLocation.z = Math.round(eventLocation.z);

						if (!eventList.contains(eventLocation)) {
							// check both conditions, since I don't know which one is duplicated...
							if(alreadyMatched.containsKey(node1.relatedLabelInOtherFrame))
							{
								Vector3D storedVector = alreadyMatched.get(node1.relatedLabelInOtherFrame);
								eventList.remove(storedVector);
								eventLocation = eventLocation.add(storedVector).multiply(0.5); // take running average to create new event location
								System.out.println("\t UPDATED EVENT LOCATION " + (node1.relatedLabelInOtherFrame+1) + " to " + (node1.relatedLabelInThisFrame +1));
								
							} else if(alreadyMatched.containsKey(node2.relatedLabelInOtherFrame))
							{
								Vector3D storedVector = alreadyMatched.get(node2.relatedLabelInOtherFrame);
								eventList.remove(storedVector);
								eventLocation = eventLocation.add(storedVector).multiply(0.5); // take running average to create new event location
								System.out.println("\t UPDATED EVENT LOCATION2 " + (node2.relatedLabelInOtherFrame+1) + " to " + (node2.relatedLabelInThisFrame +1));
								
							}
								
							alreadyMatched.put(node1.relatedLabelInOtherFrame, eventLocation);
							eventList.add(eventLocation);
							System.out.println("EVENT LOCATION " + eventLocation);
						}
						
						// THIS SECTION IS REPEATED FOR THE ABOVE, but uses node 2's related label. I don't know which one is the appropriate one
						// This methodology can clearly be refined...
						eventLocation = getMinHalfwayPoint(node2.location, labelsSkeletonGraphs.get(node2.relatedLabelInOtherFrame)); // (node1.location.add(node2.location).multiply(0.5));
						eventLocation.x = Math.round(eventLocation.x);
						eventLocation.y = Math.round(eventLocation.y);
						eventLocation.z = Math.round(eventLocation.z);

						if (!eventList.contains(eventLocation)) {
							// check both conditions, since I don't know which one is duplicated...
							if(alreadyMatched.containsKey(node1.relatedLabelInOtherFrame))
							{
								Vector3D storedVector = alreadyMatched.get(node1.relatedLabelInOtherFrame);
								eventList.remove(storedVector);
								eventLocation = eventLocation.add(storedVector).multiply(0.5); // take running average to create new event location
								System.out.println("\t UPDATED EVENT LOCATION " + (node1.relatedLabelInOtherFrame+1) + " to " + (node1.relatedLabelInThisFrame +1));
								
							} else if(alreadyMatched.containsKey(node2.relatedLabelInOtherFrame))
							{
								Vector3D storedVector = alreadyMatched.get(node2.relatedLabelInOtherFrame);
								eventList.remove(storedVector);
								eventLocation = eventLocation.add(storedVector).multiply(0.5); // take running average to create new event location
								System.out.println("\t UPDATED EVENT LOCATION2 " + (node2.relatedLabelInOtherFrame+1) + " to " + (node2.relatedLabelInThisFrame +1));
								
							}
								
							alreadyMatched.put(node1.relatedLabelInOtherFrame, eventLocation);
							eventList.add(eventLocation);
							System.out.println("EVENT LOCATION " + eventLocation);
						}
					}
				}
			}
		}

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
//				System.out.println(i);
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
//					System.out.println("ADDED AVERAGE " + averageEventLocation);
//				}
//			}
		}

		long endTime = System.currentTimeMillis();
		System.out.println("findEvents() - Total execution time: " + (endTime - startTime) + "ms");

		return eventList;
	}
	
	public Vector3D getMinHalfwayPoint(Vector3D point, Graph<Vector3D, DefaultEdge> labeledGraph)
	{
		Vector3D halfwayPoint = null;
		float distance = Float.POSITIVE_INFINITY;
		
		for(Vector3D vertex: labeledGraph.vertexSet())
		{
			if(point.distance(vertex) < distance) {
				distance = (float) point.distance(vertex);
				halfwayPoint = point.add(vertex).multiply(0.5);
			}
		}
		
		return halfwayPoint;
	}
	
	public List<Vector3D> findDepolarisationEvents(List<List<Integer>> associatedLabelsBetweenFrames, ImageInt labeledImage) 
	{
		long startTime = System.currentTimeMillis();

		List<Vector3D> eventList = new ArrayList<Vector3D>();
		
		for(int i = 0; i < associatedLabelsBetweenFrames.size(); ++i)
		{
			if(associatedLabelsBetweenFrames.get(i).size() == 0)
			{			
				Point3D centerPoint = new Object3DVoxels(labeledImage, (i + 1)).getCenterAsPoint();
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
						System.out.println("REMOVED duplicate " + toRemove.get(0));
						int count;
						for (count = 1; count < toRemove.size(); count++) {
							Vector3D removeVertex = toRemove.get(count);
							averageEventLocation = averageEventLocation.add(removeVertex);
							
							eventsGraph.removeVertex(removeVertex);
							System.out.println("REMOVED " + removeVertex);
						}

						averageEventLocation = averageEventLocation.multiply(1.0 / count);
						averageEventLocation.x = Math.round(averageEventLocation.x);
						averageEventLocation.y = Math.round(averageEventLocation.y);
						averageEventLocation.z = Math.round(averageEventLocation.z);

						eventsGraph.addVertex(averageEventLocation);
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
				System.out.println("Succesfully removed " + location);
				break;
			}
		}

		if (initialSize == eventList.size())
			System.out.println("NO DUPLICATES DEMOVED for " + location);

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