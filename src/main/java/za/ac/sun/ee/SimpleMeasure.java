/* Modified by Rensu Theart - 2021
 * 
 * 
 * Use can use the SimpleMeasure like this:
 * Print table of all parameters:
 * 	
 * 		labels_F1_measure.printStats(labels_F1_measure.compactnessStats);
 * 
 * Get a sub element in the table:
 * 
 * 		System.out.println("Label at 0: " + labels_F1_measure.compactnessStats.get(0)[SimpleMeasure.Compactness.LABEL.ordinal()]);
 */

package za.ac.sun.ee;

import ij.IJ;
import ij.ImagePlus;

import java.lang.reflect.InvocationTargetException;
import java.util.Arrays;
import java.util.Collections;
import java.util.Comparator;
import java.util.LinkedList;
import java.util.List;

import ij.plugin.Duplicator;
import mcib3d.geom.Point3D;
import mcib3d.geom2.Objects3DIntPopulation;
import mcib3d.geom2.measurements.*;
import mcib3d.image3d.ImageHandler;

/**
 *
 **
 * /** Copyright (C) 2008- 2012 Thomas Boudier and others
 *
 *
 *
 * This file is part of mcib3d
 *
 * mcib3d is free software; you can redistribute it and/or modify it under the
 * terms of the GNU General Public License as published by the Free Software
 * Foundation; either version 3 of the License, or (at your option) any later
 * version.
 *
 * This program is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
 * details.
 *
 * You should have received a copy of the GNU General Public License along with
 * this program. If not, see <http://www.gnu.org/licenses/>.
 *
 * @author thomas
 */
public class SimpleMeasure {
	private Objects3DIntPopulation population;

	public enum Centroid {
		LABEL, CX_PIX, CY_PIX, CZ_PIX, CX_UNIT, CY_UNIT, CZ_UNIT
	};

	public enum Volume {
		LABEL, VOLUME_PIX, VOLUME_UNIT
	};

	public enum Surface {
		LABEL, SURFACE_PIX, SURFACE_UNIT, SURFACE_CORRECTED, SURFACE_NB_VOXELS
	};

	public enum Compactness {
		LABEL, COMP_UNIT, COMP_PIX, COMP_CORRECTED, COMP_DISCRETE, SPHER_PIX, SPHER_UNIT, SPHER_CORRECTED, SPHER_DISCRETE
	};

	public enum Elipsoid {
		LABEL, ELL_VOL_UNIT, ELL_SPARENESS, ELL_MAJOR_RADIUS_UNIT, ELL_ELONGATION, ELL_FLATNESS
	};

	public enum DistanceCentreContour {
		LABEL, DIST_CENTER_MIN_PIX, DIST_CENTER_MAX_PIX, DIST_CENTER_AVG_PIX, DIST_CENTER_SD_PIX, DIST_CENTER_MIN_UNIT, DIST_CENTER_MAX_UNIT, DIST_CENTER_AVG_UNIT, DIST_CENTER_SD_UNIT
	};

	public enum Feret {
		LABEL, FERET_UNIT, FERET1_X, FERET1_Y, FERET1_Z, FERET2_X, FERET2_Y, FERET2_Z
	};

	/**
	 * Based on class MeasureCentroid
	 * 
	 * Indices: 0 LABEL, 1 CX_PIX, 2 CY_PIX, 3 CZ_PIX, 4 CX_UNIT, 5 CY_UNIT, 6
	 * CZ_UNIT
	 * 
	 */
	public List<Double[]> centroidStats;
	/**
	 * Based on class MeasureVolume
	 * 
	 * Indices: 0 LABEL, 1 VOLUME_PIX, 2 VOLUME_UNIT
	 */
	public List<Double[]> volumeStats;
	/**
	 * Based on class MeasureSurface
	 * 
	 * Indices: 0 LABEL, 1 SURFACE_PIX, 2 SURFACE_UNIT, 3 SURFACE_CORRECTED, 4
	 * SURFACE_NB_VOXELS
	 */
	public List<Double[]> surfaceStats;
	public List<Double[]> intensityStats;
	/**
	 * Based on class MeasureCompactness
	 * 
	 * Indices: 0 LABEL, 1 COMP_UNIT, 2 COMP_PIX, 3 COMP_CORRECTED, 4 COMP_DISCRETE,
	 * 5 SPHER_PIX, 6 SPHER_UNIT, 7 SPHER_CORRECTED, 8 SPHER_DISCRETE
	 */
	public List<Double[]> compactnessStats;
	/**
	 * Based on class MeasureEllipsoid
	 * 
	 * Indices: 0 LABEL, 1 ELL_VOL_UNIT, 2 ELL_SPARENESS, 3 ELL_MAJOR_RADIUS_UNIT, 4
	 * ELL_ELONGATION, 5 ELL_FLATNESS
	 */
	public List<Double[]> elipsoidStats;
	/**
	 * Based on class MeasureDistancesCenter
	 * 
	 * Indices: 0 LABEL, 1 DIST_CENTER_MIN_PIX, 2 DIST_CENTER_MAX_PIX, 3
	 * DIST_CENTER_AVG_PIX, 4 DIST_CENTER_SD_PIX, 5 DIST_CENTER_MIN_UNIT, 6
	 * DIST_CENTER_MAX_UNIT, 7 DIST_CENTER_AVG_UNIT, 8 DIST_CENTER_SD_UNIT
	 */
	public List<Double[]> distanceCentreContourStats;
	/**
	 * Based on class MeasureFeret
	 * 
	 * Indices: 0 LABEL, 1 FERET_UNIT, 2 FERET1_X, 3 FERET1_Y, 4 FERET1_Z, 5
	 * FERET2_X, 6 FERET2_Y, 7 FERET2_Z
	 */
	public List<Double[]> feretStats;

	/**
	 * The number of structures for which parameters exist
	 */
	public int length = 0;

	public SimpleMeasure(ImagePlus in) {
		population = new Objects3DIntPopulation(ImageHandler.wrap(in));
		calculateAll(in);
	}

	public SimpleMeasure(ImageHandler in) {
		population = new Objects3DIntPopulation(in);
		calculateAll(in.getImagePlus());
	}

	private void calculateAll(ImagePlus in) {
		centroidStats = getMeasuresCentroid();
		volumeStats = getMeasuresVolume();
		surfaceStats = getMeasuresSurface();
		intensityStats = getMeasuresStats(in);
		compactnessStats = getMeasuresCompactness();
		elipsoidStats = getMeasuresEllipsoid();
		distanceCentreContourStats = getMeasuresDistanceCentreContour();
		feretStats = getMeasuresFeret();

		length = centroidStats.size();
	}

	private List<Double[]> sortList(List<Double[]> inList) {
		Collections.sort(inList, new Comparator<Double[]>() {
			public int compare(Double[] doubles, Double[] otherDoubles) {
				return doubles[0].compareTo(otherDoubles[0]);
			}
		});

		return inList;
	}

	public void printStats(List<Double[]> stats) {
		for (Double[] da : stats) {
			System.out.println(Arrays.toString(da));
		}
	}

	private List<Double[]> getMeasuresCentroid() {
		List<Double[]> results = new LinkedList<>();
		try {
			return sortList(population.getMeasurements(MeasureCentroid.class));
		} catch (NoSuchMethodException e) {
			e.printStackTrace();
		} catch (InvocationTargetException e) {
			e.printStackTrace();
		} catch (InstantiationException e) {
			e.printStackTrace();
		} catch (IllegalAccessException e) {
			e.printStackTrace();
		}

		return results;
	}

	private List<Double[]> getMeasuresVolume() {
		List<Double[]> results = new LinkedList<>();
		try {
			return sortList(population.getMeasurements(MeasureVolume.class));
		} catch (NoSuchMethodException e) {
			e.printStackTrace();
		} catch (InvocationTargetException e) {
			e.printStackTrace();
		} catch (InstantiationException e) {
			e.printStackTrace();
		} catch (IllegalAccessException e) {
			e.printStackTrace();
		}
		return results;
	}

	private List<Double[]> getMeasuresSurface() {
		List<Double[]> results = new LinkedList<>();
		try {
			return sortList(population.getMeasurements(MeasureSurface.class));
		} catch (NoSuchMethodException e) {
			e.printStackTrace();
		} catch (InvocationTargetException e) {
			e.printStackTrace();
		} catch (InstantiationException e) {
			e.printStackTrace();
		} catch (IllegalAccessException e) {
			e.printStackTrace();
		}
		return results;
	}

	private List<Double[]> getMeasuresStats(ImageHandler raw) {
		return sortList(population.getMeasurementsIntensity(raw));
	}

	private List<Double[]> getMeasuresStats(ImagePlus myPlus) {
		return getMeasuresStats(ImageHandler.wrap(myPlus));
	}

	private List<Double[]> getMeasuresCompactness() {
		List<Double[]> results = new LinkedList<>();
		try {
			results.addAll(population.getMeasurements(MeasureCompactness.class));
		} catch (NoSuchMethodException e) {
			e.printStackTrace();
		} catch (InvocationTargetException e) {
			e.printStackTrace();
		} catch (InstantiationException e) {
			e.printStackTrace();
		} catch (IllegalAccessException e) {
			e.printStackTrace();
		}
		return sortList(results);
	}

	private List<Double[]> getMeasuresEllipsoid() {
		List<Double[]> results = new LinkedList<>();
		try {
			results.addAll(population.getMeasurements(MeasureEllipsoid.class));
		} catch (NoSuchMethodException e) {
			e.printStackTrace();
		} catch (InvocationTargetException e) {
			e.printStackTrace();
		} catch (InstantiationException e) {
			e.printStackTrace();
		} catch (IllegalAccessException e) {
			e.printStackTrace();
		}
		return sortList(results);
	}

	private List<Double[]> getMeasuresDistanceCentreContour() {
		List<Double[]> results = new LinkedList<>();
		try {
			results.addAll(population.getMeasurements(MeasureDistancesCenter.class));
		} catch (NoSuchMethodException e) {
			e.printStackTrace();
		} catch (InvocationTargetException e) {
			e.printStackTrace();
		} catch (InstantiationException e) {
			e.printStackTrace();
		} catch (IllegalAccessException e) {
			e.printStackTrace();
		}
		return sortList(results);
	}

	private List<Double[]> getMeasuresFeret() {
		List<Double[]> results = new LinkedList<>();
		try {
			results.addAll(population.getMeasurements(MeasureFeret.class));
		} catch (NoSuchMethodException e) {
			e.printStackTrace();
		} catch (InvocationTargetException e) {
			e.printStackTrace();
		} catch (InstantiationException e) {
			e.printStackTrace();
		} catch (IllegalAccessException e) {
			e.printStackTrace();
		}
		return sortList(results);
	}

	public static ImagePlus extractCurrentStack(ImagePlus plus) {
		// check dimensions
		int[] dims = plus.getDimensions();// XYCZT
		int channel = plus.getChannel();
		int frame = plus.getFrame();
		ImagePlus stack;
		// crop actual frame
		if ((dims[2] > 1) || (dims[4] > 1)) {
			IJ.log("Hyperstack found, extracting current channel " + channel + " and frame " + frame);
			Duplicator duplicator = new Duplicator();
			stack = duplicator.run(plus, channel, channel, 1, dims[3], frame, frame);
		} else
			stack = plus.duplicate();

		return stack;
	}

	public Point3D[] getCentroidPoints() {
		Point3D[] pointArray = new Point3D[length];

		for (int i = 0; i < length; i++) {
			pointArray[i] = new Point3D(centroidStats.get(i)[1], centroidStats.get(i)[2], centroidStats.get(i)[3]);
		}

		return pointArray;
	}

	public int[] getVolumes() {
		int[] volumeArray = new int[length];

		for (int i = 0; i < length; i++) {
			volumeArray[i] = volumeStats.get(i)[1].intValue();
		}

		return volumeArray;
	}
}
