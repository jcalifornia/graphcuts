package ucla.chou.graphcutshapes;

import ij.IJ;
import ij.ImagePlus;
import ij.ImageStack;
import ij.gui.Roi;
import ij.process.ByteProcessor;
import ij.process.ImageProcessor;

/**
 * This class segments a given ImageProcessor. We will interpret the foreground
 * as source, and the background as the sink
 * 
 * TODO conversion to 3D!!!
 * 
 * @author Josh Chang joshchang@ucla.edu
 * 
 */
public class GraphCutSegmenter3D {
    int width;
    int height;
    int depth;
    ImageStack is;
    public GraphCut gc;
    public double lengthPenalty = 20;
    public float Energy;
    private IntensityModel intensityModel;
    private ShapePrior3D shapeKernelDensityEstimate;

    public GraphCutSegmenter3D(ImageStack is, double MU) {
	this(is);
	this.lengthPenalty = MU;

    }

    public GraphCutSegmenter3D(ImageStack is) {
	this.is = is;
	this.width = is.getWidth();
	this.height = is.getHeight();
	this.depth = is.getSize();
	/**
	 * Initialize eight-neighbor graph cut neighborhood should be
	 * 4*(width-1)*(height-1)+width+height-2
	 */
	this.gc = new GraphCut(width * height * depth, 4 * (width - 1)
		* (height - 1) * (depth - 1) + width + height + depth - 1);
    }

    /**
     * Return an Roi corresponding to the current graph cut This function uses
     * the LevelSet class for laziness, should write a method for this guy but
     * whatever
     * 
     * @param slice
     *            0...depth-1
     * 
     * @return Roi corresponding to the boundary
     */
    // TODO: Don't be so lazy
    public Roi returnRoi(int slice) {
	boolean[][] mask = new boolean[width][height];
	for (int y = 0; y < height; y++) {
	    for (int x = 0; x < width; x++) {
		mask[x][y] = gc.getTerminal(slice * width * height + y * width
			+ x) == Terminal.FOREGROUND ? true : false;
	    }
	}
	ImplicitShape2D rp = new ImplicitShape2D(mask);

	return rp.getRoi(false);

    }

    public float fractionInner() {
	int inner = 0;
	for (int slice = 0; slice < depth; slice++) {
	    for (int y = 0; y < height; y++) {
		for (int x = 0; x < width; x++) {
		    if (gc.getTerminal(slice * width * height + y * width + x) == Terminal.FOREGROUND)
			inner++;
		}
	    }
	}
	return (float) (1.0 * inner / width / height);
    }

    public boolean[][][] returnMask() {
	boolean[][][] mask = new boolean[depth][width][height];
	for (int s = 0; s < depth; s++) {
	    IJ.showStatus("Retrieving segmentation mask");
	    for (int y = 0; y < height; y++) {
		for (int x = 0; x < width; x++) {

		    mask[s][x][y] = gc.getTerminal(s * width * height + y
			    * width + x) == Terminal.FOREGROUND ? true : false;

		}
	    }
	    IJ.showProgress(s + 1, depth);
	}
	IJ.showStatus("Done retrieving segmentation mask");
	return mask;
    }

    public ImplicitShape3D getLevelSet() {
	return new ImplicitShape3D(returnMask());
    }

    /**
     * 
     * @param slice
     *            0...depth-1
     * @return
     */
    public ImageProcessor returnMaskProcessor(int slice) {
	int[][] mask = new int[width][height];
	ByteProcessor ip1 = new ByteProcessor(width, height, new byte[width
		* height], null);

	for (int y = 0; y < height; y++) {
	    for (int x = 0; x < width; x++) {
		int nodeid = slice * width * height + y * width + x;
		// IJ.log("getting" + nodeid);
		mask[x][y] = gc.getTerminal(nodeid) == Terminal.FOREGROUND ? 0
			: 255;
		ip1.putPixel(x, y, mask[x][y]);
	    }
	}
	//IJ.log("ok");
	return (ip1);
    }

    public void setLengthPenalty(double MU) {
	this.lengthPenalty = MU;
    }

    /**
     * Set edge weights with no prior information other than length penalty
     * could just set a pointer maybe
     */
    public void setEdgeWeights() {
	for (int slice = 0; slice < depth; slice++) {
	    // IJ.log("alive" + gc.getEdgeNum());
	    IJ.showStatus("Setting n-link weights");
	    for (int y = 0; y < height - 1; y++) {
		for (int x = 0; x < width - 1; x++) {
		    gc.setEdgeWeight(slice * width * height + y * width + x,
			    slice * width * height + y * width + x + 1,
			    (float) (lengthPenalty * Math.PI / 8));
		    gc.setEdgeWeight(slice * width * height + y * width + x,
			    slice * width * height + (y + 1) * width + x,
			    (float) (lengthPenalty * Math.PI / 8));

		    /*
		     * These are for 8-neighbor gc.setEdgeWeight(y * width + x,
		     * (y + 1) * width + x + 1, (float) (lengthPenalty * Math.PI
		     * / 8 / Math.sqrt(2))); gc.setEdgeWeight(y * width + x + 1,
		     * (y + 1) * width + x, (float) (lengthPenalty * Math.PI / 8
		     * / Math.sqrt(2)));
		     */
		    if (slice < depth - 1)
			gc.setEdgeWeight(
				slice * width * height + y * width + x,
				(slice + 1) * width * height + y * width + x,
				(float) (lengthPenalty * Math.PI / 8));

		}
		gc.setEdgeWeight(
			slice * width * height + y * width + width - 1, slice
				* width * height + (y + 1) * width + width - 1,
			(float) (lengthPenalty * Math.PI / 8));
		if (slice < depth - 1)
		    gc.setEdgeWeight(slice * width * height + y * width + width
			    - 1, (slice + 1) * width * height + y * width
			    + width - 1, (float) (lengthPenalty * Math.PI / 8));
	    }
	    for (int x = 0; x < width - 1; x++) {
		gc.setEdgeWeight(slice * width * height + (height - 1) * width
			+ x, slice * width * height + (height - 1) * width + x
			+ 1, (float) (lengthPenalty * Math.PI / 8));
		if (slice < depth - 1)
		    gc.setEdgeWeight(slice * width * height + (height - 1)
			    * width + x, (slice + 1) * width * height
			    + (height - 1) * width + x, (float) (lengthPenalty
			    * Math.PI / 8));
	    }

	    IJ.showProgress(slice + 1, depth);

	}
    }

    /**
     * The prior shape is a reference shape, or rather a prior guess at the true
     * segmentation that is used for determining weighting according to the MM
     * algorithm. Setting these weights should be parallelizable
     * 
     * @param skde
     *            Shape kernel density estimation
     * @param priorShape
     *            The segmentation guess
     */
    public void setEdgeWeights(ShapePrior skde, ImplicitShape2D priorShape) {
	// what I would like to do is reset the shape weights rather
	// than recompute and reset edge weights
	gc.resetEdgeNum();
	float weight;
	/**
	 * Calculate contributions from the density estimation
	 */

	double[] kernelWeights = new double[skde.shapes.size()];

	int j = 0;
	for (; j < skde.shapes.size(); j++) {
	    kernelWeights[j] = Math
		    .exp(-0.5
			    * Math.log(2 * Math.PI / skde.getTauSquared(true))
			    - 0.5
			    * (skde.computeDistance(priorShape,
				    skde.shapes.get(j)) * skde.getTauSquared(true)))
		    * skde.weights[j];

	}
	kernelWeights = normalize(kernelWeights);
	ImplicitShape2D currentshape;
	// weights are length penalty + shape penalty
	float[] weights = new float[] { 0, 0, 0, 0 };

	for (int y = 0; y < height - 1; y++) {
	    for (int x = 0; x < width - 1; x++) {
		weights = new float[] { (float) (lengthPenalty * Math.PI / 8),
			(float) (lengthPenalty * Math.PI / 8),
			(float) (lengthPenalty * Math.PI / 8 / Math.sqrt(2)),
			(float) (lengthPenalty * Math.PI / 8 / Math.sqrt(2)) };

		for (j = 0; j < skde.shapes.size(); j++) {
		    currentshape = skde.shapes.get(j);
		    weights[0] += 0.5
			    * kernelWeights[j]
			    * Math.pow(
				    Math.abs(0.5 * currentshape.get(x, y) + 0.5
					    * currentshape.get(x + 1, y)),
				    skde.shapes.get(0).lambda)
			    * skde.getTauSquared(true);
		    weights[1] += 0.5
			    * kernelWeights[j]
			    * Math.pow(
				    Math.abs(0.5 * currentshape.get(x, y) + 0.5
					    * currentshape.get(x, y + 1))
					    * Math.abs(0.5
						    * currentshape.get(x, y)
						    + 0.5
						    * currentshape
							    .get(x, y + 1)),
				    skde.shapes.get(0).lambda)
			    * skde.getTauSquared(true);
		    weights[2] += kernelWeights[j]
			    * Math.pow(
				    Math.abs(0.5 * currentshape.get(x, y) + 0.5
					    * currentshape.get(x + 1, y + 1))
					    * Math.abs(0.5
						    * currentshape.get(x, y)
						    + 0.5
						    * currentshape
							    .get(x, y + 1)),
				    skde.shapes.get(0).lambda)
			    * skde.getTauSquared(true) / Math.sqrt(2);
		    weights[3] += kernelWeights[j]
			    * Math.pow(
				    Math.abs(0.5 * currentshape.get(x + 1, y)
					    + 0.5 * currentshape.get(x, y + 1))
					    * Math.abs(0.5
						    * currentshape.get(x, y)
						    + 0.5
						    * currentshape
							    .get(x, y + 1)),
				    skde.shapes.get(0).lambda)
			    * skde.getTauSquared(true) / Math.sqrt(2);
		}
		gc.setEdgeWeight(y * width + x, y * width + x + 1, weights[0]);
		gc.setEdgeWeight(y * width + x, (y + 1) * width + x, weights[1]);
		gc.setEdgeWeight(y * width + x, (y + 1) * width + x + 1,
			weights[2]);
		gc.setEdgeWeight(y * width + x + 1, (y + 1) * width + x,
			weights[3]);

	    }

	    weight = (float) (lengthPenalty * Math.PI / 8);
	    for (j = 0; j < skde.shapes.size(); j++) {

		weight += 0.5
			* kernelWeights[j]
			* Math.pow(
				Math.abs(0.5
					* skde.shapes.get(j).get(width - 1, y)
					+ 0.5
					* skde.shapes.get(j).get(width - 1,
						y + 1)),
				skde.shapes.get(0).lambda)
			* skde.getTauSquared(true);

	    }
	    gc.setEdgeWeight(y * width + width - 1,
		    (y + 1) * width + width - 1, weight);
	}
	for (int x = 0; x < width - 1; x++) {
	    weight = (float) (lengthPenalty * Math.PI / 8);
	    for (j = 0; j < skde.shapes.size(); j++) {

		weight += 0.5
			* kernelWeights[j]
			* Math.pow(
				Math.abs(0.5
					* skde.shapes.get(j).get(x, height - 1)
					+ 0.5
					* skde.shapes.get(j).get(x + 1,
						height - 1)),
				skde.shapes.get(0).lambda)
			* skde.getTauSquared(true);

	    }
	    gc.setEdgeWeight((height - 1) * width + x, (height - 1) * width + x
		    + 1, weight);
	}
    }

    /**
     * Set the graph cut weights
     * 
     * @param prior
     *            The prior probability of being in foreground
     * @param intensityModel
     *            The image statistics
     * @param skde 
     *            Shape kernel density estimate
     */

    public void setNodeWeights(ShapePrior3D skde, ImplicitShape3D priorShape,
	    IntensityModel s) {

	double[] kernelWeights = new double[skde.shapes.size()];

	/**
	 * Compute weights
	 */
	for (int j = 0; j < skde.shapes.size(); j++) {
	    kernelWeights[j] = skde.weights[j]
		    * Math.exp(-0.5
			    * Math.log(2 * Math.PI / skde.getTauSquared())
			    - 0.5
			    * (skde.computeDistance(priorShape,
				    skde.shapes.get(j)) * skde.getTauSquared()));

	    if (Double.isNaN(kernelWeights[j])) {
		kernelWeights[j] = 0;
	    }
	}

	final double[] weights = normalize(kernelWeights);

	float source, sink;

	/**
	 * 
	 * Add shape edge weights
	 * 
	 * 
	 * 
	 ****/
	for (int slice = 0; slice < depth; slice++) {
	    for (int x = 0; x < width; x++) {
		for (int y = 0; y < height; y++) {

		    source = 0;
		    sink = 0;

		    if (Double.isNaN(source) || Double.isNaN(sink))
			IJ.log("Kernel object tau^2: " + skde.getTauSquared());

		    for (int j = 0; j < skde.shapes.size(); j++) {
			if (skde.shapes.get(j).get(x, y, slice) < 0.0) {
			    source += weights[j]
				    * Math.pow(
					    Math.abs(skde.shapes.get(j).get(x,
						    y, slice)),
					    skde.shapes.get(0).lambda) / 2
				    * skde.getTauSquared();
			} else
			    sink += weights[j]
				    * Math.pow(
					    Math.abs(skde.shapes.get(j).get(x,
						    y, slice)),
					    skde.shapes.get(0).lambda) / 2
				    * skde.getTauSquared();
		    }

		    // Near image edges, adjust

		    if (x == 0 || y == 0 || x == width - 1 || y == height - 1) {
			sink += this.lengthPenalty * Math.PI
				* (1 + 2 * Math.pow(2, -0.5)) / 8;
			for (int j = 0; j < skde.shapes.size(); j++) {
			    sink += weights[j]
				    * Math.pow(
					    Math.abs(skde.shapes.get(j).get(x,
						    y, slice) + 0.5),
					    skde.shapes.get(0).lambda) / 2
				    * skde.getTauSquared();

			}
		    }

		    if (Double.isNaN(source) || Double.isNaN(sink)) {
			source = 0;
			sink = 0;
		    }
		    gc.setTerminalWeights(y * width + x, source, sink);
		}
	    }
	}

    }

    public double isGetSafeValue(int x, int y, int slice) {
	x = x < 0 ? 0 : x;
	x = x >= width ? width - 1 : x;
	y = y < 0 ? 0 : y;
	y = y >= height ? height - 1 : y;
	return is.getProcessor(slice).getPixelValue(x, y);

    }

    /**
     * Chan - Vese node weights
     */
    public void setNodeWeights(IntensityModel likelihood) {

	float source, sink;
	for (int s = 1; s <= this.depth; s++) {
	    ImageProcessor currentIP = is.getProcessor(s).duplicate();
	    IJ.showStatus("Setting terminal links");
	    // IJ.log("alive" +s);
	    IJ.showProgress(s, this.depth);
	    for (int y = 0; y < height; y++) {
		for (int x = 0; x < width; x++) {

		    source = 0;
		    sink = 0;

		    double pixelval = currentIP.getPixelValue(x, y);
		    if (Double.isNaN(source) || Double.isNaN(sink))
			IJ.showMessage("BLAH1");
		    source += likelihood.logpOut(pixelval);
		    sink += likelihood.logpIn(pixelval);

		    if (Double.isNaN(source) || Double.isNaN(sink)) {
			source = 0;
			sink = 0;

		    }

		    gc.setTerminalWeights((s - 1) * width * height + y * width
			    + x, source, sink);

		}
	    }
	}

    }

    public float relaxEnergy() {
	this.Energy = gc.computeMaximumFlow(true, null);
	return (this.Energy);
    }

    public double[] normalize(double[] numbers) {
	if (numbers.length == 1)
	    return new double[] { 1 };
	double sum = 0;
	double[] num = new double[numbers.length];
	for (int i = 0; i < numbers.length; i++) {
	    sum += numbers[i];
	}
	for (int i = 0; i < numbers.length; i++) {
	    num[i] = sum == 0 ? (double) 1 / numbers.length : numbers[i] / sum;
	}
	return num;
    }

    public void setIntensityModel(IntensityModel s) {
	this.intensityModel = s;
    }

    public IntensityModel getS() {
	return intensityModel;
    }

    public void setShapeKernelDensityEstimate(ShapePrior3D shapep) {
	this.shapeKernelDensityEstimate = shapep;
    }

    public ShapePrior3D getShapeKernelDensityEstimate() {
	return shapeKernelDensityEstimate;
    }

    public ImagePlus returnMaskedImage() {
	ImageStack resultStack = new ImageStack(width, height, null);
	IJ.log("Retrieving segmentation image");
	for (int s = 0; s < depth; s++) {
	    IJ.showStatus("Retrieving segmentation image...");

	    resultStack.addSlice(this.returnMaskProcessor(s).duplicate());
	    IJ.showProgress(s + 1, depth);
	}
	return new ImagePlus("", resultStack);
    }

    public double getVolume() {
	int in = 0;
	int vol = 0;
	for (int s = 0; s < depth; s++) {
	    for (int x = 0; x < width; x++) {
		for (int y = 0; y < height; y++) {
		    vol++;
		    in += gc.getTerminal(s * width * height + y * width + x) == Terminal.FOREGROUND ? 1
			    : 0;
		}
	    }
	}
	return (double) in / vol;
    }

}
