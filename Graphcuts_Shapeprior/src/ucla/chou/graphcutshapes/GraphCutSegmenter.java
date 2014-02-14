package ucla.chou.graphcutshapes;

import ij.IJ;
import ij.ImagePlus;
import ij.gui.Roi;
import ij.process.ByteProcessor;
import ij.process.ImageProcessor;

/**
 * This class segments a given ImageProcessor. We will interpret the foreground
 * as source, and the background as the sink
 * 
 * @author Josh Chang joshchang@ucla.edu
 * 
 */
public class GraphCutSegmenter {
    int width;
    int height;
    ImageProcessor ip;
    public GraphCut gc;
    public float lengthPenalty = 20;
    public float Energy;
    private IntensityModel intensityModel;
    private ShapePrior shapeKernelDensityEstimate;

    public boolean[][] innermask;
    public boolean[][] outermask; 
    private double edgeweight=10;
    final boolean testLinear=true;

    public GraphCutSegmenter(ImageProcessor ip, float mu) {
	this(ip);
	this.lengthPenalty = mu;

    }

    public void setInnerMask(boolean[][] im) {
	this.innermask = im;
    }

    public GraphCutSegmenter(ImageProcessor ip) {
	this.ip = ip;
	this.width = ip.getWidth();
	this.height = ip.getHeight();
	/**
	 * Initialize eight-neighbor graph cut neighborhood should be
	 * 4*(width-1)*(height-1)+width+height-2
	 */
	this.gc = new GraphCut(width * height, 4 * (width - 1) * (height - 1)
		+ width + height - 1);
    }

    /**
     * Return an Roi corresponding to the current graph cut This function uses
     * the LevelSet class for laziness, should write a method for this guy but
     * whatever
     * 
     * @return Roi corresponding to the boundary
     */
    public Roi returnRoi() {
	boolean[][] mask = new boolean[width][height];
	for (int y = 0; y < height; y++) {
	    for (int x = 0; x < width; x++) {
		mask[x][y] = gc.getTerminal(y * width + x) == Terminal.FOREGROUND ? true
			: false;
	    }
	}
	ImplicitShape2D rp = new ImplicitShape2D(mask);

	return rp.getRoi(true);

    }

    public float fractionInner() {
	int inner = 0;
	for (int y = 0; y < height; y++) {
	    for (int x = 0; x < width; x++) {
		if (gc.getTerminal(y * width + x) == Terminal.FOREGROUND)
		    inner++;
	    }
	}
	return (float) (1.0 * inner / width / height);
    } 

    public boolean[][] returnMask() {
	boolean[][] mask = new boolean[width][height];
	for (int y = 0; y < height; y++) {
	    for (int x = 0; x < width; x++) {
		mask[x][y] = gc.getTerminal(y * width + x) == Terminal.FOREGROUND ? true
			: false;
	    }
	}
	return mask;
    }

    public ImplicitShape2D getLevelSet() {
	return new ImplicitShape2D(returnMask());
    }

    public ImageProcessor returnMaskProcessor() {
	int[][] mask = new int[width][height];
	ByteProcessor ip1 = new ByteProcessor(width, height, new byte[width
		* height], null);
	for (int y = 0; y < height; y++) {
	    for (int x = 0; x < width; x++) {
		mask[x][y] = gc.getTerminal(y * width + x) == Terminal.FOREGROUND ? 0
			: 255;
		ip1.putPixel(x, y, mask[x][y]);
	    }
	}
	return (ip1);
    }

    public ImagePlus returnMaskedImage() {

	return (new ImagePlus("", returnMaskProcessor()));
    }

    public void setLengthPenalty(float MU) {
	this.lengthPenalty = MU;
    }

    /**
     * Set edge weights with no prior information other than length penalty
     * could just set a pointer maybe
     */
    public void setEdgeWeights() {
	gc.resetEdgeNum();

	for (int y = 0; y < height - 1; y++) {
	    for (int x = 0; x < width - 1; x++) {
		gc.setEdgeWeight(y * width + x, y * width + x + 1,
			(float) (lengthPenalty * Math.PI / 8));
		gc.setEdgeWeight(y * width + x, (y + 1) * width + x,
			(float) (lengthPenalty * Math.PI / 8));
		gc.setEdgeWeight(y * width + x, (y + 1) * width + x + 1,
			(float) (lengthPenalty * Math.PI / 8 / Math.sqrt(2)));
		gc.setEdgeWeight(y * width + x + 1, (y + 1) * width + x,
			(float) (lengthPenalty * Math.PI / 8 / Math.sqrt(2)));

	    }
	    gc.setEdgeWeight(y * width + width - 1,
		    (y + 1) * width + width - 1, (float) (lengthPenalty
			    * Math.PI / 8));
	}
	for (int x = 0; x < width - 1; x++) {
	    gc.setEdgeWeight((height - 1) * width + x, (height - 1) * width + x
		    + 1, (float) (lengthPenalty * Math.PI / 8));
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

	int j;
	double[] kernelWeights = getShapeWeights(skde, priorShape);
	
	if(testLinear){
	    for( j =0; j<skde.shapes.size();j++){
		kernelWeights[j]=(double) 1.0/skde.shapes.size();
	    }
	}
	
	for(j=0; j<skde.shapes.size();j++){
	    IJ.log("Shape  "+j + " weight " + kernelWeights[j]);
	}
	
	
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
			    * skde.getTauSquared(true) * Math.PI / 8 ;
		    weights[0]+= this.edgeweight* Math.PI / 8 * Math.exp(-Math.pow( (ip.getPixelValue(x, y)-ip.getPixelValue(x+1,y)),2)/(40000));
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
			    * skde.getTauSquared(true) * Math.PI / 8 ;
		    weights[1]+= this.edgeweight* Math.PI / 8* Math.exp(-Math.pow(ip.getPixelValue(x, y)-ip.getPixelValue(x,y+1),2)/40000);
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
			    * skde.getTauSquared(true) / Math.sqrt(2) * Math.PI
			    / 8 ;
		    weights[2]+= this.edgeweight* Math.PI / 8 / Math.sqrt(2)* Math.exp(-Math.pow(ip.getPixelValue(x, y)-ip.getPixelValue(x+1,y+1),2)/40000);
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
			    * skde.getTauSquared(true) / Math.sqrt(2) * Math.PI
			    / 8 ;
		    weights[3]+= this.edgeweight* Math.PI / 8 / Math.sqrt(2)* Math.exp(-Math.pow(ip.getPixelValue(x+1, y)-ip.getPixelValue(x,y+1),2)/40000);
		}
		gc.setEdgeWeight(y * width + x, y * width + x + 1,
			Float.isNaN(weights[0]) ? Float.MAX_VALUE : weights[0]);
		gc.setEdgeWeight(y * width + x, (y + 1) * width + x,
			Float.isNaN(weights[1]) ? Float.MAX_VALUE : weights[1]);
		gc.setEdgeWeight(y * width + x, (y + 1) * width + x + 1,
			Float.isNaN(weights[2]) ? Float.MAX_VALUE : weights[2]);
		gc.setEdgeWeight(y * width + x + 1, (y + 1) * width + x,
			Float.isNaN(weights[3]) ? Float.MAX_VALUE : weights[3]);

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
			* skde.getTauSquared(true) * Math.PI / 8 ;

	    }
	    gc.setEdgeWeight(y * width + width - 1,
		    (y + 1) * width + width - 1,
		    Float.isNaN(weight) ? Float.MAX_VALUE : weight);
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
				skde.shapes.get(j).lambda)
			* skde.getTauSquared(true) * Math.PI / 8 ;

	    }
	    gc.setEdgeWeight((height - 1) * width + x, (height - 1) * width + x
		    + 1, Float.isNaN(weight) ? Float.MAX_VALUE : weight);
	}
    }

    /**
     * Get the shape weights
     * 
     * @param skde
     * @param shape
     * @return
     */
    public double[] getShapeWeights(ShapePrior skde, ImplicitShape2D shape) {
	double[] distances = new double[skde.shapes.size()];
	double[] weights = new double[skde.shapes.size()];
	double mindistance = Double.MAX_VALUE;
	double[] kernelWeights = skde.weights;
	
	for (int j = 0; j < skde.shapes.size(); j++) {
	    // there is some sort of buffer overflow

	    distances[j] = skde.computeDistance(shape, skde.shapes.get(j))
		    * skde.getTauSquared(true);
	
	    /**
	     * If NaN, set to infinite distance
	     */
	    if (Double.isNaN(kernelWeights[j])) {
		distances[j] = Double.MAX_VALUE;
	    }
	}

	IJ.log("min distance " + mindistance);
	
	
	for(int j=0; j<skde.shapes.size();j++){
	    weights[j]=0;
	    for(int k=0;k<skde.shapes.size();k++){
		weights[j]+=kernelWeights[k]*Math.exp(distances[j]-distances[k]);
	    }
	    weights[j]=kernelWeights[j]/weights[j];
	}
	weights = normalize(weights);
	IJ.log("shape weights " + weights[0] + " , " + weights[1] + " , "
		+ weights[2] + " , " + weights[weights.length - 1]);
	IJ.log("shape distances " + distances[0] + " , " + distances[1] + " , "
		+ distances[2] + " , " + distances[distances.length - 1]);
	return weights;
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
 
    public void setNodeWeights(ShapePrior skde, ImplicitShape2D priorShape,
	    IntensityModel likelihood) {

	double[] weights = getShapeWeights(skde, priorShape);
	if(testLinear){
	    for(int  j =0; j<skde.shapes.size();j++){
		weights[j]=1/skde.shapes.size();
	    }
	}
	float source, sink;

	/**
	 * 
	 * Add shape edge weights
	 * 
	 * 
	 * 
	 ****/
	for (int x = 0; x < width; x++) {
	    for (int y = 0; y < height; y++) {

		source = 0;
		sink = 0;

		float pixelval = ip.getPixelValue(x, y);
		if (Float.isNaN(source) || Float.isNaN(sink))
		    IJ.showMessage("NaN error for n-link at (" + x + "," + y
			    + ")");
		source += likelihood.logpOut(pixelval);
		sink += likelihood.logpIn(pixelval);

		if (Float.isNaN(source)) {
		    source = Float.MAX_VALUE;
		}
		if (Float.isNaN(sink)) {
		    sink = Float.MAX_VALUE;
		}

		for (int j = 0; j < skde.shapes.size(); j++) {
		    if (skde.shapes.get(j).in(x, y)) {
			source += weights[j]
				* skde.shapes.get(j).rho(x, y) / 2
				* skde.getTauSquared(true);
		    } else
			sink += weights[j]
				* skde.shapes.get(j).rho(x, y) / 2
				* skde.getTauSquared(true);
		}

		// Near image edges, adjust

		if (x == 0 || y == 0 || x == width - 1 || y == height - 1) {
		    sink += this.lengthPenalty * Math.PI * 2
			    * (1 + 2 * Math.pow(2, -0.5)) / 8;
		    for (int j = 0; j < skde.shapes.size(); j++) {
			sink += weights[j]
				* Math.pow(Math.abs(skde.shapes.get(j)
					.get(x, y) + 0.5),
					skde.shapes.get(0).lambda) / 2
				* skde.getTauSquared(true) * Math.PI * 2
				* (1 + 2 * Math.pow(2, -0.5)) / 8;
			;

		    }
		}

		gc.setTerminalWeights(y * width + x,
			Float.isNaN(source) ? Float.MAX_VALUE : source,
			Float.isNaN(sink) ? Float.MAX_VALUE : sink);
	    }
	}

    }

    public float ipGetSafeValue(int x, int y) {
	x = x < 0 ? 0 : x;
	x = x >= width ? width - 1 : x;
	y = y < 0 ? 0 : y;
	y = y >= height ? height - 1 : y;
	return ip.getPixelValue(x, y);

    }

    /**
     * Chan - Vese node weights
     */
    public void setNodeWeights(IntensityModel likelihood) {

	float source, sink;
	for (int x = 0; x < width; x++) {
	    for (int y = 0; y < height; y++) {
		source = 0;
		sink = 0;

		float pixelval = ip.getPixelValue(x, y);
		source += likelihood.logpOut(pixelval);
		sink += likelihood.logpIn(pixelval);

		gc.setTerminalWeights(y * width + x,
			Float.isNaN(source) ? Float.MAX_VALUE : source,
			Float.isNaN(sink) ? Float.MAX_VALUE : sink);
	    }
	}

    }

    public float relaxEnergy() {
	this.Energy = gc.computeMaximumFlow(true, null);
	return (this.Energy);
    }

    public float overlap(boolean[][] mask1, boolean[][] mask2) {
	int w1 = mask1.length;
	int w2 = mask2.length;
	int h1 = mask1[0].length;
	int h2 = mask2[0].length;
	if (w1 != w2 || h1 != h2)
	    return 0;
	float area = 0;
	int over = 0;
	int in = 0;
	for (int x = 0; x < w1; x++) {
	    for (int y = 0; y < h1; y++) {
		if (!mask1[x][y]) {
		    in++;
		    if (!mask2[x][y]) {
			over++;
		    }

		}
	    }
	}
	area = (float) over / in;
	return area;

    }

    public float[] normalize(float[] kernelWeights) {
	if (kernelWeights.length == 1)
	    return new float[] { 1 };
	float sum = 0;
	float[] num = new float[kernelWeights.length];
	for (int i = 0; i < kernelWeights.length; i++) {
	    sum += kernelWeights[i];
	}
	for (int i = 0; i < kernelWeights.length; i++) {
	    num[i] = sum == 0 ? (float) 1 / kernelWeights.length
		    : kernelWeights[i] / sum;
	}
	return num;
    }

    public double[] normalize(double[] kernelWeights) {
	if (kernelWeights.length == 1)
	    return new double[] { 1 };
	double sum = 0;
	double[] num = new double[kernelWeights.length];
	for (int i = 0; i < kernelWeights.length; i++) {
	    sum += kernelWeights[i];
	}
	for (int i = 0; i < kernelWeights.length; i++) {
	    num[i] = sum == 0 ? (double) 1 / kernelWeights.length
		    : kernelWeights[i] / sum;
	}
	return num;
    }

    public void setIntensityModel(IntensityModel s) {
	this.intensityModel = s;
    }

    public IntensityModel getS() {
	return intensityModel;
    }

    public void setShapeKernelDensityEstimate(
	    ShapePrior shapeKernelDensityEstimate) {
	this.shapeKernelDensityEstimate = shapeKernelDensityEstimate;
    }

    public ShapePrior getShapeKernelDensityEstimate() {
	return shapeKernelDensityEstimate;
    }

}
