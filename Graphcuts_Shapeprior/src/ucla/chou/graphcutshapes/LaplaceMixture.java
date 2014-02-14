package ucla.chou.graphcutshapes;

import java.util.AbstractList;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

import ij.IJ;
import ij.ImageStack;
import ij.process.ImageProcessor;

public class LaplaceMixture implements IntensityModel {

    double[] posteriorparams; // [\mu_\Omega,\mu_\Delta,b_\Omega,b^2_\Delta]

    public LaplaceMixture() {

    }

    /*
     * Infer by adding data
     * 
     * @see ucla.chou.graphcutshapes.IntensityModel#Infer(double[][],
     * boolean[][])
     */
    public void Infer(ImageProcessor ip, boolean[][] labels, boolean updatePrior) {

	List<Float> inPixels = new ArrayList<Float>();
	List<Float> outPixels = new ArrayList<Float>();

	int w = ip.getWidth();
	int h = ip.getHeight();

	double medin, medout, bin, bout;

	for (int y = 0; y < h; y++) {
	    for (int x = 0; x < w; x++) {
		float pixval = ip.getPixelValue(x, y);
		if (labels[x][y]) {
		    inPixels.add((float) pixval);
		} else {
		    outPixels.add(pixval);

		}
	    }
	}

	Collections.sort(inPixels);
	Collections.sort(outPixels);

	if (inPixels.size() % 2 == 1)
	    medin = inPixels.get((int) Math.floor(inPixels.size() / 2));
	else
	    medin = 0.5 * inPixels.get(inPixels.size() / 2 - 1) + 0.5
		    * inPixels.get(inPixels.size() / 2);
	if (outPixels.size() % 2 == 1)
	    medout = outPixels.get((int) Math.floor(outPixels.size() / 2));
	else
	    medout = 0.5 * outPixels.get(outPixels.size() / 2 - 1) + 0.5
		    * outPixels.get(outPixels.size() / 2);

	double boutsum = 0;
	double binsum = 0;

	for (float s : inPixels) {
	    binsum += Math.abs(s - medin);
	}

	for (float s : outPixels) {
	    boutsum += Math.abs(s - medout);
	}

	bin = binsum / inPixels.size();
	bout = boutsum / outPixels.size();

	this.posteriorparams = new double[] { medin, medout, 1 / bin, 1 / bout };

    }

    public double getPosteriorMean(boolean in) {
	if (in) {
	    return posteriorparams[0];
	}
	return posteriorparams[1]; 
    }

    public double getPosteriorPrecision(boolean in) {
	if (in) {
	    return posteriorparams[2]/4;
	}
	return posteriorparams[3]/4;
    }

    public double pIn(double value) {
	return Math.exp(-logpIn(value));
    }

    public double logpIn(double value) {
	return -(Math.log(getPosteriorPrecision(true) / 2) - getPosteriorPrecision(true)
		* Math.abs(value - getPosteriorMean(true)));
    }

    public double pOut(double value) {
	return Math.exp(-logpOut(value));
    }

    public double logpOut(double value) {
	return -(Math.log(getPosteriorPrecision(false) / 2) - getPosteriorPrecision(false)
		* Math.abs(value - getPosteriorMean(false)));
    }

    public double tauin() {
	return posteriorparams[2];
    }

    public double tauout() {
	return posteriorparams[3];
    }

    /**
     * Infer using the imageStack
     */
    public void Infer(ImageStack is, boolean[][][] labels) {
	Infer(is, labels, true);

    }

    /**
     * 
     * @param is
     * @param labels
     *            [slice][x][y]
     * @param updatePrior
     */
    public void Infer(ImageStack is, boolean[][][] labels, boolean updatePrior) {
	int slices = is.getSize();
	List<Float> inPixels = new ArrayList<Float>();
	List<Float> outPixels = new ArrayList<Float>();
	double medin, medout, bin, bout;

	for (int s = 1; s <= slices; s+=2) {
	    IJ.showProgress(s, slices);
	    boolean[][] mask = labels[s - 1];

	    for (int y = 0; y < mask[0].length; y+=2) {
		for (int x = 0; x < mask.length; x+=2) {

		    if (mask[x][y])
			inPixels.add(is.getProcessor(s).getPixelValue(x, y));
		    else
			outPixels.add(is.getProcessor(s).getPixelValue(x, y));
		}
	    }

	}

	Collections.sort(inPixels);
	Collections.sort(outPixels);

	if (inPixels.size() % 2 == 1)
	    medin = inPixels.get((int) Math.floor(inPixels.size() / 2));
	else
	    medin = 0.5 * inPixels.get(inPixels.size() / 2 - 1) + 0.5
		    * inPixels.get(inPixels.size() / 2);
	if (outPixels.size() % 2 == 1)
	    medout = outPixels.get((int) Math.floor(outPixels.size() / 2));
	else
	    medout = 0.5 * outPixels.get(outPixels.size() / 2 - 1) + 0.5
		    * outPixels.get(outPixels.size() / 2);

	double boutsum = 0;
	double binsum = 0;

	for (float s : inPixels) {
	    binsum += Math.abs(s - medin);
	}

	for (float s : outPixels) {
	    boutsum += Math.abs(s - medout);
	}

	bin = binsum / inPixels.size();
	bout = boutsum / outPixels.size();

	this.posteriorparams = new double[] { medin, medout, 1 / bin, 1 / bout };

    }

    /**
     * Infer using single image
     */
    public void Infer(ImageProcessor ip, boolean[][] labels) {
	Infer(ip, labels, true);

    }

    private class SortedList<E> extends AbstractList<E> {

	private ArrayList<E> internalList = new ArrayList<E>();

	// Note that add(E e) in AbstractList is calling this one
	@Override
	public void add(int position, E e) {
	    internalList.add(e);
	    Collections.sort(internalList, null);
	}

	@Override
	public E get(int i) {
	    return internalList.get(i);
	}

	@Override
	public int size() {
	    return internalList.size();
	}

    }

    @Override
    public void multiplyPrecision(double d) {
	this.posteriorparams[2]*=d;
	this.posteriorparams[3]*=d;
	
    }

}
