package ucla.chou.graphcutshapes;

import ij.IJ;
import ij.ImageStack;
import ij.process.ImageProcessor;

public class NonParametricMixture implements IntensityModel {

    double[] priorparams; // [\bar{\mu}_\Omega, \bar{\mu}_\Delta,a,b ]
    double[] originalpriorparams;
    double[] posteriorparams; // [\mu_\Omega,\mu_\Delta,\tau^2_\Omega,\tau^2_\Delta,
			      // a_\Omega, a_\Delta, b_\Omega, b_\Delta]

    public NonParametricMixture() {
	this.priorparams = new double[8];
	this.priorparams[0]=0;
	this.priorparams[1]=0;
	this.priorparams[2]=0;
	this.priorparams[3]=0;
	this.priorparams[4]=1e-4;
	this.priorparams[5]=1e-4;
	this.priorparams[6]=1e-4; 
	this.priorparams[7]=1e-4;
	

	this.originalpriorparams = this.priorparams.clone();
	//IJ.log("hihi");
    }

    public NonParametricMixture(double[] priorparams) {
	this.priorparams = priorparams;
    }

    public void reset() {
	this.priorparams = this.originalpriorparams.clone();
    }

    

    /*
     * Infer by adding data
     * 
     * @see ucla.chou.graphcutshapes.IntensityModel#Infer(double[][],
     * boolean[][])
     */
    public void Infer(ImageProcessor ip, boolean[][] labels,boolean updatePrior) {
	double sumin = 0;
	double sum2in = 0;
	long areain = 0;

	double sumout = 0;
	double sum2out = 0;
	long areaout = 0;

	int w = ip.getWidth();
	int h = ip.getHeight();

	//IJ.log("WTF IS GOING ON???");
	// Compute sufficient statistics for the data
	for (int y = 0; y < h; y++) {
	    for (int x = 0; x < w; x++) {
		double pixval = ip.getPixelValue(x, y);
		if (labels[x][y]) {
		    sumin += pixval;
		    sum2in += Math.pow(pixval, 2);
		    areain++;
		} else {
		    sumout += pixval;
		    sum2out += Math.pow(pixval, 2);
		    areaout++;
		}
	    }
	}

	double posteriormeanin = (priorparams[0] + sumin) / (1 + areain);
	double posteriormeanout = (priorparams[1] + sumout) / (1 + areaout);

	double posteriorvariancein = (2 * priorparams[2] + sum2in - sumin
		* sumin / areain)
		/ (2 * priorparams[4] + areain - 2)
		+ areain
		/ (2 * priorparams[4] + areain - 2)
		* Math.pow(sumin / areain - priorparams[0], 2) / (1 + areain);
	double posteriorprecisionin = Math.pow(posteriorvariancein, -1);

	double posteriorvarianceout = (2 * priorparams[3] + sum2out - sumout
		* sumout / areaout)
		/ (2 * priorparams[4] + areaout - 2)
		+ areaout
		/ (2 * priorparams[4] + areaout - 2)
		* Math.pow(sumout / areaout - priorparams[1], 2)
		/ (1 + areaout);
	double posteriorprecisionout = Math.pow(posteriorvarianceout, -1);

	double ain = priorparams[4] + areain / 2;
	double aout = priorparams[4] + areaout / 2;
	double bin = priorparams[5] + (sum2in - Math.pow(sumin, 2) / areain)
		/ 2 + areain / 2 * Math.pow(sumin / areain - priorparams[0], 2)
		/ (1 + areain);
	double bout = priorparams[5]
		+ (sum2out - Math.pow(sumout, 2) / areaout) / 2 + areaout / 2
		* Math.pow(sumout / areaout - priorparams[0], 2)
		/ (1 + areaout);
	
	//IJ.log("inner: " + posteriormeanin + "+-" + Math.sqrt(posteriorvariancein) );

	this.posteriorparams = new double[] { posteriormeanin,
		posteriormeanout, posteriorprecisionin, posteriorprecisionout,
		ain, aout, bin, bout };
	if(updatePrior) this.priorparams=this.posteriorparams.clone(); // update the prior with posterior
    }
    
    public double getPosteriorMean(boolean in){
	if(in){
	    return posteriorparams[0];
	}
	return posteriorparams[1];
    }
    
    public double getPosteriorPrecision(boolean in){ 
	if(in){
	    return posteriorparams[2];
	}
	return posteriorparams[3];
    }

    public double pIn(double value) {
	return Math.exp(-logpIn(value));
    }

    public double logpIn(double value) {
	return 0.5*(Math.log(2*Math.PI/getPosteriorPrecision(true))+getPosteriorPrecision(true)*Math.pow(value-getPosteriorMean(true),2));
    }

    public double pOut(double value){
	return Math.exp(-logpOut(value));
    }
    public double logpOut(double value){
	return 0.5*(Math.log(2*Math.PI/getPosteriorPrecision(false))+getPosteriorPrecision(false)*Math.pow(value-getPosteriorMean(false),2));
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
	Infer(is,labels,true);
	
    }

    /**
     * 
     * @param is
     * @param labels  [slice][x][y]
     * @param updatePrior
     */
    public void Infer(ImageStack is, boolean[][][] labels, boolean updatePrior) {
	int slices = is.getSize();
	for(int s=1; s<=slices; s++){
	    IJ.showProgress(s, slices);
	    boolean[][] mask =labels[s-1];
	    Infer(is.getProcessor(s),mask,true);
	}
	
	if(updatePrior) this.priorparams = this.posteriorparams.clone();
	else this.priorparams=this.originalpriorparams.clone();
	
    }

    /**
     * Infer using single image
     */
    public void Infer(ImageProcessor ip, boolean[][] labels) {
	Infer(ip,labels,true);
	
    }

    @Override
    public void multiplyPrecision(double d) {
	posteriorparams[2]*=d;
	posteriorparams[3]*=d;
	
    }

}
