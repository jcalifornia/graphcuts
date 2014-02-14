package ucla.chou.graphcutshapes;

import ij.IJ;
import ij.gui.Roi;
import ij.process.ImageProcessor;

import java.util.Vector;

import ucla.chou.graphcutshapes.util.Ellipse;

import Jama.Matrix;

/** This class fits an ellipse to an ROI. */
public class SuperEllipseFitter {

    static final double HALFPI = 1.5707963267949;

    /** Initialized by makeRoi() */
    public Vector<Integer> xCoordinates;
    /** Initialized by makeRoi() */
    public Vector<Integer> yCoordinates;
    /** Initialized by makeRoi() */

    public double xc, yc, a, b, theta, eps;

    /**
     * Weigting for the objective functions
     */
    private double w1 = 0;
    private double w2 = 0;

    private int width;
    private int height;

    private ImageProcessor ip;

    private ImplicitShape2D inputShape;

    private Roi inputRoi;
    private Roi fittedRoi;

    /**
     * Fits an ellipse to the current ROI. The 'stats' argument, currently not
     * used, can be null. The fit parameters are returned in public fields.
     */
    public void fit() {

	// First get starting values from the inertial
	eps = 1;
	double[] params = ellipseFitter();
	xc = params[0];
	yc = params[1];
	a = params[2];
	b = params[3];
	theta = params[4];

//	for (int iter = 0; iter < 10; iter++) {
//
//	    double[][] step = getLMAstep(new double[] { xc, yc, a, b, theta,
//		    eps });
//	    xc += step[0][0];
//	    yc += step[1][0];
//	    IJ.log(step[0][0] + "");
//
//	}

	this.fittedRoi = (new ImplicitShape2D(getMask())).getRoi(false);

    }

    /**
     * Populate the lists of coordinates
     * @param roi
     * @param width
     * @param height
     */
    public SuperEllipseFitter(Roi roi, int width, int height) {
	this.setInputRoi(roi);
	this.inputShape = new ImplicitShape2D(roi, width, height);
	this.width = width;
	this.height = height;
	xCoordinates = new Vector<Integer>();
	yCoordinates = new Vector<Integer>();
	for (int x = 0; x < inputShape.getWidth(); x++) {
	    for (int y = 0; y < inputShape.getHeight(); y++) {
		if (inputShape.delta(x, y) > 0 && x>10 && y>10 && width-x>10 && height-y>10) {
		    xCoordinates.add(x);
		    yCoordinates.add(y);
		}
	    }
	}
	IJ.log("pts " + xCoordinates.size());
	IJ.log("ok");
	fit();
    }

    
    public SuperEllipseFitter(Vector<Integer> xCoordinates, Vector<Integer> yCoordinates, int width, int height){
	this.xCoordinates=xCoordinates;
	this.yCoordinates=yCoordinates;
	this.width=width;
	this.height=height;
    }
    /**
     * Populates the coordinate array using a mask
     * @param mask
     */
    public SuperEllipseFitter(boolean[][] mask) {
	this.inputShape = new ImplicitShape2D(mask);
	this.setInputRoi(inputShape.getRoi(true));
	this.width = mask.length;
	this.height = mask[0].length;
	xCoordinates = new Vector<Integer>();
	yCoordinates = new Vector<Integer>();
	for (int x = 0; x < inputShape.getWidth(); x++) {
	    for (int y = 0; y < inputShape.getHeight(); y++) {
		if (inputShape.delta(x, y) > 0) {
		    xCoordinates.add(x);
		    yCoordinates.add(y);
		}
	    }
	}
	fit();
    }

    /**
     * Set weight based on curvature
     * @param x
     * @param y
     * @return
     */
    public double weight(int x, int y) {
	if (Math.abs(x) < 8 || Math.abs(y) < 8 || Math.abs(x - width) < 8
		|| Math.abs(y - height) < 8)
	    return 0;
	return Math.exp(-this.inputShape.getCurvature(x, y));
    }

    /**
     * Get the gradient of the energy with respect to the current ellipse
     * parameters
     * 
     */
    public double[] GradF(int x, int y, double[] params) {
	double X = X(x, y, params);
	double Y = Y(x, y, params);
	double xc = params[0];
	double yc = params[1];
	double a = params[2];
	double b = params[3];
	double theta = params[4];
	double eps = params[5];
	double factor = (1 - w1 - w2) * Q0(x, y, params);
	double factorx = 1
		/ eps
		* Math.pow(
			((x - xc) * Math.cos(theta) - (y - yc)
				* Math.sin(theta))
				/ a, (1 - eps) / eps);
	double factory = 1
		/ eps
		* Math.pow(
			((y - yc) * Math.cos(theta) + (x - xc)
				* Math.sin(theta))
				/ b, (1 - eps) / eps);
	// xc yc a b theta eps

	double dxc = factor
		* (X * factorx * (-Math.cos(theta)) / a - Y * factory
			* Math.sin(theta) / b);
	double dyc = factor
		* (X * factorx * (Math.sin(theta) / a) - Y * factory
			* Math.cos(theta) / b);
	double da = factor * X * -X / a / eps;
	double db = factor * Y * -Y / b / eps;
	double dtheta = factor
		* (-X
			* factorx
			* ((x - xc) * Math.sin(theta) + (y - yc)
				* Math.cos(theta)) / a + Y
			* factory
			/ b
			* (-(y - yc) * Math.sin(theta + (x - xc)
				* Math.cos(theta))));
	double deps = factor
		* (-X * X / eps * Math.log(X + Double.MIN_VALUE) - Y * Y / eps
			* Math.log(Y + Double.MIN_VALUE));
	return new double[] { dxc, dyc, da, db, dtheta, deps };
    }

    /**
     * Total gradient
     * 
     * @param params
     * @return
     */
    public double[] GradF(double[] params) {
	double[] grad = new double[] { 0, 0, 0, 0, 0, 0 };
	int N = xCoordinates.size();
	double sumweights = 0;

	double weight;
	for (int i = 0; i < N; i++) {
	    weight = weight(xCoordinates.get(i), yCoordinates.get(i));
	    double[] ptgrad = GradF(xCoordinates.get(i), yCoordinates.get(i),
		    params);
	    sumweights += weight;
	    grad[0] += weight * ptgrad[0];
	    grad[1] += weight * ptgrad[1];
	    grad[2] += weight * ptgrad[2];
	    grad[3] += weight * ptgrad[3];
	    grad[4] += weight * ptgrad[4];
	    grad[5] += weight * ptgrad[5];
	}

	grad[0] = grad[0] / sumweights;
	grad[1] = grad[1] / sumweights;
	grad[2] = grad[2] / sumweights;
	grad[3] = grad[3] / sumweights;
	grad[4] = grad[4] / sumweights;
	grad[5] = grad[5] / sumweights;

	return grad;
    }

    // Just get the entire Hessian matrix
    public double[][] HesF(double[] params) {
	int N = xCoordinates.size();
	double xc = params[0];
	double yc = params[1];
	double a = params[2];
	double b = params[3];
	double theta = params[4];
	double eps = params[5];

	double xcxc = 0;
	double xcyc = 0;
	double xca = 0;
	double xcb = 0;
	double xctheta = 0;
	double xceps = 0;
	double ycyc = 0;
	double yca = 0;
	double ycb = 0;
	double yctheta = 0;
	double yceps = 0;
	double aa = 0;
	double ab = 0;
	double atheta = 0;
	double aeps = 0;
	double bb = 0;
	double btheta = 0;
	double beps = 0;
	double thetatheta = 0;
	double thetaeps = 0;
	double epseps = 0;

	double sumweights = 0;

	double weight;

	for (int i = 0; i < N; i++) {
	    weight = weight(xCoordinates.get(i), yCoordinates.get(i));
	    sumweights += weight;
	    int x = xCoordinates.get(i);
	    int y = yCoordinates.get(i);
	    double X = X(x, y, params);
	    double Y = Y(x, y, params);

	    double factorx = 1
		    / eps
		    * Math.pow(
			    ((x - xc) * Math.cos(theta) - (y - yc)
				    * Math.sin(theta))
				    / a, (1 - eps) / eps);
	    double factory = 1
		    / eps
		    * Math.pow(
			    ((y - yc) * Math.cos(theta) + (x - xc)
				    * Math.sin(theta))
				    / b, (1 - eps) / eps);
	    double dxc = X * factorx * -Math.cos(theta) / a + Y * factory
		    * Math.sin(theta) / b;
	    double dyc = X * factorx * Math.sin(theta) / a + Y * factory
		    * -Math.cos(theta) / b;
	    double da = X * -X / a / eps;
	    double db = Y * -Y / b / eps;
	    double dtheta = (-X * factorx
		    * ((x - xc) * Math.sin(theta) + (y - yc) * Math.cos(theta))
		    / a + Y
		    * factory
		    / b
		    * (-(y - yc) * Math.sin(theta + (x - xc) * Math.cos(theta))));
	    double deps = (-X * X / eps * Math.log(X + Double.MIN_VALUE) - Y
		    * Y / eps * Math.log(Y + Double.MIN_VALUE));

	    xcxc += weight * (1 - w1 - w2) * dxc * dxc;
	    xcyc += weight * (1 - w1 - w2) * dxc * dyc;
	    xca += weight * (1 - w1 - w2) * dxc * da;
	    xcb += weight * (1 - w1 - w2) * dxc * db;
	    xctheta += weight * (1 - w1 - w2) * dxc * dtheta;
	    xceps += weight * (1 - w1 - w2) * dxc * deps;

	    ycyc += weight * (1 - w1 - w2) * dyc * dyc;
	    yca += weight * (1 - w1 - w2) * dyc * da;
	    ycb += weight * (1 - w1 - w2) * dyc * db;
	    yctheta += weight * (1 - w1 - w2) * dyc * dtheta;
	    yceps += weight * (1 - w1 - w2) * dyc * deps;

	    aa += weight * (1 - w1 - w2) * da * da;
	    ab += weight * (1 - w1 - w2) * da * db;
	    atheta += weight * (1 - w1 - w2) * da * dtheta;
	    aeps += weight * (1 - w1 - w2) * da * deps;

	    bb += weight * (1 - w1 - w2) * db * db;
	    btheta += weight * (1 - w1 - w2) * db * dtheta;
	    beps += weight * (1 - w1 - w2) * db * deps;

	    thetatheta += weight * (1 - w1 - w2) * dtheta * dtheta;
	    thetaeps += weight * (1 - w1 - w2) * dtheta * deps;

	    epseps += weight * (1 - w1 - w2) * deps * deps;

	}
	return new double[][] {
		new double[] { xcxc / sumweights, xcyc / sumweights,
			xca / sumweights, xcb / sumweights,
			xctheta / sumweights, xceps / sumweights },
		new double[] { xcyc / sumweights, ycyc / sumweights,
			yca / sumweights, ycb / sumweights,
			yctheta / sumweights, yceps / sumweights },
		new double[] { xca / sumweights, yca / sumweights,
			aa / sumweights, ab / sumweights, atheta / sumweights,
			aeps / sumweights },
		new double[] { xcb / sumweights, ycb / sumweights,
			ab / sumweights, bb / sumweights, btheta / sumweights,
			beps / sumweights },
		new double[] { xctheta / sumweights, yctheta / sumweights,
			atheta / sumweights, btheta / sumweights,
			thetatheta / sumweights, thetaeps / sumweights },
		new double[] { xceps / sumweights, yceps / sumweights,
			aeps / sumweights, beps / sumweights,
			thetaeps / sumweights, epseps / sumweights } };
    }

    /**
     * F(x,y)
     * 
     * @param x
     * @param y
     * @param params
     * @return
     */
    public double F(int x, int y, double[] params) {
	return (1 - w1 - w2) * Math.pow(Q0(x, y, params), 2) + w1
		* Math.pow(Q1(x, y, params), 2) + w2
		* Math.pow(Q2(x, y, params), 2);
    }

    /**
     * Use the Ellipse class to fit an Ellipse. This can be used to get intial parameters for our fitting
     * @return [xc yc a b theta]
     */
    public double[] ellipseFitter() {
	if(xCoordinates.size()==0 || yCoordinates.size()==0)
	    return new double[]{0,0,0,0,0};
	Ellipse el = new Ellipse(xCoordinates, yCoordinates);
	xc = inputShape.inertialCenter[0];
	yc = inputShape.inertialCenter[1];
	return new double[] { el.getCenterX(), el.getCenterY(), el.getMajor(), el.getMinor(),
		-el.getAngle() };
    }

    /**
     * 
     * @param x
     * @param y
     * @param params
     *            xc yc a b theta eps
     * @return
     */
    public double Q0(int x, int y, double[] params) {
	double xc = params[0];
	double yc = params[1];
	double a = params[2];
	double b = params[3];
	double theta = params[4];
	double eps = params[5];
	return Math.pow(
		((x - xc) * Math.cos(theta) - (y - yc) * Math.sin(theta)) / a,
		2 / eps)
		+ Math.pow(
			((x - xc) * Math.sin(theta) + (y - yc)
				* Math.cos(theta))
				/ b, 2 / eps) - 1;

    }

    /**
     * Q1 evaluated at a single point x y
     * 
     * @param x
     * @param y
     * @param params
     * @return
     */
    public double Q1(int x, int y, double[] params) {
	double xc = params[0];
	double yc = params[1];
	double a = params[2];
	double b = params[3];
	double theta = params[4];
	double eps = params[5];
	return 0;
    }

    public double Q2(int x, int y, double[] params) {
	double xc = params[0];
	double yc = params[1];
	double a = params[2];
	double b = params[3];
	double theta = params[4];
	double eps = params[5];
	return 0;
    }

    public double X(int x, int y, double[] params) {
	double xc = params[0];
	double yc = params[1];
	double a = params[2];
	double theta = params[4];
	double eps = params[5];

	return Math
		.pow((x - xc) * Math.cos(theta) / a - (y - yc)
			* Math.sin(theta) / a, 1 / eps);

    }

    public double Y(int x, int y, double[] params) {
	double xc = params[0];
	double yc = params[1];
	// double a = params[2];
	double b = params[3];
	double theta = params[4];
	double eps = params[5];

	return Math
		.pow((x - xc) * Math.sin(theta) / b + (y - yc)
			* Math.cos(theta) / b, 1 / eps);

    }

    /**
     * LMA algorithm step
     * 
     * @param params
     * @return
     */
    public double[][] getLMAstep(double[] params) {
	double[][] hes = this.HesF(params);
	double[] grad = this.GradF(params);

	Matrix H = new Matrix(hes);
	Matrix b = new Matrix(grad, 6);

	Matrix d = H.solve(b);

	return d.getArray();
    }

    public boolean[][] getMask() {
	boolean[][] mask = new boolean[width][height];
	for (int x = 0; x < mask.length; x++) {
	    for (int y = 0; y < mask[0].length; y++) {
		mask[x][y] = Math.pow((x - xc) * Math.cos(theta) / a - (y - yc)
			* Math.sin(theta) / a, 2 / eps)
			+ Math.pow((y - yc) * Math.cos(theta) / b + (x - xc)
				* Math.sin(theta) / b, 2 / eps) <= 1 ? true
			: false;
	    }
	}
	return mask;
    }

    /**
     * Get level set embedding of the fitted superellipse
     * 
     * @return
     */
    public ImplicitShape2D getLevelSets() {
	return new ImplicitShape2D(getMask());
    }

    /** Draws the ellipse on the specified image. */
    public void drawEllipse(ImageProcessor ip) {
	ImplicitShape2D shape = new ImplicitShape2D(getMask());
	Roi r = shape.getRoi(false);
	ip.draw(r);

    }

    public Roi getRoi() {

	return this.fittedRoi;
    }

    public Roi getInputRoi() {
	return inputRoi;
    }

    public void setInputRoi(Roi inputRoi) {
	this.inputRoi = inputRoi;
    }

    public ImageProcessor getIp() {
	return ip;
    }

    public void setIp(ImageProcessor ip) {
	this.ip = ip;
    }

    public double[] getParams() {
	if(Double.isNaN(xc)||Double.isNaN(yc)||Double.isNaN(a)||Double.isNaN(b))
	    return new double[] {0,0,0,0,0,1};
		
	return new double[] { xc, yc, a, b, theta, eps };
    }

}