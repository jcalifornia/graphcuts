package ucla.chou.graphcutshapes;

import java.awt.Polygon;

import ij.IJ;
import ij.ImagePlus;
import ij.gui.PolygonRoi;
import ij.gui.Roi;
import ij.gui.Wand;
import ij.process.ByteProcessor;
import ij.process.FloatProcessor;
import ij.process.ImageProcessor;

/*
 * This class describes a shape object stored implicitly as
 * a signed Euclidean distance function
 */

public class ImplicitShape2D {
    public float[][] signedDistance;

    double inertialScale;
    double[] inertialOrientation;
    double inertialAngle; // inertialAngle of the inertialOrientation
    double[] inertialCenter;
    double lambda = 2;
    double alighmentlambda = 2; // use this for alignment
    final int d = 2; // Dimension... need to rewrite a few things for 3D
    public double mass;
    int width;
    int height;

    double epsilon = 0.5; // for delta and heaviside function approx

    public void setAlignmentlambda(double val) {
	this.alighmentlambda = val;
    }

    public void setlambda(double val) {
	this.lambda = val;
    }

    public ImplicitShape2D(Roi r, int width, int height) {
	boolean[][] mask = new boolean[width][height];
	for (int x = 0; x < width; x++) {
	    for (int y = 0; y < height; y++) {
		mask[x][y] = r.contains(x, y);
	    }
	}
	this.signedDistance = masktoSignedDistance(mask);
	calculateInertialNormalization();
    }

    public ImplicitShape2D(boolean[][] mask) {
	this.signedDistance = masktoSignedDistance(mask);
	this.width = mask.length;
	this.height = mask[0].length;

	calculateInertialNormalization();

    }

    /*
     * Calculate inertialCenter of mass
     */
    double[] centerofMass() {
	int w = signedDistance.length;
	int h = signedDistance[0].length;
	double xc = 0;
	double yc = 0;
	double normal = 0;
	for (int i = 0; i < w; i++) {
	    for (int j = 0; j < h; j++) {
		if (in(i, j)) {
		    xc += rho(i, j) * i;
		    yc += rho(i, j) * j;
		    normal += rho(i, j);

		}
	    }
	}
	double[] center = { xc / normal, yc / normal };
	this.mass = normal;
	return center;
    }

    /*
     * Calculate inertialOrientation
     */
    public double[] inertialOrientation() {
	double[] orientation = new double[2];
	// Compute the second moment about the inertialCenter
	int w = signedDistance.length;
	int h = signedDistance[0].length;
	double xc = this.inertialCenter[0];
	double yc = this.inertialCenter[1];
	double q1 = 0;
	double q2 = 0;
	// double q3=0;
	double q4 = 0;

	double[] normGradPhi = new double[] { 0, 0 };
	for (int x = 0; x < w; x++) {
	    for (int y = 0; y < h; y++) {
		if (in(x, y)) {
		    // get local curvature

		    double[] gradphi = this.getNormGradient(x, y);
		    normGradPhi[0] += gradphi[0] * rho(x, y);
		    normGradPhi[1] += gradphi[1] * rho(x, y);

		    // Get matrix entries for 2nd moment
		    q1 += (x - xc) * (x - xc) * rho(x, y);
		    q4 += (y - yc) * (y - yc) * rho(x, y);
		    q2 += (x - xc) * (y - yc) * rho(x, y);
		}
	    }
	}

	// double theta = 0.5*Math.atan2(q2,(q1-q4));

	double Tr = q1 + q4;
	double det = q1 * q4 - q2 * q2;

	double l1 = 0.5 * (Tr + Math.sqrt(Tr * Tr - 4 * det));
	double l2 = 0.5 * (Tr - Math.sqrt(Tr * Tr - 4 * det));

	// find the two eigenvectors
	double x21 = q1 - l1;
	double norm1 = Math.sqrt(q2 * q2 + (l1 - q1) * (l1 - q1));

	double x22 = q1 - l2;
	double norm2 = Math.sqrt(q2 * q2 + (l2 - q1) * (l2 - q1));
	// which one is larger?

	if (Math.abs(l1) >= Math.abs(l2)) {

	    orientation[0] = q2 / norm1;
	    orientation[1] = -x21 / norm1;
	} else {
	    orientation[0] = q2 / norm2;
	    orientation[1] = -x22 / norm2;
	}

	double massup = 0;
	double massdown = 0;

	for (int x = 0; x < width; x++) {
	    for (int y = 0; y < height; y++) {
		// project point onto the orientation line, see which side it is
		// above or below the center
		if (orientation[0] * (x - this.inertialCenter[0])
			+ orientation[1] * (y - this.inertialCenter[1]) < 0) {
		    massdown += Math.pow(Math.abs(this.get(x, y)), 0);
		} else {
		    massup += Math.pow(Math.abs(this.get(x, y)), 0);
		}
	    }
	}

	if (massdown < massup) {
	    orientation[1] *= -1;
	    orientation[0] *= -1;
	}

	// compare it to normGradPhi to determine which direction it should go
	// in!

	double theta = Math.atan2(orientation[1], orientation[0]);
	theta += theta < 0 ? Math.PI * 2 : 0;
	double theta2 = Math.atan2(normGradPhi[1], normGradPhi[0]);
	theta2 += theta2 < 0 ? Math.PI * 2 : 0;

	// reverse it if inertialOrientation is not in same half-circle as
	// normGradPhi

	if (theta2 < Math.PI / 2) {
	    if (theta > theta2 + Math.PI / 2 && theta < theta + 3 * Math.PI / 2) {
		orientation[1] *= -1;
		orientation[0] *= -1;
	    }
	} else if (theta2 < 3 * Math.PI / 2) {
	    if (theta < theta2 - Math.PI / 2 || theta > theta2 + Math.PI / 2) {
		orientation[1] *= -1;
		orientation[0] *= -1;
	    }

	} else {

	    if (theta < theta2 - Math.PI / 2 && theta > Math.PI / 2 - theta2) {
		orientation[0] *= -1;
		orientation[1] *= -1;
	    }
	}
	// this.inertialAngle=theta;
	this.inertialAngle = Math.atan2(orientation[1], orientation[0]);
	return orientation;
    }

    /*
     * Calculate inertialScale
     */
    double inertialScale() {
	double sum = 0;
	for (int x = 0; x < this.signedDistance.length; x++) {
	    for (int y = 0; y < this.signedDistance[0].length; y++) {
		sum += rho(x, y);
		;
	    }
	}
	return Math.pow(sum, 1.0 / (d + lambda));

    }

    public double rho(int x, int y) {

	return this.lambda == 0 ? (in(x, y) ? 1 : 0) : (in(x, y) ? Math.pow(
		this.get(x, y), this.lambda) : 0);

    }

    boolean[][] getMask() {
	int l = this.signedDistance.length;
	int w = this.signedDistance[0].length;
	boolean[][] mask = new boolean[l][w];
	for (int x = 0; x < l; x++) {
	    for (int y = 0; y < w; y++) {
		mask[x][y] = in(x, y);
	    }
	}
	return mask;

    }

    /**
     * Use this is the two templates have already been aligned
     * 
     * @param Lambda
     * @return
     */
    double untransformedDistance(ImplicitShape2D Lambda) {
	double d = 0;

	double mass = 0;
	for (int x = 0; x < width; x++) {
	    for (int y = 0; y < height; y++) {

		d += (Math.abs(heaviside(x, y) - Lambda.heaviside(x, y))
			* this.rho(x, y) + Lambda.delta(x, y))
			* this.rho(x, y);
		mass += this.rho(x, y);

	    }

	}

	return d / Math.pow(this.inertialScale,this.lambda);

    }

    /**
     * Distance between this shape and other shape computes first the transform
     * 
     * @param Lambda
     * @param alpha
     * @param c
     * @param omega
     * @return
     */
    double distance(ImplicitShape2D Lambda) {
	// First step is to get the transform coords for the other shape
	double[] alignment = alignByNewton(Lambda);
	double[][] R = new double[][] {
		new double[] { Math.cos(alignment[3]), -Math.sin(alignment[3]) },
		new double[] { Math.sin(alignment[3]), Math.cos(alignment[3]) } };
	// now calculate the distance
	double[] sprime;
	double d = 0;

	double mass = 0;
	for (int x = 0; x < width; x++) {
	    for (int y = 0; y < height; y++) {
		sprime = new double[] {
			alignment[0]
				* (R[0][0] * (x - alignment[1]) + R[0][1]
					* (y - alignment[2])),
			alignment[0]
				* (R[1][0] * (x - alignment[1]) + R[1][1]
					* (y - alignment[2])) };

		d += (Math.abs(heaviside(x, y)
			- Lambda.heaviside((int) Math.round(sprime[0]),
				(int) Math.round(sprime[1]))) + delta(x, y))
			* Lambda.rho((int) sprime[0], (int) sprime[1]);
		mass += Lambda.rho((int) sprime[0], (int) sprime[1]);

	    }

	}

	return d / Math.pow(this.inertialScale,this.lambda);

    }

    void calculateInertialNormalization() {
	this.inertialCenter = centerofMass();
	this.inertialScale = inertialScale();
	this.inertialOrientation = inertialOrientation();

    }

    public float[][] masktoSignedDistance(boolean[][] mask) {
	int width = mask.length;
	int height = mask[0].length;
	float[][] grid = new float[width][height];
	for (int x = 0; x < width; x++) {
	    for (int y = 0; y < height; y++) {
		if (mask[x][y])
		    grid[x][y] = Float.MAX_VALUE;
		else
		    grid[x][y] = -Float.MAX_VALUE;
	    }

	}
	return (ExactDistanceTransform(grid));
    }

    /**
     * Exact distance transformation in 1D taken from the Fiji project
     * 
     * @param f
     * @return
     */
    public static float[] EDTTransform1D(float[] f) {
	float d[] = new float[f.length];
	float[] fNeg = new float[f.length + 1];
	float[] zNeg = new float[f.length + 1];
	int[] yNeg = new int[f.length + 1];
	float[] fPos = new float[f.length + 1];
	float[] zPos = new float[f.length + 1];
	int[] yPos = new int[f.length + 1];

	fNeg[0] = -Float.MAX_VALUE;
	yNeg[0] = -1;
	zNeg[0] = Float.MAX_VALUE;
	int kNeg = 0, kPos = 0;
	fPos[0] = Float.MAX_VALUE;
	yPos[0] = -1;
	zPos[0] = Float.MAX_VALUE;

	float fval, fvalNeg, fvalPos, s;

	for (int q = 0; q < f.length; q++) {
	    fval = f[q];
	    fvalNeg = fval < 0 ? fval : 0f;
	    s = ((fvalNeg - q * q) - (fNeg[kNeg] - yNeg[kNeg] * yNeg[kNeg]))
		    / -2 / (q - yNeg[kNeg]);
	    for (;;) {
		// calculate the intersection
		s = ((fvalNeg - q * q) - (fNeg[kNeg] - yNeg[kNeg] * yNeg[kNeg]))
			/ -2 / (q - yNeg[kNeg]);
		if (s > zNeg[kNeg])
		    break;
		if (--kNeg < 0)
		    break;
	    }

	    kNeg++;
	    yNeg[kNeg] = q;
	    fNeg[kNeg] = fvalNeg;
	    zNeg[kNeg] = s;
	    fvalPos = fval > 0 ? fval : 0f;
	    for (;;) {
		// calculate the intersection
		s = ((fvalPos + q * q) - (fPos[kPos] + yPos[kPos] * yPos[kPos]))
			/ 2 / (q - yPos[kPos]);
		if (s > zPos[kPos])
		    break;
		if (--kPos < 0)
		    break;
	    }
	    kPos++;
	    yPos[kPos] = q;
	    fPos[kPos] = fvalPos;
	    zPos[kPos] = s;
	}
	zNeg[++kNeg] = Float.MAX_VALUE;
	zPos[++kPos] = Float.MAX_VALUE;

	int iNeg = 0, iPos = 0;
	for (int q = 0; q < f.length; q++) {
	    while (zNeg[iNeg + 1] < q)
		iNeg++;
	    while (zPos[iPos + 1] < q)
		iPos++;

	    d[q] = f[q] < 0 ? -(q - yNeg[iNeg]) * (q - yNeg[iNeg]) + fNeg[iNeg]
		    : (q - yPos[iPos]) * (q - yPos[iPos]) + fPos[iPos];
	    // d[q] = d[q]<0? 0.5f - (float)Math.sqrt(-d[q]): -0.5f + (float)
	    // Math.sqrt(d[q]);
	}

	return d;

    }

    /**
     * Perform the exact distance transformation given a signed distance
     * 
     * @param signedDistance
     * @return
     */
    public static float[][] ExactDistanceTransform(float[][] grid) {
	// Restore the signed distance

	float[][] newgrid = new float[grid.length][grid[0].length];

	float[] c = new float[grid.length];
	float[] r = new float[grid[0].length];
	for (int x = 0; x < grid.length; x++) {
	    for (int y = 0; y < grid[0].length; y++) {
		r[y] = grid[x][y] < 0 ? -Float.MAX_VALUE : Float.MAX_VALUE;
	    }
	    float[] d1 = EDTTransform1D(r);
	    for (int y = 0; y < grid[0].length; y++) {
		newgrid[x][y] = d1[y];
	    }
	}

	for (int y = 0; y < grid[0].length; y++) {
	    for (int x = 0; x < grid.length; x++) {
		c[x] = newgrid[x][y];
	    }
	    float[] d2 = EDTTransform1D(c);
	    for (int x = 0; x < grid.length; x++) {
		newgrid[x][y] = d2[x];
	    }
	}

	for (int x = 0; x < grid.length; x++) {
	    for (int y = 0; y < grid[0].length; y++) {
		newgrid[x][y] = (float) (newgrid[x][y] < 0 ? 0.5f - Math
			.sqrt(-newgrid[x][y]) : -0.5f
			+ Math.sqrt(newgrid[x][y]));
	    }
	}

	return newgrid;

    }

    public double getFast(double x, double y) {
	return get((int) Math.round(x), (int) Math.round(y));
    }

    /**
     * Get the value here using bilinear interpolation!
     * 
     * @param x
     * @param y
     * @return
     */
    public double get(double x, double y) {
	int x1 = (int) Math.floor(x);
	int x2 = x1 + 1;
	int y1 = (int) Math.floor(y);
	int y2 = y1 + 2;

	double fr1 = (x2 - x) / (x2 - x1) * get(x1, y1) + (x - x1) / (x2 - x1)
		* get(x2, y1);
	double fr2 = (x2 - x) / (x2 - x1) * get(x1, y2) + (x - x1) / (x2 - x2)
		* get(x2, y2);

	return (y2 - y) / (y2 - y1) * fr1 + (y - y1) / (y2 - y1) * fr2;

    }

    /**
     * Get the value of the signed distance at the specified coordinates
     * 
     * @param x
     * @param y
     * @return
     */
    public double get(int x, int y) {
	// Should really do some interpolation here
	return this.signedDistance[Math.max(0,
		Math.min(x, signedDistance.length - 1))][Math.max(0,
		Math.min(y, signedDistance[0].length - 1))];
    }

    int getWidth() {
	return width;
    }

    int getHeight() {
	return height;
    }

    boolean in(int x, int y) {
	if (x < 0 || y < 0 || x >= signedDistance.length
		|| y >= signedDistance[0].length)
	    return false;
	return get(x, y) >= 0 ? true : false;
    }

    public ByteProcessor getMaskProcessor() {
	ByteProcessor ip2 = new ByteProcessor(width, height);
	for (int x = 0; x < width; x++) {
	    for (int y = 0; y < height; y++) {
		ip2.set(x, y, signedDistance[x][y] >= 0 ? 0 : 255);
	    }
	}
	return (ip2);
    }

    /**
     * 
     * @return FloatProcessor with pixels set to the signed distance
     */
    public FloatProcessor getLSFloatProcessor() {
	FloatProcessor ip2 = new FloatProcessor(width, height);
	for (int x = 0; x < width; x++) {
	    for (int y = 0; y < height; y++) {
		ip2.setf(x, y, signedDistance[x][y]);
	    }
	}
	ip2.resetMinAndMax();
	ip2.crop();
	return (ip2);
    }

    public ImagePlus getLSImagePlus() {
	FloatProcessor ip2 = getLSFloatProcessor();
	ImagePlus imp = new ImagePlus("Signed distance", ip2);
	return imp;

    }

    /**
     * Get the ROI of the current level set
     * 
     * @param reset
     *            Reset the signed distance function after obtaining the ROI
     * @return
     */
    public PolygonRoi getRoi(boolean reset) {
	// Returns the ROI corresponding to the zero level set
	ImageProcessor ip2 = this.getLSFloatProcessor();
	int maxpt[] = { 0, 0 };
	/*
	 * Find location of maximum point
	 */
	for (int x = 0; x < width; x++) {
	    for (int y = 0; y < height; y++) {
		if (signedDistance[x][y] > signedDistance[maxpt[0]][maxpt[1]]) {
		    maxpt[0] = x;
		    maxpt[1] = y;
		}
	    }
	}
	Wand w = new Wand(ip2);
	// IJ.log(maxpt[0] + "" + signedDistance[maxpt[0]][maxpt[1]]);
	w.autoOutline(maxpt[0], maxpt[1], 0,
		signedDistance[maxpt[0]][maxpt[1]], Wand.FOUR_CONNECTED);
	if (w.npoints == 0) {
	    return null;
	}

	Polygon p = new Polygon(w.xpoints, w.ypoints, w.npoints);
	PolygonRoi poly = new PolygonRoi(p, Roi.POLYGON);
	if (reset) {
	    for (int x = 0; x < width; x++) {
		for (int y = 0; y < height; y++) {
		    if (poly.contains(x, y))
			signedDistance[x][y] = Float.MAX_VALUE;
		    else
			signedDistance[x][y] = -Float.MAX_VALUE;
		}

	    }

	    this.signedDistance = ExactDistanceTransform(this.signedDistance);
	    this.calculateInertialNormalization();
	}

	return poly;
    }

    /**
     * Translate current shape to the characteristics of a new shape also
     * reproduce its width and height TODO
     * 
     * @param shape
     */
    public void translateInertially(ImplicitShape2D shape) {
	// Obtain new coordinate by rotation
	boolean[][] mask = new boolean[shape.width][shape.height];
	/*
	 * need to multiply this.shape by this factor in order to match shape
	 */
	double scalingfactor = Math.pow(shape.inertialScale
		/ this.inertialScale, -1);// ,-(this.d+this.lambda));
	// double scalingfactor = 1;
	// how much we need to rotate this shape to match the reference
	double rotation = shape.inertialAngle + this.inertialAngle;
	double newangle;
	double newdistance;
	/*
	 * In new shape coordinates
	 */
	for (int x = 0; x < shape.width; x++) {
	    for (int y = 0; y < shape.height; y++) {
		newangle = Math.atan2(y - shape.inertialCenter[1], x
			- shape.inertialCenter[0])
			+ rotation;
		newdistance = Math.sqrt(Math
			.pow(x - shape.inertialCenter[0], 2)
			+ Math.pow(y - shape.inertialCenter[1], 2))
			* scalingfactor;
		// mask[x][y] = this.in(x-shape.center[0]+this.center[0],
		// y+this.center[1]-shape.center[1]);
		mask[x][y] = this.in(this.inertialCenter[0] + newdistance
			* Math.cos(newangle), this.inertialCenter[1]
			+ newdistance * Math.sin(newangle));
	    }
	}
	this.signedDistance = this.masktoSignedDistance(mask);
	this.width = shape.width;
	this.height = shape.height;
	this.calculateInertialNormalization();

    }

    /**
     * Should do something better here but whatever
     * 
     * @param e
     * @param f
     * @return
     */
    public boolean in(double e, double f) {
	return this.getFast(e, f) >= 0;
    }

    /**
     * Get central curvature
     * 
     * @param x
     * @param y
     * @return
     */
    public double getCurvature(int x, int y) {
	double dphiX = (get(x + 1, y) - get(x - 1, y)) / 2;
	double dphiY = (get(x, y + 1) - get(x, y - 1)) / 2;
	double dphiXX = ((get(x + 1, y) + get(x - 1, y) - (2 * get(x, y)))) / 4;
	double dphiYY = ((get(x, y + 1) + get(x, y - 1) - (2 * get(x, y)))) / 4;
	double dphiXY = (get(x + 1, y + 1) - get(x + 1, y - 1)
		- get(x - 1, y + 1) + get(x - 1, y - 1)) / 4;

	double curvature = (dphiXX * dphiY * dphiY + dphiYY * dphiX * dphiX - 2
		* dphiX * dphiY * dphiXY)
		/ (Math.pow(dphiX * dphiX + dphiY * dphiY, 1.5) + 1e-10);

	return curvature; // * deltaPhi;
    }

    /**
     * Central differences gradient at point (x,y)
     * 
     * @param x
     * @param y
     * @return
     */
    public double[] getGradient(int x, int y) {
	double dphiX = (get(x + 1, y) - get(x - 1, y)) / 2;
	double dphiY = (get(x, y + 1) - get(x, y - 1)) / 2;
	return new double[] { dphiX, dphiY };
    }

    /**
     * Normalized central differences gradient at point (x,y)
     * 
     * @param x
     * @param y
     * @return
     */
    public double[] getNormGradient(int x, int y) {
	double dphiX = (get(x + 1, y) - get(x - 1, y)) / 2;
	double dphiY = (get(x, y + 1) - get(x, y - 1)) / 2;
	double norm = Math.sqrt(dphiY * dphiY + dphiX * dphiX);
	if (norm == 0)
	    return new double[] { 0, 0 };
	return new double[] { dphiX / norm, dphiY / norm };
    }

    public double length() {
	double tot = 0;
	for (int x = 0; x < width; x++)
	    for (int y = 0; y < height; y++)
		tot += delta(x, y);
	return tot;
    }

    public double delta(int x, int y) {
	if (Math.abs(get(x, y)) >= 2)
	    return 0;
	return (1 / Math.PI * this.epsilon / (this.epsilon * this.epsilon + get(
		x, y) * get(x, y)));
    }

    public double heaviside(int x, int y) {
	if (get(x, y) < -2)
	    return 0;
	else if (get(x, y) > 2)
	    return 1;
	else
	    return (0.5 + Math.atan(get(x, y) / this.epsilon) / Math.PI);
    }

    /**
     * Align using Newton-Raphson to minimize the energy
     * 
     * @param Lambda
     *            The shape to align to this shape. It finds the parameters for
     *            the transformation alpha R(s-c)
     */
    public double[] alignByNewton(ImplicitShape2D Lambda) {
	double alpha = Lambda.inertialScale / this.inertialScale;//
	// how much we need to rotate this shape to match the reference
	// double omega = Lambda.inertialAngle - this.inertialAngle;
	double omega = 0;

	if (Math.abs(omega) > Math.PI) {
	    if (omega > Math.PI)
		omega -= 2 * Math.PI;
	    else
		omega += 2 * Math.PI;
	}
	double[][] R = new double[][] {
		new double[] { Math.cos(omega), -Math.sin(omega) },
		new double[] { Math.sin(omega), Math.cos(omega) } };
	double[] c = new double[] {
		-(R[0][0] * Lambda.inertialCenter[0] - R[0][1]
			* Lambda.inertialCenter[1])
			/ alpha + this.inertialCenter[0],
		-(-R[1][0] * Lambda.inertialCenter[0] + R[1][1]
			* Lambda.inertialCenter[1])
			/ alpha + this.inertialCenter[1] };
	return (alignByNewton(Lambda, alpha, c, omega));

    }

    public double[] alignByNewton(ImplicitShape2D Lambda, double alpha0,
	    double[] c0, double omega0) {

	// Maybe use L\infty
	double gap = Double.MAX_VALUE;
	// Initialize using the inertial values
	double[] grad;
	double[][] invhes;
	double[] change;
	int iter = 0;

	double alpha = alpha0;
	double omega = omega0;
	double[] c = new double[] { c0[0], c0[1] };
	while (iter++ < 15) {
	    grad = getTransformGradient(Lambda, alpha, c, omega);
	    invhes = getInverseTransformHessian(Lambda, alpha, c, omega);
	    // IJ.log("Parameters: alpha " + alpha + "  sx  " + c[0] + "  sy  "
	    // + c[1] + " omega " + omega);

	    // Hgrad -- the problem is that this becomes singular near the
	    // minimum!
	    change = new double[] {
		    invhes[0][0] * grad[0] + invhes[0][1] * grad[1]
			    + invhes[0][2] * grad[2] + invhes[0][3] * grad[3],
		    invhes[1][0] * grad[0] + invhes[1][1] * grad[1]
			    + invhes[1][2] * grad[2] + invhes[1][3] * grad[3],
		    invhes[2][0] * grad[0] + invhes[2][1] * grad[1]
			    + invhes[2][2] * grad[2] + invhes[2][3] * grad[3],
		    invhes[3][0] * grad[0] + invhes[3][1] * grad[1]
			    + invhes[3][2] * grad[2] + invhes[3][3] * grad[3] };

	    //
	    gap = Math.sqrt(change[0] * change[0] * 100 * 100 + change[1]
		    * change[1] + change[2] * change[2] + change[3] * change[3]
		    * 100 * 100);

	    /**
	     * The Newton method jumps are too large, do some adhoc
	     * regularization should really use LMA here instead
	     */
	    if (gap > 3) {
		change[0] /= 10;
		change[1] /= 10;
		change[2] /= 10;
		change[3] /= 10;
	    }

	    alpha -= change[0];
	    c[0] -= change[1];
	    c[1] -= change[2];
	    omega -= change[3];

	    if (Math.abs(omega) > Math.PI) {
		if (omega > Math.PI)
		    omega -= 2 * Math.PI;
		else
		    omega += 2 * Math.PI;
	    }

	    /**
	     * meh, close enough - maybe use objective function instead here?
	     */
	    if (Math.abs(change[0] / alpha) < 0.001)
		break;
	    // if(Math.abs(c[0])>10) c[0]=0;
	    // else if(Math.abs(c[1])>10) c[1]=0;
	    // if (Math.abs(alpha) >1.2 || alpha<0.9) alpha = 1; // reset it;
	    if (Math.abs(alpha) < 0.1)
		alpha *= -1; // maybe need to flip it around

	}

	// check for sanity... if the alignment failed, ie, the center of the
	// shape is now out of the image,
	// just return params so that

	IJ.log("Completed shape aignment in " + iter + " iterations");
	// IJ.log("Parameters: alpha  " + alpha + "  sx  " + c[0] + "  sy  "
	// + c[1] + "  omega  " + omega);
	return new double[] { alpha, c[0], c[1], omega };

    }

    /**
     * Return the Hessian Matrix
     * 
     * @param x
     * @param y
     * @return
     */
    public double[][] HessianMatrix(int x, int y) {
	// wrt x
	double d2x2 = (-get(x + 2, y) + 16 * get(x + 1, y) - 30 * get(x, y)
		+ 16 * get(x - 1, y) - get(x - 2, y)) / 12;
	double d2y2 = (-get(x, y + 2) + 16 * get(x, y + 1) - 30 * get(x, y)
		+ 16 * get(x, y - 1) - get(x, y - 2)) / 12;

	double d2xy = 0.25 * (get(x + 1, y + 1) - get(x + 1, y - 1)
		- get(x - 1, y + 1) + get(x - 1, y - 1));

	return new double[][] { new double[] { d2x2, d2xy },
		new double[] { d2xy, d2y2 } }; // row major
    }

    /**
     * Implicit Affine transform the current shape - computes the new shape
     * under transformed coordinates of the current shape
     * 
     * @param inertialScale
     * @param inertialCenter
     * @param rotation
     */
    public void affineTransform(double alpha, double[] c, double omega) {
	// Obtain new coordinate by rotation
	boolean[][] mask = new boolean[width][height];
	/**
	 * Rotation matrix
	 */
	double[][] R = new double[][] {
		new double[] { Math.cos(omega), -Math.sin(omega) },
		new double[] { Math.sin(omega), Math.cos(omega) } };
	for (int x = 0; x < width; x++) {
	    for (int y = 0; y < height; y++) {
		int newx = (int) Math.round(alpha
			* (R[0][0] * (x - c[0]) + R[0][1] * (y - c[1])));
		int newy = (int) Math.round(alpha
			* (R[1][0] * (x - c[0]) + R[1][1] * (y - c[1])));

		mask[x][y] = this.in(newx, newy);

	    }
	}
	this.signedDistance = masktoSignedDistance(mask);
	this.calculateInertialNormalization();
    }

    /**
     * Inverse transform
     * 
     * @param alpha
     * @param c
     * @param omega
     */
    public void affineTransform2(double alpha, double[] c, double omega) {
	affineTransform2(alpha, c, omega, this.width, this.height);

    }

    public void affineTransform2(double alpha, double[] c, double omega,
	    int width, int height) {
	this.width=width;
	this.height=height;
	boolean[][] mask = new boolean[width][height];
	/**
	 * Rotation matrix
	 */
	double[][] R = new double[][] {
		new double[] { Math.cos(omega), Math.sin(omega) },
		new double[] { -Math.sin(omega), Math.cos(omega) } };
	for (int x = 0; x < width; x++) {
	    for (int y = 0; y < height; y++) {
		int newx = (int) Math.round(1 / alpha
			* (R[0][0] * (x) + R[0][1] * (y)) + c[0]);
		int newy = (int) Math.round(1 / alpha
			* (R[1][0] * (x) + R[1][1] * (y)) + c[1]);

		mask[x][y] = this.in(newx, newy);

	    }
	}
	this.signedDistance = masktoSignedDistance(mask);
	this.calculateInertialNormalization();

    }

    public void affineTransform(double alpha, double[] c, double omega,
	    int width, int height) {
	this.width = width;
	this.height = height;
	affineTransform(alpha, c, omega);
    }

    /**
     * Stuff for alignment This shape should be Omega, other shape is Lambda
     */

    /**
     * Return the gradient of shape energy wrt transformation params Don't
     * assume that the other shape is transformed.
     * 
     * @param Lambda
     * @return
     */
    public double[] getTransformGradient(ImplicitShape2D Lambda, double alpha,
	    double[] c, double omega) {
	double gradalpha = 0;
	double gradc1 = 0;
	double gradc2 = 0;
	double gradomega = 0;
	double lambda = this.alighmentlambda;

	/**
	 * Rotation matrix
	 */
	double[][] R = new double[][] {
		new double[] { Math.cos(omega), -Math.sin(omega) },
		new double[] { Math.sin(omega), Math.cos(omega) } };

	for (int x = 0; x < width; x++) {
	    for (int y = 0; y < height; y++) {
		// First, transform the coordinate to get coordinates in the
		// reference frame of omega
		// alpha R(s-c)

		int newx = (int) Math.round(alpha
			* (R[0][0] * (x - c[0]) + R[0][1] * (y - c[1])));
		int newy = (int) Math.round(alpha
			* (R[1][0] * (x - c[0]) + R[1][1] * (y - c[1])));

		double[] gradPhiLambda = Lambda.getGradient(newx, newy);

		double multfactor = lambda
			* (Math.pow(
				this.heaviside(x, y)
					- Lambda.heaviside(newx, newy), 2) + this
				.delta(x, y))
			* Math.pow(Math.abs(Lambda.get(newx, newy)), lambda - 1);
		multfactor *= Lambda.in(newx, newy) ? 1 : -1;

		gradalpha += multfactor
			* ((newx / alpha) * gradPhiLambda[0] + (newy / alpha)
				* gradPhiLambda[1]);

		/**
		 * WRT c -alpha*R
		 */

		gradc1 += -alpha
			* (R[0][0] * gradPhiLambda[0] + R[0][1]
				* gradPhiLambda[1]) * multfactor;
		gradc2 += -alpha
			* (R[1][0] * gradPhiLambda[0] + R[1][1]
				* gradPhiLambda[1]) * multfactor;

		/**
		 * wrt rotation omega R' (s-c)
		 */
		gradomega += alpha
			* ((-Math.sin(omega) * (x - c[0]) - Math.cos(omega)
				* (y - c[1]))
				* gradPhiLambda[0] + (Math.cos(omega)
				* (x - c[0]) - Math.sin(omega) * (y - c[1]))
				* gradPhiLambda[1]) * multfactor;

	    }
	}

	if (Double.isNaN(gradalpha)) {
	    gradalpha = 0;
	    // IJ.log("grad isNAN");
	}

	if (Double.isNaN(gradc1)) {
	    gradc1 = 0;
	    // IJ.log("grad isNAN");
	}
	if (Double.isNaN(gradc2)) {
	    gradc2 = 0;
	    // IJ.log("grad isNAN");
	}
	if (Double.isNaN(gradomega)) {
	    gradomega = 0;
	    // IJ.log("grad isNAN");
	}

	return new double[] { gradalpha, gradc1, gradc2, gradomega };
    }

    /**
     * Get the Hessian Matrix
     * 
     * @param Lambda
     * @param alpha
     * @param c
     * @param omega
     * @return
     */
    public double[][] getTransformHessian(ImplicitShape2D Lambda, double alpha,
	    double[] c, double omega) {
	double d2alpha = 0;
	double d2omega = 0;
	double d2alphaomega = 0;
	double d2alphac1 = 0;
	double d2alphac2 = 0;
	double d2omegac1 = 0;
	double d2omegac2 = 0;
	double d2cc11 = 0;
	double d2cc12 = 0;
	double d2cc22 = 0;
	double lambda = this.alighmentlambda;
	/**
	 * Rotation matrix
	 */
	double[][] R = new double[][] {
		new double[] { Math.cos(omega), -Math.sin(omega) },
		new double[] { Math.sin(omega), Math.cos(omega) } };
	double mult1, mult2, mult3, mult4;
	double mult5 = 0;
	for (int x = 0; x < width; x++) {
	    for (int y = 0; y < height; y++) {

		int newx = (int) Math.round(alpha
			* (R[0][0] * (x - c[0]) + R[0][1] * (y - c[1])));
		int newy = (int) Math.round(alpha
			* (R[1][0] * (x - c[0]) + R[1][1] * (y - c[1])));

		mult1 = (delta(x, y) + Math.pow(
			heaviside(x, y) - Lambda.heaviside(newx, newy), 2))
			* Math.pow(Math.abs(Lambda.get(newx, newy)), lambda - 2)
			* lambda * (lambda - 1);
		mult2 = (delta(x, y) + Math.pow(
			heaviside(x, y) - Lambda.heaviside(newx, newy), 2))
			* Math.pow(Math.abs(Lambda.get(newx, newy)), lambda - 1)
			* lambda * 2 * Lambda.delta(newx, newy);
		mult3 = (delta(x, y) + Math.pow(
			heaviside(x, y) - Lambda.heaviside(newx, newy), 2))
			* Math.pow(Math.abs(Lambda.get(newx, newy)), lambda - 1)
			* lambda * (Lambda.in(newx, newy) ? 1 : -1);
		mult4 = mult3;

		mult5 = -2 * lambda * Lambda.delta(newx, newy)
			* (this.heaviside(x, y) - Lambda.heaviside(newx, newy))
			* (Lambda.in(newx, newy) ? 1 : -1);

		double[] gradphiLambda = Lambda.getGradient(newx, newy);
		double[][] hes = Lambda.HessianMatrix(newx, newy);

		double[] Rprimesc = new double[] {
			-Math.sin(omega) * (x - c[0]) - Math.cos(omega)
				* (y - c[1]),
			Math.cos(omega) * (x - c[0]) - Math.sin(omega)
				* (y - c[1]) };

		/**
		 * first and second terms are outer products wrt first
		 * derivatives
		 */

		// alpha alpha
		d2alpha += (mult1 + mult2)
			* Math.pow(
				((newx / alpha) * gradphiLambda[0] + (newy / alpha)
					* gradphiLambda[1]), 2);
		// omega omega
		d2omega += (mult1 + mult2)
			* alpha
			* alpha
			* Math.pow(Rprimesc[0] * gradphiLambda[0] + Rprimesc[1]
				* gradphiLambda[1], 2);

		// The derivatives wrt c are an outer product of Rgradphi. The
		// negatives cancel

		// ( alpha R gradphi )^2 <- outer product

		d2cc11 += (mult1 + mult2)
			* alpha
			* alpha
			* Math.pow(R[0][0] * gradphiLambda[0] + R[0][1]
				* gradphiLambda[1], 2);
		d2cc12 += (mult1 + mult2)
			* alpha
			* alpha
			* (R[0][0] * gradphiLambda[0] + R[0][1]
				* gradphiLambda[1])
			* (R[1][0] * gradphiLambda[0] + R[1][1]
				* gradphiLambda[1]);

		d2cc22 += (mult1 + mult2)
			* alpha
			* alpha
			* Math.pow(R[1][0] * gradphiLambda[0] + R[1][1]
				* gradphiLambda[1], 2);

		// alpha - omega

		d2alphaomega += (mult1 + mult2)
			* alpha
			* ((newx / alpha) * gradphiLambda[0] + (newy / alpha)
				* gradphiLambda[1])
			* (Rprimesc[0] * gradphiLambda[0] + Rprimesc[1]
				* gradphiLambda[1]);

		// alpha c

		d2alphac1 += -(mult1 + mult2)
			* (newx * gradphiLambda[0] + newy * gradphiLambda[1])
			* (R[0][0] * gradphiLambda[0] + R[0][1]
				* gradphiLambda[1]);
		d2alphac2 += -(mult1 + mult2)
			* (newx * gradphiLambda[0] + newy * gradphiLambda[1])
			* (R[1][0] * gradphiLambda[0] + R[1][1]
				* gradphiLambda[1]);

		// omega c

		d2omegac1 += (mult1 + mult2)
			* -alpha
			* alpha
			* (Rprimesc[0] * gradphiLambda[0] + Rprimesc[1]
				* gradphiLambda[1])
			* (R[0][0] * gradphiLambda[0] + R[0][1]
				* gradphiLambda[1]);

		d2omegac2 += (mult1 + mult2)
			* -alpha
			* alpha
			* (Rprimesc[0] * gradphiLambda[0] + Rprimesc[1]
				* gradphiLambda[1])
			* (R[1][0] * gradphiLambda[0] + R[1][1]
				* gradphiLambda[1]);

		/**
		 * Third term involves the hessian matrix of phi_Lambda
		 */

		d2alpha += mult3
			* (hes[0][0] * newx / alpha + hes[0][1] * newy / alpha)
			* newx / alpha
			+ (hes[1][0] * newx / alpha + hes[1][1] * newy / alpha)
			* newy / alpha;
		d2omega += mult3
			* alpha
			* alpha
			* ((hes[0][0] * Rprimesc[0] + hes[0][1] * Rprimesc[1])
				* Rprimesc[0] + (hes[1][0] * Rprimesc[0] + hes[1][1]
				* Rprimesc[1])
				* Rprimesc[1]);

		d2alphaomega += mult3
			* ((hes[0][0] * Rprimesc[0] + hes[0][1] * Rprimesc[1])
				* newx + (hes[1][0] * Rprimesc[0] + hes[1][1]
				* Rprimesc[1])
				* newy);

		d2alphac1 += -mult3
			* (hes[0][0] * (R[0][0] * newx + R[0][1] * newy) + hes[0][1]
				* (R[1][0] * newx + R[1][1] * newy));
		d2alphac2 += -mult3
			* (hes[1][0] * (R[0][0] * newx + R[0][1] * newy) + hes[1][1]
				* (R[1][0] * newx + R[1][1] * newy));

		d2omegac1 += -mult3
			* alpha
			* alpha
			* ((hes[0][0] * Rprimesc[0] + hes[0][1] * Rprimesc[1])
				* R[0][0] + (hes[1][0] * Rprimesc[0] + hes[1][1]
				* Rprimesc[1])
				* R[0][1]);
		d2omegac2 += -mult3
			* alpha
			* alpha
			* ((hes[0][0] * Rprimesc[0] + hes[0][1] * Rprimesc[1])
				* R[0][1] + (hes[1][0] * Rprimesc[0] + hes[1][1]
				* Rprimesc[1])
				* R[1][1]);

		d2cc11 += mult3
			* alpha
			* alpha
			* (R[0][0]
				* (hes[0][0] * R[0][0] + hes[0][1] * R[1][0]) + R[0][1]
				* (hes[1][0] * R[0][0] + hes[1][1] * R[1][0]));
		d2cc12 += mult3
			* alpha
			* alpha
			* (R[0][0]
				* (hes[0][0] * R[0][0] + hes[0][1] * R[0][1]) + R[0][1]
				* (hes[1][0] * R[0][1] + hes[1][1] * R[1][1]));
		d2cc22 += mult3
			* alpha
			* alpha
			* (R[1][0]
				* (hes[0][0] * R[0][0] + hes[0][1] * R[0][1]) + R[1][1]
				* (hes[1][0] * R[0][1] + hes[1][1] * R[1][1]));

		/**
		 * Fourth term
		 */

		// d2alpha +=0; // no contribution
		// d2cc += 0 no contribution

		d2omega += -(newx * gradphiLambda[0] + newy * gradphiLambda[1])
			* mult4;
		d2alphaomega += mult4
			* (Rprimesc[0] * gradphiLambda[0] + Rprimesc[1]
				* gradphiLambda[1]);

		d2alphac1 += -mult4
			* (R[0][0] * gradphiLambda[0] + R[0][1]
				* gradphiLambda[1]);
		d2alphac2 += -mult4
			* (R[1][0] * gradphiLambda[0] + R[1][1]
				* gradphiLambda[1]);

		d2omegac1 += -mult4
			* alpha
			* (-Math.sin(omega) * gradphiLambda[0] - Math
				.cos(omega) * gradphiLambda[1]);
		d2omegac2 += -mult4
			* alpha
			* (Math.cos(omega) * gradphiLambda[0] - Math.sin(omega)
				* gradphiLambda[1]);

	    }
	}

	// IJ.log("Hessian at: " + alpha + " " + c[0] + " " + c[1] + " " +
	// omega);
	// IJ.log("           ["+d2alpha + " " + d2alphac1 + " " + d2alphac2 +
	// " " + d2alphaomega+"]");
	// IJ.log("           ["+d2alphac1 + " " + d2cc11 + " " + d2cc12 + " " +
	// d2omegac1+"]");
	// IJ.log("           ["+d2alphac2 + " " + d2cc12 + " " + d2alphac2 +
	// " " + d2omegac2+"]");
	// IJ.log("           ["+d2alphaomega + " " + d2omegac1 + " " +
	// d2omegac2 + " " + d2omega+"]");

	double[][] hes = new double[][] {
		new double[] { d2alpha, d2alphac1, d2alphac2, d2alphaomega },
		new double[] { d2alphac1, d2cc11, d2cc12, d2omegac1 },
		new double[] { d2alphac2, d2cc12, d2cc22, d2omegac2 },
		new double[] { d2alphaomega, d2omegac1, d2omegac2, d2omega } };

	for (int i = 0; i < 4; i++) {
	    for (int j = 0; j < 4; j++) {
		hes[i][j] = Double.isNaN(hes[i][j]) ? Double.MAX_VALUE
			: hes[i][j];
	    }
	}

	/**
	 * Hessian is 4x4 matrix
	 */
	return hes;
    }

    /**
     * Inverse of the Hessian matrix
     * 
     * @param Lambda
     * @param alpha
     * @param c
     * @param omega
     * @return
     */
    public double[][] getInverseTransformHessian(ImplicitShape2D Lambda,
	    double alpha, double[] c, double omega) {
	double[][] hes = getTransformHessian(Lambda, alpha, c, omega);

	double det = hes[0][3] * hes[1][2] * hes[2][1] * hes[3][0] - hes[0][2]
		* hes[1][3] * hes[2][1] * hes[3][0] - hes[0][3] * hes[1][1]
		* hes[2][2] * hes[3][0] + hes[0][1] * hes[1][3] * hes[2][2]
		* hes[3][0] + hes[0][2] * hes[1][1] * hes[2][3] * hes[3][0]
		- hes[0][1] * hes[1][2] * hes[2][3] * hes[3][0] - hes[0][3]
		* hes[1][2] * hes[2][0] * hes[3][1] + hes[0][2] * hes[1][3]
		* hes[2][0] * hes[3][1] + hes[0][3] * hes[1][0] * hes[2][2]
		* hes[3][1] - hes[0][0] * hes[1][3] * hes[2][2] * hes[3][1]
		- hes[0][2] * hes[1][0] * hes[2][3] * hes[3][1] + hes[0][0]
		* hes[1][2] * hes[2][3] * hes[3][1] + hes[0][3] * hes[1][1]
		* hes[2][0] * hes[3][2] - hes[0][1] * hes[1][3] * hes[2][0]
		* hes[3][2] - hes[0][3] * hes[1][0] * hes[2][1] * hes[3][2]
		+ hes[0][0] * hes[1][3] * hes[2][1] * hes[3][2] + hes[0][1]
		* hes[1][0] * hes[2][3] * hes[3][2] - hes[0][0] * hes[1][1]
		* hes[2][3] * hes[3][2] - hes[0][2] * hes[1][1] * hes[2][0]
		* hes[3][3] + hes[0][1] * hes[1][2] * hes[2][0] * hes[3][3]
		+ hes[0][2] * hes[1][0] * hes[2][1] * hes[3][3] - hes[0][0]
		* hes[1][2] * hes[2][1] * hes[3][3] - hes[0][1] * hes[1][0]
		* hes[2][2] * hes[3][3] + hes[0][0] * hes[1][1] * hes[2][2]
		* hes[3][3];

	// if(det==0){
	// IJ.log(""+hes[0][0] + " " + hes[0][1] + " " + hes[0][2] + " " +
	// hes[0][3]);
	// IJ.log(""+hes[1][0] + " " + hes[1][1] + " " + hes[1][2] + " " +
	// hes[1][3]);
	// IJ.log(""+hes[2][0] + " " + hes[2][1] + " " + hes[2][2] + " " +
	// hes[2][3]);
	// IJ.log(""+hes[3][0] + " " + hes[3][1] + " " + hes[3][2] + " " +
	// hes[3][3]);
	// }

	double b00 = hes[1][2] * hes[2][3] * hes[3][1] - hes[1][3] * hes[2][2]
		* hes[3][1] + hes[1][3] * hes[2][1] * hes[3][2] - hes[1][1]
		* hes[2][3] * hes[3][2] - hes[1][2] * hes[2][1] * hes[3][3]
		+ hes[1][1] * hes[2][2] * hes[3][3];
	double b01 = hes[0][3] * hes[2][2] * hes[3][1] - hes[0][2] * hes[2][3]
		* hes[3][1] - hes[0][3] * hes[2][1] * hes[3][2] + hes[0][1]
		* hes[2][3] * hes[3][2] + hes[0][2] * hes[2][1] * hes[3][3]
		- hes[0][1] * hes[2][2] * hes[3][3];
	double b02 = hes[0][2] * hes[1][3] * hes[3][1] - hes[0][3] * hes[1][2]
		* hes[3][1] + hes[0][3] * hes[1][1] * hes[3][2] - hes[0][1]
		* hes[1][3] * hes[3][2] - hes[0][2] * hes[1][1] * hes[3][3]
		+ hes[0][1] * hes[1][2] * hes[3][3];
	double b03 = hes[0][3] * hes[1][2] * hes[2][1] - hes[0][2] * hes[1][3]
		* hes[2][1] - hes[0][3] * hes[1][1] * hes[2][2] + hes[0][1]
		* hes[1][3] * hes[2][2] + hes[0][2] * hes[1][1] * hes[2][3]
		- hes[0][1] * hes[1][2] * hes[2][3];
	double b11 = hes[0][2] * hes[2][3] * hes[3][0] - hes[0][3] * hes[2][2]
		* hes[3][0] + hes[0][3] * hes[2][0] * hes[3][2] - hes[0][0]
		* hes[2][3] * hes[3][2] - hes[0][2] * hes[2][0] * hes[3][3]
		+ hes[0][0] * hes[2][2] * hes[3][3];
	double b12 = hes[0][3] * hes[1][2] * hes[3][0] - hes[0][2] * hes[1][3]
		* hes[3][0] - hes[0][3] * hes[1][0] * hes[3][2] + hes[0][0]
		* hes[1][3] * hes[3][2] + hes[0][2] * hes[1][0] * hes[3][3]
		- hes[0][0] * hes[1][2] * hes[3][3];
	double b13 = hes[0][2] * hes[1][3] * hes[2][0] - hes[0][3] * hes[1][2]
		* hes[2][0] + hes[0][3] * hes[1][0] * hes[2][2] - hes[0][0]
		* hes[1][3] * hes[2][2] - hes[0][2] * hes[1][0] * hes[2][3]
		+ hes[0][0] * hes[1][2] * hes[2][3];
	double b22 = hes[0][1] * hes[1][3] * hes[3][0] - hes[0][3] * hes[1][1]
		* hes[3][0] + hes[0][3] * hes[1][0] * hes[3][1] - hes[0][0]
		* hes[1][3] * hes[3][1] - hes[0][1] * hes[1][0] * hes[3][3]
		+ hes[0][0] * hes[1][1] * hes[3][3];
	double b23 = hes[0][3] * hes[1][1] * hes[2][0] - hes[0][1] * hes[1][3]
		* hes[2][0] - hes[0][3] * hes[1][0] * hes[2][1] + hes[0][0]
		* hes[1][3] * hes[2][1] + hes[0][1] * hes[1][0] * hes[2][3]
		- hes[0][0] * hes[1][1] * hes[2][3];
	double b33 = hes[0][1] * hes[1][2] * hes[2][0] - hes[0][2] * hes[1][1]
		* hes[2][0] + hes[0][2] * hes[1][0] * hes[2][1] - hes[0][0]
		* hes[1][2] * hes[2][1] - hes[0][1] * hes[1][0] * hes[2][2]
		+ hes[0][0] * hes[1][1] * hes[2][2];

	// If Hessian is zero, we are also at a fixed point, except when lambda
	// = 1!!!! wtf!!!! if Determinant is zero, return zero matrix

	if (det < 0) {
	    // IJ.log("Determinant " + det);
	    return new double[][] { new double[] { 0, 0, 0, 0 },
		    new double[] { 0, 0, 0, 0 }, new double[] { 0, 0, 0, 0 },
		    new double[] { 0, 0, 0, 0 } };
	}

	return new double[][] {
		new double[] { b00 / det, b01 / det, b02 / det, b03 / det },
		new double[] { b01 / det, b11 / det, b12 / det, b13 / det },
		new double[] { b02 / det, b12 / det, b22 / det, b23 / det },
		new double[] { b03 / det, b13 / det, b23 / det, b33 / det } };

    }

    public double getMeanEdgeCurvature() {
	double s = 0;
	double norm = 0;
	for (int x = 0; x < width; x++) {
	    for (int y = 0; y < height; y++) {
		s += this.getCurvature(x, y) * this.delta(x, y);
		norm += this.delta(x, y);
	    }

	}
	return s / norm;
    }

    public double getMaxEdgeCurvature() {
	double max = Double.MIN_VALUE;
	for (int x = 0; x < width; x++) {
	    for (int y = 0; y < height; y++) {
		if (delta(x, y) > 0)
		    max = getCurvature(x, y) > max ? getCurvature(x, y) : max;
	    }
	}
	return max;
    }

    /**
     * Find the location with the largest value in
     * 
     * @return (x,y)
     */
    public int[] getPeak() {
	int maxx = 0;
	int maxy = 0;
	for (int x = 0; x < this.signedDistance.length; x++) {
	    for (int y = 0; y < this.signedDistance[0].length; y++) {
		if (this.in(x, y)) {
		    if (get(x, y) > get(maxx, maxy)) {
			maxx = x;
			maxy = y;
		    }
		}
	    }
	}
	return new int[] { maxx, maxy };
    }

}
