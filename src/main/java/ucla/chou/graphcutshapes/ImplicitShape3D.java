package ucla.chou.graphcutshapes;

import java.awt.Polygon;

import ij.IJ;
import ij.ImagePlus;
import ij.gui.PolygonRoi;
import ij.gui.Roi;
import ij.gui.Wand;
import ij.process.FloatProcessor;
import ij.process.ImageProcessor;

/*
 * This class describes a shape object stored implicitly as
 * a signed euclidean distance function
 */

public class ImplicitShape3D {
    /**
     * [slice][x][y]
     */
    public float[][][] signedDistance;
    public boolean[][][] mask;

    // keep this just in case we need it
    public float[][][] standardScaledDistance; // A copy of the signed distance
					       // that has been scaled to a
					       // standard mesh (VGA?)
    double scale;
    double[] orientation;
    double[] center;
    final double lambda = 1;
    final int d = 3;
    int width;
    int height;
    int slices;

    private double epsilon=0.5;

    public ImplicitShape3D(boolean[][][] mask) {
	IJ.log("Computing signed distance in 3D");
	this.signedDistance = masktoSignedDistance(mask);
	this.width = mask.length;
	this.height = mask[0].length;
	normalize();

	this.mask=mask;
	// Save a normalized copy of the shape
	// width = 320, height = 240 half-VGA

	// &TODO

    }

    /*
     * Calculate inertialCenter of mass
     */
    double[] centerofMass() {
	int d = signedDistance.length;
	int w = signedDistance[0].length;
	int h = signedDistance[0][0].length;
	double xc = 0;
	double yc = 0;
	double zc = 0;
	double normal = 0;

	for (int z = 0; z < slices; z++) {
	    for (int i = 0; i < w; i++) {
		for (int j = 0; j < h; j++) {
		    if (signedDistance[z][i][j] >= 0) {
			xc += signedDistance[z][i][j] * i;
			yc += signedDistance[z][i][j] * j;
			zc += signedDistance[z][i][j] * z;
			normal += signedDistance[z][i][j];
		    }
		}
	    }
	}
	double[] center = { xc / normal, yc / normal, zc / normal };
	return center;
    }

    /*
     * Calculate inertialOrientation TODO
     */
    double[] orientation() {
	double[] orientation = new double[2];
	// Compute the second moment about the inertialCenter
	int d = signedDistance.length;
	int w = signedDistance[0].length;
	int h = signedDistance[0][0].length;
	double xc = this.center[0];
	double yc = this.center[1];
	double zc = this.center[2];
	double q1 = 0;
	double q2 = 0;
	// double q3=0;
	double q4 = 0;
	for (int z = 0; z < slices; z++) {
	    for (int x = 0; x < w; x++) {
		for (int y = 0; y < h; y++) {
		    if (signedDistance[z][x][y] >= 0) {
			q1 += (x - xc)
				* (x - xc)
				* Math.pow(Math.abs(signedDistance[z][x][y]),
					this.lambda);
			q4 += (y - yc)
				* (y - yc)
				* Math.pow(Math.abs(signedDistance[z][x][y]),
					this.lambda);
			q2 += (x - xc)
				* (y - yc)
				* Math.pow(Math.abs(signedDistance[z][x][y]),
					this.lambda);
		    }
		}
	    }
	}

	double Tr = q1 + q4;
	double det = q1 * q4 - q2 * q2;

	double l1 = 0.5 * (Tr + Math.sqrt(Tr * Tr - 4 * det));
	// double l2 = 0.5 * (Tr - Math.sqrt(Tr * Tr - 4 * det));

	double norm = Math.sqrt(q2 * q2 + (q1 - l1) * (q1 - l1));

	orientation[0] = (q2) / norm;
	orientation[1] = (q1 - l1) / norm;

	// The inertialOrientation is the null space
	return orientation;
	// Should also set the sign of this vector but worry about that later
    }

    /*
     * Calculate inertialScale
     */
    double scale() {
	double sum = 0;
	for (int z = 0; z < slices; z++) {
	    for (int x = 0; x < this.signedDistance.length; x++) {
		for (int y = 0; y < this.signedDistance[0].length; y++) {
		    sum += this.signedDistance[z][x][y] > 0 ? Math
			    .pow(Math.abs(this.signedDistance[z][x][y]),
				    this.lambda) : 0;
		}
	    }
	}
	return Math.pow(sum, 1.0 / (d + lambda));

    }

    // re-inertialCenter and re-align to the reference shape

    /**
     * Rescale the current implicit shape to match the reference shape
     * 
     * @param reference
     */
    public void rescale(ImplicitShape3D reference) {
	ImplicitShape3D rescaled = scaleTo(reference);
	this.signedDistance = rescaled.signedDistance;
	this.width = rescaled.width;
	this.height = rescaled.height;
	this.slices = rescaled.slices;
	this.center = rescaled.center;
	this.scale = rescaled.scale;
	this.orientation = rescaled.orientation;

    }

    /**
     * Return a shape that is the current shape rescaled to the reference shape
     * 
     * @param reference
     * @return
     */
    ImplicitShape3D scaleTo(ImplicitShape3D reference) {

	int w = reference.getWidth();
	int h = reference.getHeight();
	int slices = reference.getSlices();
	reference.getScale();
	boolean mask[][][] = new boolean[w][h][slices];

	double r, rold;
	double theta = Math.acos(this.orientation[0] * reference.orientation[0]
		+ this.orientation[1] * reference.orientation[1]);
	double xold;
	double yold;
	for (int x = 0; x < w; x++) {
	    for (int y = 0; y < h; y++) {
		// This is how far this point is from the desired reference
		// frame
		r = Math.sqrt(Math.pow(x - reference.center[0], 2)
			+ Math.pow(y - reference.center[0], 2));
		rold = Math.pow(this.scale / reference.scale * r,
			1 / this.lambda);
	    }
	}

	return new ImplicitShape3D((mask));

    }

    private int getSlices() {
	return this.slices;
    }

    private double getScale() {
	return (this.scale);
    }

    boolean[][][] getMask() {
	int l = this.signedDistance.length;
	int w = this.signedDistance[0].length;
	boolean[][][] mask = new boolean[l][w][slices];
	for (int z = 0; z < slices; z++) {
	    for (int x = 0; x < l; x++) {
		for (int y = 0; y < w; y++) {
		    mask[z][x][y] = in(x, y, z);
		}
	    }
	}
	return mask;

    }

    double distance(ImplicitShape3D reference) {
	// compute d(self,reference)
	ImplicitShape3D transformed = scaleTo(reference);

	int w = this.getWidth();
	int h = this.getHeight();
	int slices = this.getSlices();
	double d = 0;
	for (int z = 0; z < slices; z++) {
	    for (int x = 0; x < w; x++) {
		for (int y = 0; y < h; y++) {
		    if (reference.get(x, y, z) * this.get(x, y, z) < 0) {
			d += Math.pow(Math.abs(reference.get(x, y, z)), lambda);
		    }

		}
	    }
	}

	return d;

    }

    void normalize() {
	this.center = centerofMass();
	this.scale = scale();
	this.orientation = orientation();

    }

    /**
     * 
     * @param mask
     *            [slice][x][y]
     * @return
     */
    public float[][][] masktoSignedDistance(boolean[][][] mask) {

	int slices = mask.length;
	int width = mask[0].length;
	int height = mask[0][0].length;
	float[][][] grid = new float[slices][width][height];
	for (int slice = 0; slice < slices; slice++) {
	    for (int x = 0; x < width; x++) {
		for (int y = 0; y < height; y++) {
		    if (mask[slice][x][y])
			grid[slice][x][y] = -Float.MAX_VALUE;
		    else
			grid[slice][x][y] = Float.MAX_VALUE;
		}

	    }
	}
	IJ.log("Passing grid to signed-distance transformation");
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
     * @param grid
     *            [slice][x][y]
     * @return
     */
    public static float[][][] ExactDistanceTransform(float[][][] grid) {
	// Restore the signed distance

	int depth = grid.length;
	int width = grid[0].length;
	int height = grid[0][0].length;

	float[][][] newgrid = new float[depth][width][height];
	//IJ.log("still ok");
	float[] c = new float[width];
	float[] r = new float[height];
	float[] p = new float[depth];
	IJ.log("Traversing slices");
	for (int s = 0; s < depth; s++) {
	    IJ.showProgress(s, depth);
	    for (int x = 0; x < width; x++) {
		for (int y = 0; y < height; y++) {
		    r[y] = grid[s][x][y] < 0 ? -Float.MAX_VALUE
			    : Float.MAX_VALUE;
		}
		float[] d1 = EDTTransform1D(r);
		for (int y = 0; y < height; y++) {
		    newgrid[s][x][y] = d1[y];
		}
	    }
	//    IJ.log("still ok 1 "+s);
	}
	for (int s = 0; s < depth; s++) {
	   
	    for (int y = 0; y < height; y++) {
		for (int x = 0; x < width; x++) {
		    c[x] = newgrid[s][x][y];
		}
		float[] d2 = EDTTransform1D(c);
		for (int x = 0; x < grid.length; x++) {
		    newgrid[s][x][y] = d2[x];
		}
	    }
	  //  IJ.log("still ok 2"+s);
	}
	for (int x = 0; x < width; x++) {
	    IJ.showProgress(x, width);
	    for (int y = 0; y < height; y++) {

		for (int s = 0; s < depth; s++) {
		    p[s] = newgrid[s][x][y];
		}

		float[] d3 = EDTTransform1D(p);
		for (int z = 0; z < depth; z++) {
		    newgrid[z][x][y] = d3[z];
		}
	    }
	}

	for (int z = 0; z < depth; z++)
	    for (int x = 0; x < width; x++) {
		for (int y = 0; y < height; y++) {
		    newgrid[z][x][y] = (float) (newgrid[z][x][y] < 0 ? 0.5f - Math
			    .sqrt(-newgrid[z][x][y]) : -0.5f
			    + Math.sqrt(newgrid[z][x][y]));
		}
	    }

	IJ.log("Completed distance transformation");
	return newgrid;

    }

    /**
     * Get the value of the signed distance at the specified coordinates
     * 
     * @param x
     * @param y
     * @return
     */
    public double get(int x, int y, int z) {
	// Should really do some interpolation here
	return this.signedDistance[Math.min(z, slices)][Math.max(0,
		Math.min(x, signedDistance.length - 1))][Math.max(0,
		Math.min(y, signedDistance[0].length - 1))];
    }

    int getWidth() {
	return width;
    }

    int getHeight() {
	return height;
    }

    boolean in(int x, int y, int z) {
	return this.signedDistance[z][x][y] <= 0 ? true : false;
    }

    /**
     * 
     * @return FloatProcessor with pixels set to the signed distance
     */
    public FloatProcessor getLSFloatProcessor(int z) {
	FloatProcessor ip2 = new FloatProcessor(width, height);
	for (int x = 0; x < width; x++) {
	    for (int y = 0; y < height; y++) {
		ip2.setf(x, y, signedDistance[z][x][y]);
	    }
	}
	ip2.resetMinAndMax();
	ip2.crop();
	return (ip2);
    }

    public ImagePlus getLSImagePlus() {
	FloatProcessor ip2 = getLSFloatProcessor(1);
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
    public PolygonRoi getRoi(int slice) {
	// Returns the ROI corresponding to the zero level set
	ImageProcessor ip2 = this.getLSFloatProcessor(slice);
	int minpt[] = { 0, 0 };

	for (int x = 0; x < width; x++) {
	    for (int y = 0; y < height; y++) {
		if (signedDistance[slice][x][y] < signedDistance[slice][minpt[0]][minpt[1]]) {
		    minpt[0] = x;
		    minpt[1] = y;
		}
	    }
	}
	Wand w = new Wand(ip2);
	w.autoOutline(minpt[0], minpt[1],
		signedDistance[slice][minpt[0]][minpt[1]], 0,
		Wand.EIGHT_CONNECTED);
	if (w.npoints == 0) {
	    return null;
	}

	Polygon p = new Polygon(w.xpoints, w.ypoints, w.npoints);
	PolygonRoi poly = new PolygonRoi(p, Roi.POLYGON);
	return poly;
    }

    
    public double delta(int x, int y,int z){
	if(Math.abs(get(x,y,z))>1) return 0;
	return(1/Math.PI*this.epsilon/(this.epsilon*this.epsilon+get(x,y,z)*get(x,y,z)));
    }
    
    public double heaviside(int x, int y,int z){
	if(get(x,y,z)<-1) return 0;
	else if(get(x,y,z)>1) return 1;
	else return(0.5*Math.atan(get(x,y,z)/this.epsilon)/Math.PI);
    }
     
    public double[] getGradient(int x, int y, int z){
	double dphiX = (get(x + 1, y,z) - get(x - 1, y,z)) / 2;
	double dphiY = (get(x, y + 1,z) - get(x, y - 1,z)) / 2;
	double dphiZ = (get(x,y,z+1) - get(x,y,z-1))/2;
	return new double[] { dphiX, dphiY, dphiZ };
    }
    
    public double[][] HessianMatrix(int x, int y, int z){
	//double[] grad0 = getGradient(x,y,z);
	
	double[][] hes = new double[3][3];
	hes[0][0] = get(x+1,y,z)-2*get(x,y,z)+get(x-1,y,z);
	hes[0][1] = (get(x+1,y+1,z)-get(x+1,y-1,z)-get(x-1,y+1,z)+get(x-1,y-1,z))/4;
	hes[0][2] = (get(x+1,y,z+1)-get(x+1,y,z-1)-get(x-1,y,z+1)+get(x-1,y,z-1))/4;
	hes[1][0] = hes[0][1];
	hes[1][1] = get(x,y+1,z)-2*get(x,y,z)+get(x,y-1,z);
	hes[1][2] =  (get(x,y+1,z+1)-get(x,y+1,z-1)-get(x,y-1,z+1)+get(x,y-1,z-1))/4;
	hes[2][0] = hes[0][2];
	hes[2][1] = hes[1][2];
	hes[2][2] = get(x,y,z+1)-2*get(x,y,z) +get(x,y,z-1);
	
	
	return hes;
    }
    
    public double gaussianCurvature(int x, int y, int z){
	double[][] hes = this.HessianMatrix(x, y, z);
	double[] grad = this.getGradient(x, y, z);
	// need determinant of 3x3 matrix
	double curv = hes[0][0]*(hes[1][1]*hes[2][2]-hes[1][2]*hes[2][1])
		-hes[0][1]*(hes[1][0]*hes[2][2]-hes[1][2]*hes[2][0])
		+hes[0][2]*(hes[1][0]*hes[2][1]-hes[1][1]*hes[2][0]);
	return curv*Math.pow(norm(grad),-1.5);
	
    }
    
    public double meanCurvature(int x,int y,int z){
	double[][] hes = this.HessianMatrix(x, y, z);
	double[] grad = this.getGradient(x, y, z);
	double curv = hes[0][0]*(grad[1]*grad[1]+grad[2]*grad[2]) + hes[1][1]*(grad[0]*grad[0]
		+grad[2]*grad[2])+hes[2][2]*(grad[0]*grad[0]+grad[1]*grad[1])-2*hes[0][1]*grad[0]*grad[1]
			-2*hes[0][2]*grad[0]*grad[2]-2*hes[1][2]*grad[1]*grad[2];
	return curv*Math.pow(norm(grad),-1.5);
    }
    
    public double norm(double[] v){
	int s=0;
	for(int i=0;i<v.length;i++){
	    s+=v[i]*v[i];
	}
	return Math.sqrt(s);
    }
}
