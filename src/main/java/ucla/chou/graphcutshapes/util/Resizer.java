package ucla.chou.graphcutshapes.util;

import ij.IJ;
import ij.ImagePlus;
import ij.ImageStack;
import ij.measure.Calibration;
import ij.process.ImageProcessor;
import ij.process.StackProcessor;

/** This plugin implements the Edit/Crop and Image/Adjust/Size commands. */
public class Resizer {
    public static final int IN_PLACE = 16, SCALE_T = 32;

    private static boolean averageWhenDownsizing = true;
    private static int interpolationMethod = ImageProcessor.BILINEAR;
    private String[] methods = ImageProcessor.getInterpolationMethods();
    private double origWidth, origHeight;
    private ImagePlus imp;

    public Resizer(ImagePlus imp) {
	this.imp = imp;
    }

    public void resize(int newWidth, int newHeight, int newDepth) {
	ImageProcessor ip = imp.getProcessor();

	int stackSize = imp.getStackSize();

	ip.setInterpolationMethod(interpolationMethod);

	try {
	    StackProcessor sp = new StackProcessor(imp.getStack(), ip);
	    ImageStack s2 = sp.resize(newWidth, newHeight,
		    averageWhenDownsizing);
	    int newSize = s2.getSize();
	    if (s2.getWidth() > 0 && newSize > 0) {

		Calibration cal = imp.getCalibration();
		if (cal.scaled()) {
		    cal.pixelWidth *= origWidth / newWidth;
		    cal.pixelHeight *= origHeight / newHeight;
		}

		imp.setStack(null, s2);

	    }
	    if (stackSize > 1 && newSize < stackSize)
		System.out.println("ImageJ ran out of memory causing \nthe last "
			+ (stackSize - newSize) + " slices to be lost.");
	} catch (OutOfMemoryError o) {
	    IJ.outOfMemory("Resize");
	}
	this.zScale(newDepth, interpolationMethod);
	imp.changes = true;

	//zScale(newHeight, interpolationMethod + IN_PLACE + SCALE_T);

    }

    public void zScale(int newDepth, int interpolationMethod) {
	if (imp.isHyperStack())
	    zScaleHyperstack(newDepth, interpolationMethod);
	else {
	    boolean inPlace = (interpolationMethod & IN_PLACE) != 0;
	    interpolationMethod = interpolationMethod & 15;
	    int stackSize = imp.getStackSize();
	    int bitDepth = imp.getBitDepth();
	    if (newDepth <= stackSize / 2
		    && interpolationMethod == ImageProcessor.NONE)
		shrinkZ(newDepth, inPlace);
	    else
		resizeZ(newDepth, interpolationMethod);
	    ImageProcessor ip = imp.getProcessor();
	    double min = ip.getMin();
	    double max = ip.getMax();
	    if (bitDepth == 16 || bitDepth == 32)
		imp.getProcessor().setMinAndMax(min, max);
	}

	Calibration cal = imp.getCalibration();
	if (cal.scaled())
	    cal.pixelDepth *= (double) imp.getNSlices() / imp.getNSlices();
	Object info = imp.getProperty("Info");
	if (info != null)
	    imp.setProperty("Info", info);
	if (imp.isHyperStack())
	    imp.setOpenAsHyperStack(imp.isHyperStack());
    }

    private void zScaleHyperstack(int depth2, int interpolationMethod) {
	boolean inPlace = (interpolationMethod & IN_PLACE) != 0;
	boolean scaleT = (interpolationMethod & SCALE_T) != 0;
	interpolationMethod = interpolationMethod & 15;
	int channels = imp.getNChannels();
	int slices = imp.getNSlices();
	int frames = imp.getNFrames();
	int slices2 = slices;
	int frames2 = frames;
	int bitDepth = imp.getBitDepth();
	if (slices == 1 && frames > 1)
	    scaleT = true;
	if (scaleT)
	    frames2 = depth2;
	else
	    slices2 = depth2;
	double scale = (double) (depth2 - 1) / slices;
	if (scaleT)
	    scale = (double) (depth2 - 1) / frames;
	if (scale <= 0.5 && interpolationMethod == ImageProcessor.NONE)
	    shrinkHyperstack(depth2, inPlace, scaleT);
	ImageStack stack1 = imp.getStack();
	int width = stack1.getWidth();
	int height = stack1.getHeight();
	ImagePlus imp2 = IJ.createImage(imp.getTitle(), bitDepth + "-bit",
		width, height, channels * slices2 * frames2);

	imp2.setDimensions(channels, slices2, frames2);
	ImageStack stack2 = imp2.getStack();
	ImageProcessor ip = imp.getProcessor();
	int count = 0;
	if (scaleT) {
	    IJ.showStatus("T Scaling...");
	    ImageProcessor xtPlane1 = ip.createProcessor(width, frames);
	    xtPlane1.setInterpolationMethod(interpolationMethod);
	    ImageProcessor xtPlane2;
	    Object xtpixels1 = xtPlane1.getPixels();
	    int last = slices * channels * height - 1;
	    for (int z = 1; z <= slices; z++) {
		for (int c = 1; c <= channels; c++) {
		    for (int y = 0; y < height; y++) {
			IJ.showProgress(count++, last);
			for (int t = 1; t <= frames; t++) {
			    int index = imp.getStackIndex(c, z, t);
			    // IJ.log("1: "+c+" "+z+" "+t+" "+index+" "+xzPlane1);
			    Object pixels1 = stack1.getPixels(index);
			    System.arraycopy(pixels1, y * width, xtpixels1,
				    (t - 1) * width, width);
			}
			xtPlane2 = xtPlane1.resize(width, depth2,
				averageWhenDownsizing);
			Object xtpixels2 = xtPlane2.getPixels();
			for (int t = 1; t <= frames2; t++) {
			    int index = imp2.getStackIndex(c, z, t);
			    // IJ.log("2: "+c+" "+z+" "+t+" "+index+" "+xzPlane2);
			    Object pixels2 = stack2.getPixels(index);
			    System.arraycopy(xtpixels2, (t - 1) * width,
				    pixels2, y * width, width);
			}
		    }
		}
	    }
	} else {
	    IJ.showStatus("Z Scaling...");
	    ImageProcessor xzPlane1 = ip.createProcessor(width, slices);
	    xzPlane1.setInterpolationMethod(interpolationMethod);
	    ImageProcessor xzPlane2;
	    Object xypixels1 = xzPlane1.getPixels();
	    int last = frames * channels * height - 1;
	    for (int t = 1; t <= frames; t++) {
		for (int c = 1; c <= channels; c++) {
		    for (int y = 0; y < height; y++) {
			IJ.showProgress(count++, last);
			for (int z = 1; z <= slices; z++) {
			    int index = imp.getStackIndex(c, z, t);
			    Object pixels1 = stack1.getPixels(index);
			    System.arraycopy(pixels1, y * width, xypixels1,
				    (z - 1) * width, width);
			}
			xzPlane2 = xzPlane1.resize(width, depth2,
				averageWhenDownsizing);
			Object xypixels2 = xzPlane2.getPixels();
			for (int z = 1; z <= slices2; z++) {
			    int index = imp2.getStackIndex(c, z, t);
			    Object pixels2 = stack2.getPixels(index);
			    System.arraycopy(xypixels2, (z - 1) * width,
				    pixels2, y * width, width);
			}
		    }
		}
	    }
	}
	imp2.setDimensions(channels, slices2, frames2);
	imp.setStack(imp2.getStack());
    }

    private void shrinkHyperstack(int newDepth, boolean inPlace, boolean scaleT) {
	int channels = imp.getNChannels();
	int slices = imp.getNSlices();
	int frames = imp.getNFrames();
	int factor = (int) Math.round((double) slices / newDepth);
	if (scaleT)
	    factor = frames / newDepth;
	int zfactor = scaleT ? 1 : factor;
	int tfactor = scaleT ? factor : 1;
	ImageStack stack = imp.getStack();
	ImageStack stack2 = new ImageStack(imp.getWidth(), imp.getHeight());
	boolean virtual = stack.isVirtual();
	int slices2 = slices / zfactor + ((slices % zfactor) != 0 ? 1 : 0);
	int frames2 = frames / tfactor + ((frames % tfactor) != 0 ? 1 : 0);
	int n = channels * slices2 * frames2;
	int count = 1;
	for (int t = 1; t <= frames; t += tfactor) {
	    for (int z = 1; z <= slices; z += zfactor) {
		for (int c = 1; c <= channels; c++) {
		    int i = imp.getStackIndex(c, z, t);
		    IJ.showProgress(i, n);
		    ImageProcessor ip = stack.getProcessor(imp.getStackIndex(c,
			    z, t));
		    if (!inPlace)
			ip = ip.duplicate();
		    // IJ.log(count++ +"  "+i+" "+c+" "+z+" "+t);
		    stack2.addSlice(stack.getSliceLabel(i), ip);
		}
	    }
	}
	imp.setDimensions(channels, slices2, frames2);
	IJ.showProgress(1.0);
	imp.setStack(stack2);
    }

    private void shrinkZ(int newDepth, boolean inPlace) {
	ImageStack stack = imp.getStack();
	int factor = imp.getStackSize() / newDepth;
	boolean virtual = stack.isVirtual();
	int n = stack.getSize();
	ImageStack stack2 = new ImageStack(stack.getWidth(), stack.getHeight());
	for (int i = 1; i <= n; i += factor) {
	    if (virtual)
		IJ.showProgress(i, n);
	    ImageProcessor ip2 = stack.getProcessor(i);
	    if (!inPlace)
		ip2 = ip2.duplicate();
	    stack2.addSlice(stack.getSliceLabel(i), ip2);
	}
	imp.setStack(stack2);
    }

    private void resizeZ(int newDepth, int interpolationMethod) {
	ImageStack stack1 = imp.getStack();
	int width = stack1.getWidth();
	int height = stack1.getHeight();
	int depth = stack1.getSize(); 
	int bitDepth = imp.getBitDepth();

	 ImagePlus imp2 = IJ.createImage(imp.getTitle(), bitDepth+"-bit", width, height, newDepth);
	        if (imp2==null) return ;
	        ImageStack stack2 = imp2.getStack();
	ImageProcessor ip = imp.getProcessor();
	ImageProcessor xzPlane1 = ip.createProcessor(width, depth);
	xzPlane1.setInterpolationMethod(interpolationMethod);
	ImageProcessor xzPlane2;
	Object xypixels1 = xzPlane1.getPixels();
	System.out.println("Z Scaling...");
	for (int y = 0; y < height; y++) {
	    IJ.showProgress(y, height - 1);
	    for (int z = 0; z < depth; z++) {
		Object pixels1 = stack1.getPixels(z + 1);
		System.arraycopy(pixels1, y * width, xypixels1, z * width,
			width);
	    }
	    xzPlane2 = xzPlane1.resize(width, newDepth, averageWhenDownsizing);
	    Object xypixels2 = xzPlane2.getPixels();
	    for (int z = 0; z < newDepth; z++) {
		Object pixels2 = stack2.getPixels(z+1);
		System.arraycopy(xypixels2, z * width, pixels2, y * width,
			width);
	    }
	}
	imp.setStack(stack2);
    }

}
