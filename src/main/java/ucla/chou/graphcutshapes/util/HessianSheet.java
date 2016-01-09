package ucla.chou.graphcutshapes.util;

import ij.ImagePlus;
import ij.plugin.filter.PlugInFilter;
import ij.process.ImageProcessor;

public class HessianSheet implements PlugInFilter{

    int zDepth;
    int yDepth;
    int xDepth;
    ImagePlus imp;
    String arg;
    public void run(ImageProcessor arg0) {
	// Get the smoothness
	
	
    }

    public int setup(String arg0, ImagePlus arg1) {
	this.imp = arg1;
	this.arg=arg0;
	this.zDepth= arg1.getStackSize();
	this.yDepth = arg1.getHeight();
	this.xDepth=arg1.getWidth();
	return STACK_REQUIRED;
    }
    
}