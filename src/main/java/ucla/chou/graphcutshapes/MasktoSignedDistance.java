package ucla.chou.graphcutshapes;

import ij.IJ;
import ij.ImagePlus;
import ij.plugin.filter.PlugInFilter;
import ij.process.FloatProcessor;
import ij.process.ImageProcessor;

public class MasktoSignedDistance implements PlugInFilter {
    String arg;
    ImagePlus imp;
    ImplicitShape2D shape;
    public void run(ImageProcessor ip) {
	imp.unlock();
	boolean[][] mask = new boolean[ip.getWidth()][ip.getHeight()];
	    for(int x=0;x<mask.length;x++){
		for(int y=0;y<mask[0].length;y++){
		    mask[x][y] = ip.getPixelValue(x, y)>0 ? true : false;
		}
	    }
	    	
	try{
	    shape = new ImplicitShape2D(mask);

	    
	} catch (Exception e){
	    IJ.showStatus("Something wrong happened");
	}
	FloatProcessor fp = new FloatProcessor(shape.width,shape.height);
	for(int x=0;x<shape.width;x++){
	    for(int y=0;y<shape.height;y++){
		fp.putPixelValue(x, y, shape.get(x, y));
	    }
	}
	//fp.multiply(-1);
	ImagePlus distanceImage = new ImagePlus("Shape embedding",fp);

	distanceImage.show();
	distanceImage.updateAndDraw();
	
    }
    public int setup(String arg0, ImagePlus imp) {
	this.arg=arg0;
	this.imp = imp;
	return DOES_ALL;
    }
    
    
}