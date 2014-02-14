package ucla.chou.graphcutshapes;

import ij.ImagePlus; 
import ij.gui.Roi;
import ij.plugin.filter.PlugInFilter;
import ij.process.ImageProcessor;

public class SuperEllipse_Fitter implements PlugInFilter {
    Roi r;
    boolean[][] inputMask;
    ImagePlus imp;
    public int setup(String arg, ImagePlus imp) {
	this.r=imp.getRoi();
	inputMask = new boolean[imp.getWidth()][imp.getHeight()];
	for(int x=0;x<imp.getWidth();x++){
	    for(int y=0;y<imp.getHeight();y++){
		inputMask[x][y]=r.contains(x, y);
	    }
	}
	this.imp = imp;
	return DOES_ALL+ROI_REQUIRED;
    }

    public void run(ImageProcessor ip) {
	SuperEllipseFitter ef = new SuperEllipseFitter(inputMask);
	ef.fit();
	Roi fitted = ef.getRoi();
	//ip.draw(fitted);
	this.imp.setRoi(fitted);
	this.imp.updateAndDraw();
 
    }


}