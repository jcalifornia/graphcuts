package ucla.chou.graphcutshapes;

import fiji.util.gui.OverlayedImageCanvas;
import ij.IJ;
import ij.ImagePlus;
import ij.ImageStack;
import ij.WindowManager;
import ij.gui.Overlay;
import ij.gui.Roi;
import ij.gui.StackWindow;
import ij.plugin.PlugIn;
import ij.plugin.frame.RoiManager;
import ij.process.AutoThresholder;
import ij.process.FloatProcessor;
import ij.process.ImageProcessor;
import ij.process.LUT;

import java.awt.Color;
import java.awt.Component;
import java.awt.Dimension;
import java.awt.Graphics;
import java.awt.GridBagConstraints;
import java.awt.GridBagLayout;
import java.awt.Insets;
import java.awt.List;
import java.awt.Panel;
import java.awt.Rectangle;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.ComponentAdapter;
import java.awt.event.ComponentEvent;
import java.awt.event.ItemEvent;
import java.awt.event.ItemListener;
import java.awt.event.KeyAdapter;
import java.awt.event.KeyEvent;
import java.awt.event.KeyListener;
import java.awt.event.WindowAdapter;
import java.awt.event.WindowEvent;
import java.io.File;
import java.util.ArrayList;
import java.util.Vector;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;

import javax.swing.BorderFactory;
import javax.swing.JButton;
import javax.swing.JFileChooser;
import javax.swing.JPanel;
import javax.swing.SwingUtilities;

import ucla.chou.graphcutshapes.util.ImageOverlay;

public class GraphcutsShapeprior_Plugin implements PlugIn {

    /** maximum number of classes (labels) allowed on the GUI */
    private static final int MAX_NUM_TEMPLATES = 256;
    /** array of lists of Rois for each class */
    private Vector<Roi> templates = new Vector<Roi>();
    private Vector<Roi> segmentations = new Vector<Roi>();

    final double betaMultiplier = 0.5;
    final double lambda = 2;
    final double alignmentlambda = 2;

    /*
     * Segmentation objects
     */

    ShapePrior shapep;
    GraphCutSegmenter graphSegmenter;
    IntensityModel likelihood;

    ImplicitShape2D segmentation;

    /*
     * Image objects
     */

    // ImagePlus imageToSegment;
    ImagePlus trainingTemplates;

    ImagePlus displayImage;
    int slice = 1;

    int numTemplates = 0;

    /* GUI objects */
    RoiManager rm;
    private GuiWindow win;
    JButton loadImageButton;
    JButton loadTrainingButton;
    JButton segmentImageButton;
    JButton overlayResultButton;
    JButton optionsButton;

    JButton deleteShapeButton;

    /** available colors for available classes */
    final Color[] colors = new Color[] { Color.red, Color.green, Color.blue,
	    Color.cyan, Color.magenta };

    /** array of roi list overlays to paint the transparent rois of each class */
    RoiListOverlay[] roiOverlay;

    private java.awt.List templateList; // training templates
    private java.awt.List segmentationList; // segmentations

    /** array of buttons for adding each trace class */
    private JButton[] addTemplateButton;

    final ExecutorService exec = Executors.newFixedThreadPool(1);
    ImageOverlay resultOverlay;

    private boolean imageLoaded = false;

    LUT overlayLUT;

    public GraphcutsShapeprior_Plugin() {
	final byte[] red = new byte[256];
	final byte[] green = new byte[256];
	final byte[] blue = new byte[256];
	overlayLUT = new LUT(red, green, blue);

	loadImageButton = new JButton("Load image to segment");
	loadImageButton.setToolTipText("Opens a file selector");

	loadTrainingButton = new JButton("Load training templates");
	loadTrainingButton.setToolTipText("Select binary masked stack");

	segmentImageButton = new JButton("Segment image");
	segmentImageButton.setToolTipText("Segment the current image");
	// segmentImageButton.setEnabled(false);

	overlayResultButton = new JButton("Overlay the results");
	overlayResultButton.setToolTipText("Overlay segmentation results");

	optionsButton = new JButton("Options");
	optionsButton.setToolTipText("Advanced options");

	deleteShapeButton = new JButton("Delete selected");
	deleteShapeButton.setToolTipText("");

	templateList = new List(10, true);
	templateList.setSize(10, 10);

	segmentationList = new List(10, true);
	segmentationList.setSize(10, 10);

	this.likelihood = new LaplaceMixture(); // or Gaussian

    }

    public void run(String arg0) {
	if (WindowManager.getCurrentImage() != null) {
	    displayImage = (ImagePlus) WindowManager.getCurrentImage()
		    .duplicate();
	    this.imageLoaded = true;
	} else {
	    displayImage = new ImagePlus("Empty Image", new FloatProcessor(640,
		    480));
	}
	if (Math.max(displayImage.getWidth(), displayImage.getHeight()) > 1024)
	    if (!IJ.showMessageWithCancel("Warning",
		    "At least one dimension of the image \n"
			    + "is larger than 1024 pixels. \n"
			    + "Consider downsizing the image first \n"
			    + "Proceed?"))
		return;
	// displayImage = new ImagePlus();
	// displayImage.setProcessor("Trainable Segmentation", imageToSegment
	// .getProcessor().duplicate());

	SwingUtilities.invokeLater(new Runnable() {
	    public void run() {
		win = new GuiWindow(displayImage);
		win.pack();

	    }
	});
	IJ.log("Welcome...");
	 this.rm = RoiManager.getInstance();
	 if (rm != null)
	            rm.setVisible(true);
	        else
	            rm = new RoiManager();

	// Show the Roi Manager is it is not already showing
	
	 
    }

    public void getTrainingShapes() {
	// this.shapep = null;
	File imageFile;
	JFileChooser fileChooser = new JFileChooser();
	fileChooser.setFileSelectionMode(JFileChooser.FILES_ONLY);
	fileChooser.setDialogTitle("Open a stack of training masks");
	fileChooser.setMultiSelectionEnabled(false);

	int returnVal = fileChooser.showOpenDialog(null);
	if (returnVal == JFileChooser.APPROVE_OPTION) {
	    imageFile = fileChooser.getSelectedFile();
	} else {
	    return;
	}
	String file = imageFile.getPath();
	ImagePlus newImage = IJ.openImage(file);
	int slices = newImage.getStackSize();

	ArrayList<ImplicitShape2D> shapes = new ArrayList<ImplicitShape2D>(
		slices);
	AutoThresholder at = new AutoThresholder();

	ImageStack shapeStack = new ImageStack(newImage.getWidth(),
		newImage.getHeight());
	ImplicitShape2D userShape;
	if (displayImage.getRoi() != null) {
	    boolean[][] mask2 = new boolean[displayImage.getWidth()][displayImage
		    .getHeight()];
	    for (int y = 0; y < displayImage.getHeight(); y++) {
		for (int x = 0; x < displayImage.getWidth(); x++) {
		    mask2[x][y] = displayImage.getRoi().contains(x, y);
		}
	    }
	    userShape = new ImplicitShape2D(mask2);
	    IJ.log("User shape has inertialScale " + userShape.inertialScale
		    + " inertialOrientation "
		    + userShape.inertialOrientation[0] + ","
		    + userShape.inertialOrientation[1]);
	} else
	    userShape = null;
	for (int s = 1; s <= slices; s++) {
	    IJ.showStatus("Reading shapes...");
	    int threshold = -1;
	    ImageProcessor sliceip = newImage.getStack()
		    .getProcessor(s <= slices ? s : slices)
		    .convertToByte(false);
	    int[] data = (sliceip.getHistogram());

	    int minbin = -1, maxbin = -1;
	    for (int i = 0; i < data.length; i++) {
		if (data[i] > 0)
		    maxbin = i;
	    }
	    for (int i = data.length - 1; i >= 0; i--) {
		if (data[i] > 0)
		    minbin = i;
	    }
	    // IJ.log (""+minbin+" "+maxbin);
	    int[] data2 = new int[(maxbin - minbin) + 1];
	    for (int i = minbin; i <= maxbin; i++) {
		data2[i - minbin] = data[i];
		;
	    }

	    // Apply the selected algorithm
	    if (data2.length < 2) {
		threshold = 0;
	    }
	    threshold = at.getThreshold("Huang", data2);
	    threshold += minbin;
	    // go along one of the edge of the image and decide

	    int gl = 0;
	    int ll = 0;
	    for (int y = 0; y < sliceip.getHeight(); y++) {
		if (sliceip.getPixelValue(0, y) > threshold)
		    gl++;
		else
		    ll++;
	    }

	    boolean glin = gl < ll;

	    boolean[][] mask = new boolean[sliceip.getWidth()][sliceip
		    .getHeight()];
	    for (int y = 0; y < sliceip.getHeight(); y++) {
		for (int x = 0; x < sliceip.getWidth(); x++) {

		    mask[x][y] = sliceip.getPixelValue(x, y) > threshold ? glin
			    : !glin;
		}

	    }

	    ImplicitShape2D shape = new ImplicitShape2D(mask);

	    // If we already have an active ROI, normalize the shapes to the ROI

	    if (userShape != null) {
		// shape.translateInertially(userShape);

		double[] trans = shape.alignByNewton(userShape);
		shape.affineTransform(1 / trans[0], new double[] { -trans[1],
			-trans[2] }, -trans[3], userShape.width,
			userShape.height);
	    } else {
		if (shapes.size() > 0) {
		    // shape.translateInertially(shapes.get(0)); // get close!
		    double[] transformParams = shapes.get(0).alignByNewton(
			    shape);
		    shape.affineTransform(transformParams[0], new double[] {
			    transformParams[1], transformParams[2] },
			    transformParams[3]);
		}
	    }
	    shapes.add(shape);
	    ImageProcessor currentShapeProcessor = shape.getMaskProcessor();

	    templates.add(shape.getRoi(false));
	    shapeStack.addSlice(currentShapeProcessor);
	    IJ.showProgress(s, slices);

	}

	ImagePlus distanceImage = new ImagePlus("Training templates",
		shapeStack);
	distanceImage.show();
	distanceImage.updateAndDraw();

	// IJ.log("Constructing shape prior...");
	this.numTemplates = shapes.size();
	this.shapep = new ShapePrior(shapes);
	this.shapep.setMultiplier(this.betaMultiplier);
	// this.shapep.updateTauSquared();
	this.segmentImageButton.setEnabled(true);
	IJ.log("Completed construction of shape prior object");

    }

    /**
     * Action for button that loads an image to segment
     */
    public void loadImage() {

	java.io.File imageFile;
	JFileChooser fileChooser = new JFileChooser();
	fileChooser.setFileSelectionMode(JFileChooser.FILES_ONLY);
	fileChooser.setMultiSelectionEnabled(false);

	int oldImageHeight = displayImage.getHeight();
	int oldImageWidth = displayImage.getWidth();

	int returnVal = fileChooser.showOpenDialog(null);
	if (returnVal == JFileChooser.APPROVE_OPTION) {
	    imageFile = fileChooser.getSelectedFile();
	} else {
	    return;
	}

	displayImage = null;
	win.removeScrollBars();
	displayImage = IJ.openImage(imageFile.getPath());
	// this.imageToSegment = displayImage.duplicate();
	win.setImage(displayImage);

	win.setSize(displayImage.getWidth() - oldImageWidth + win.getWidth(),
		displayImage.getHeight() - oldImageHeight + win.getHeight());
	win.validate();
	displayImage.show();
	displayImage.updateAndDraw();
	win.getCanvas().setSize(displayImage.getWidth(),
		displayImage.getHeight());

	if (displayImage.getStackSize() > 1) {
	    win.addScrollBars();
	}
	win.repaintAll();
    }

    public void getImage() {

    }

    /**
     * Perform segmentation of currently displayed image
     */
    public void segmentImage() {
	if (this.shapep == null)
	    IJ.log("Performing segmentation without shapes...");
	else
	    IJ.log("Performing image segmentation with shapes....");
	graphSegmenter = new GraphCutSegmenter(this.displayImage.getProcessor());
	if (shapep != null)
	    graphSegmenter.setShapeKernelDensityEstimate(shapep);

	Roi userRoi = displayImage.getRoi();
	// Contains sequence of segmentations that should converge to the
	// result!!
	Vector<ImplicitShape2D> segSequence = new Vector<ImplicitShape2D>();
	ImplicitShape2D currentShape;
	boolean[][] mask = new boolean[this.displayImage.getWidth()][this.displayImage
		.getHeight()];

	if (userRoi == null || shapep == null) {

	    // Initialize without user selection, by using graph cuts with a
	    // length penalty
	    IJ.log("Initializing the segmentation algorithm using clustering...");

	    ImageProcessor bytepr = ((ImageProcessor) displayImage
		    .getProcessor().duplicate()).convertToByte(false);
	    bytepr.autoThreshold();

	    /**
	     * Traverse around edge of image and set the equivalence class of
	     * the edge to OUTSIDE of the region
	     */

	    int in = 0;
	    int out = 0;

	    for (int y = 0; y < this.displayImage.getHeight(); y++) {
		if (bytepr.getPixelValue(0, y) == 0)
		    in++;
		else
		    out++;
		if (bytepr.getPixelValue(this.displayImage.getWidth() - 1, y) == 0)
		    in++;
		else
		    out++;
	    }

	    /*
	     * ImagePlus temp = new ImagePlus("",bytepr); temp.show();
	     * temp.updateAndDraw();
	     */
	    for (int y = 0; y < this.displayImage.getHeight(); y++) {
		for (int x = 0; x < this.displayImage.getWidth(); x++) {
		    mask[x][y] = bytepr.getPixelValue(x, y) == 0 ? in < out
			    : !(in < out);
		}

	    }

	    /**
	     * The thresholded shape
	     */

	    final int iters = 1;
	    double OLDENERGY = Double.MAX_VALUE;
	    for (int i = 0; i < iters; i++) {
		IJ.log("Re-inferring the image statistics...");

		likelihood.Infer(displayImage.getProcessor(), mask, false);
		IJ.log("IM: " + likelihood.getPosteriorMean(true) + " OM: "
			+ likelihood.getPosteriorMean(false) + " IB: "
			+ likelihood.getPosteriorPrecision(true));

		IJ.log("Constructing the graph...");
		graphSegmenter = new GraphCutSegmenter(
			displayImage.getProcessor());
		graphSegmenter.setIntensityModel(likelihood);
		IJ.log("Setting edgeweights...");
		graphSegmenter.setLengthPenalty(displayImage.getWidth() / 40);
		graphSegmenter.setNodeWeights(likelihood);
		graphSegmenter.setEdgeWeights();

		IJ.log("Finding the minimum cut...");
		graphSegmenter.relaxEnergy();

		IJ.log("Found minimum cut of " + graphSegmenter.Energy);
		IJ.log("Inner proportion: " + graphSegmenter.fractionInner());

		mask = graphSegmenter.returnMask();
		OLDENERGY = graphSegmenter.Energy;
		currentShape = graphSegmenter.getLevelSet();

		displayImage.updateAndDraw();
		if (Math.abs(OLDENERGY - graphSegmenter.Energy) < 1)
		    break;
	    }

	    currentShape = graphSegmenter.getLevelSet();
	    if (graphSegmenter.fractionInner() == 0
		    || graphSegmenter.fractionInner() == 1) {
		IJ.log("Collapse of segmentation");
		return;
	    }
	} else {
	    for (int y = 0; y < this.displayImage.getHeight(); y++) {
		for (int x = 0; x < this.displayImage.getWidth(); x++) {
		    mask[x][y] = userRoi.contains(x, y);
		}
	    }
	    currentShape = new ImplicitShape2D(mask);
	}

	segSequence.add(currentShape);

	likelihood.Infer(this.displayImage.getProcessor(), mask);

	if (this.shapep == null) {
	    // We already have the answer
	    this.segmentation = currentShape;

	    return;
	}

	/**
	 * The inertial method sucks!
	 */
	shapep.setReferenceFrame(currentShape);
	currentShape = shapep.getPeakShape();

	// Weakly regularize and infer
	//shapep.setMultiplier(this.betaMultiplier*0.00001);

	Overlay ov = new Overlay(currentShape.getRoi(false));
	displayImage.setOverlay(ov);
	/**
	 * \Omega^{(0)}
	 */
	segSequence.add(currentShape);
	final int repetitions = 25; // Really the maximum number of iterations
				    // we would like

	// redo the likelihood estimation
	likelihood.Infer(this.displayImage.getProcessor(),
		currentShape.getMask());
	Double OLDENERGY = Double.MAX_VALUE;
	double ENERGY = 0;
	for (int r = 0; r < repetitions; r++) {
	    IJ.log("MM iteration " + r);
	    graphSegmenter = new GraphCutSegmenter(displayImage.getProcessor());
	    graphSegmenter.setNodeWeights(shapep, currentShape, likelihood);
	    graphSegmenter.setEdgeWeights(shapep, currentShape);
	    graphSegmenter.setLengthPenalty(displayImage.getWidth() / 160/(r+1)/(r+1));
	    // graphSegmenter.setLengthPenalty( (float)
	    // (graphSegmenter.width*graphSegmenter.height/shapep.shapes.get(0).length()));
	    IJ.log("Finding the minimum cut...");
	    ENERGY = graphSegmenter.relaxEnergy();
	    IJ.log("Found minimum cut of " + graphSegmenter.Energy);
	    IJ.log("Inner proportion: " + graphSegmenter.fractionInner());
	    float prop = graphSegmenter.fractionInner();
	    if (prop == 0 || prop == 1) {
		IJ.log("Collapse of segmentation " + prop);
		return;
	    }
	    // segment, rinse, repeat
	    currentShape = graphSegmenter.getLevelSet();
		segmentations.add(currentShape.getRoi(true));
		segmentations.get(segmentations.size()-1).setName("iter_"+r);
		rm.addRoi(segmentations.get(segmentations.size()-1));
		displayImage.setRoi(segmentations.get(segmentations.size()-1));
		
	    // do something better here
	    if (currentShape.distance(segSequence.get(segSequence.size() - 1)) > 0.00001) {
		likelihood.Infer(this.displayImage.getProcessor(),
			currentShape.getMask());
		shapep.setReferenceFrame(currentShape);

		continue;
	    }

	    // Overlay the fitted most probable shape prior here
	    // add also to Roi List

	    ov = new Overlay(shapep.getPeakShape().getRoi(false));
	    displayImage.setOverlay(ov);

	    segSequence.add(currentShape);
	    if ((OLDENERGY - ENERGY) / OLDENERGY < 0.001) {
		// break;
	    }
	    OLDENERGY = ENERGY;
	    //shapep.setMultiplier(this.betaMultiplier*r/5);
	}

	displayImage.setOverlay(null);

	
	displayImage.setRoi(currentShape.getRoi(false));
	displayImage.updateAndDraw();

    }

    /* GUI stuff below: */

    private class GuiWindow extends StackWindow implements ActionListener {
	/**
	 * 
	 */
	private static final long serialVersionUID = 8865007897342896520L;

	GridBagLayout boxAnnotation = new GridBagLayout();
	JPanel buttonsPanel = new JPanel();
	Panel all = new Panel();
	JPanel trainingJPanel = new JPanel();
	JPanel templatePanel = new JPanel();
	JPanel segmentationPanel = new JPanel();
	JPanel optionsJPanel = new JPanel();
	private CustomCanvas canvas;

	GuiWindow(ImagePlus imp) {

	    super(imp, new CustomCanvas(imp));

	    final CustomCanvas canvas = (CustomCanvas) getCanvas();
	    removeAll();
	    setTitle("MM-Graph cuts with statistical shape priors");

	    trainingJPanel.setBorder(BorderFactory
		    .createTitledBorder("Training"));
	    GridBagLayout trainingLayout = new GridBagLayout();
	    GridBagConstraints trainingConstraints = new GridBagConstraints();
	    trainingConstraints.anchor = GridBagConstraints.NORTHWEST;
	    trainingConstraints.fill = GridBagConstraints.HORIZONTAL;
	    trainingConstraints.gridwidth = 1;
	    trainingConstraints.gridheight = 1;
	    trainingConstraints.gridx = 0;
	    trainingConstraints.gridy = 0;
	    trainingConstraints.insets = new Insets(5, 5, 6, 6);
	    trainingJPanel.setLayout(trainingLayout);

	    // loadTrainingButton.setEnabled(false);
	    trainingJPanel.add(loadImageButton, trainingConstraints);
	    trainingConstraints.gridy++;
	    trainingJPanel.add(loadTrainingButton, trainingConstraints);
	    trainingConstraints.gridy++;
	    trainingJPanel.add(segmentImageButton, trainingConstraints);
	    trainingConstraints.gridy++;

	    optionsJPanel
		    .setBorder(BorderFactory.createTitledBorder("Options"));

	    GridBagLayout optionsLayout = new GridBagLayout();
	    GridBagConstraints optionsConstraints = new GridBagConstraints();
	    optionsConstraints.anchor = GridBagConstraints.NORTHWEST;
	    optionsConstraints.fill = GridBagConstraints.HORIZONTAL;
	    optionsConstraints.gridwidth = 1;
	    optionsConstraints.gridheight = 1;
	    optionsConstraints.gridx = 0;
	    optionsConstraints.gridy = 0;
	    optionsConstraints.insets = new Insets(5, 5, 6, 6);
	    optionsJPanel.setLayout(optionsLayout);

	    optionsJPanel.add(optionsButton, optionsConstraints);

	    // Buttons panel (including training and options)
	    GridBagLayout buttonsLayout = new GridBagLayout();
	    GridBagConstraints buttonsConstraints = new GridBagConstraints();
	    buttonsPanel.setLayout(buttonsLayout);
	    buttonsConstraints.anchor = GridBagConstraints.NORTHWEST;
	    buttonsConstraints.fill = GridBagConstraints.HORIZONTAL;
	    buttonsConstraints.gridwidth = 1;
	    buttonsConstraints.gridheight = 1;
	    buttonsConstraints.gridx = 0;
	    buttonsConstraints.gridy = 0;
	    buttonsPanel.add(trainingJPanel, buttonsConstraints);
	    buttonsConstraints.gridy++;
	    buttonsPanel.add(optionsJPanel, buttonsConstraints);
	    buttonsConstraints.gridy++;
	    buttonsConstraints.insets = new Insets(5, 5, 6, 6);

	    // Template list panel
	    GridBagLayout templateLayout = new GridBagLayout();
	    GridBagConstraints templateConstraints = new GridBagConstraints();
	    templatePanel.setLayout(templateLayout);

	    templateConstraints.anchor = GridBagConstraints.NORTHWEST;
	    templateConstraints.gridheight = 1;
	    templateConstraints.gridwidth = 1;
	    templateConstraints.gridx = 0;
	    templateConstraints.gridy = 0;

	    templateConstraints.fill = GridBagConstraints.HORIZONTAL;
	    templateConstraints.insets = new Insets(5, 5, 6, 6);

	    // templateConstraints.gridy++;
	    templateLayout.setConstraints(templateList, templateConstraints);
	    templatePanel.add(templateList);

	    templateConstraints.gridy++;
	    templateLayout.setConstraints(deleteShapeButton,
		    templateConstraints);
	    // templatePanel.add(deleteShapeButton);
	    templatePanel.setBorder(BorderFactory.createTitledBorder("Shapes"));
	    templatePanel.setLayout(boxAnnotation);

	    // segmentation list panel
	    GridBagLayout segmentationLayout = new GridBagLayout();
	    GridBagConstraints segmentationConstraints = new GridBagConstraints();
	    segmentationPanel.setLayout(segmentationLayout);

	    segmentationConstraints.anchor = GridBagConstraints.NORTHWEST;
	    segmentationConstraints.gridheight = 1;
	    segmentationConstraints.gridwidth = 1;
	    segmentationConstraints.gridx = 0;
	    segmentationConstraints.gridy = 0;

	    segmentationConstraints.fill = GridBagConstraints.HORIZONTAL;
	    segmentationConstraints.insets = new Insets(5, 5, 6, 6);

	    segmentationLayout.setConstraints(segmentationList,
		    segmentationConstraints);
	    segmentationPanel.add(segmentationList);

	    // Add the actions

	    loadImageButton.addActionListener(listener);
	    loadTrainingButton.addActionListener(listener);
	    segmentImageButton.addActionListener(listener);

	    GridBagLayout layout = new GridBagLayout();
	    GridBagConstraints allConstraints = new GridBagConstraints();
	    all.setLayout(layout);

	    allConstraints.anchor = GridBagConstraints.NORTHWEST;
	    allConstraints.fill = GridBagConstraints.BOTH;
	    allConstraints.gridwidth = 1;
	    allConstraints.gridheight = 1;
	    allConstraints.gridx = 0;
	    allConstraints.gridy = 0;
	    allConstraints.weightx = 0;
	    allConstraints.weighty = 0;

	    all.add(buttonsPanel, allConstraints);

	    allConstraints.gridy++;
	    allConstraints.weightx = 0;
	    allConstraints.weighty = 0;
	    // all.add(templatePanel, allConstraints);

	    allConstraints.gridy++;
	    allConstraints.weightx = 0;
	    allConstraints.weighty = 0;
	    // all.add(segmentationPanel,allConstraints);

	    // The image canvas
	    allConstraints.gridx++;
	    allConstraints.weightx = 1;
	    allConstraints.weighty = 1;
	    allConstraints.gridy = 0;
	    all.add(canvas, allConstraints);

	    allConstraints.gridx++;
	    allConstraints.anchor = GridBagConstraints.NORTHEAST;
	    allConstraints.weightx = 0;
	    allConstraints.weighty = 0;

	    GridBagLayout wingb = new GridBagLayout();
	    GridBagConstraints winc = new GridBagConstraints();
	    winc.anchor = GridBagConstraints.NORTHWEST;
	    winc.fill = GridBagConstraints.BOTH;
	    winc.weightx = 1;
	    winc.weighty = 1;
	    setLayout(wingb);
	    add(all, winc);

	    // Propagate all listeners
	    for (Component p : new Component[] { all, buttonsPanel }) {
		for (KeyListener kl : getKeyListeners()) {
		    p.addKeyListener(kl);
		}
	    }

	    addWindowListener(new WindowAdapter() {
		public void windowClosing(WindowEvent e) {
		    exec.shutdownNow();
		    loadImageButton.removeActionListener(listener);
		    loadTrainingButton.removeActionListener(listener);
		}

	    });

	    canvas.addComponentListener(new ComponentAdapter() {
		public void componentResized(ComponentEvent ce) {
		    Rectangle r = canvas.getBounds();
		    canvas.setDstDimensions(r.width, r.height);
		}
	    });
	    this.updateSliceSelector();
	    if (this.imp.getStack().getSize() > 1) {
		addScrollBars();
	    }
	}

	/**
	 * TODO WRITE THIS!!
	 */
	void addScrollBars() {

	    // IJ.log("stack");
	}

	void removeScrollBars() {
	    // IJ.log("no stack");
	}

	public void repaintAll() {
	    this.buttonsPanel.repaint();
	    this.templatePanel.repaint();
	    getCanvas().repaint();
	    this.all.repaint();
	}
    }

    private class CustomCanvas extends OverlayedImageCanvas {
	/**
	 * 
	 */
	private static final long serialVersionUID = 4865579916617812609L;

	CustomCanvas(ImagePlus imp) {
	    super(imp);
	    Dimension dim = new Dimension(Math.max(640, imp.getWidth()),
		    Math.max(480, imp.getHeight()));
	    setMinimumSize(dim);
	    setSize(dim.width, dim.height);
	    setDstDimensions(dim.width, dim.height);
	    addKeyListener(new KeyAdapter() {
		public void keyReleased(KeyEvent ke) {
		    repaint();
		}
	    });
	}

	public void setDrawingSize(int w, int h) {
	}

	public void setDstDimensions(int width, int height) {
	    super.dstWidth = width;
	    super.dstHeight = height;
	    // adjust srcRect: can it grow/shrink?
	    int w = Math.min((int) (width / magnification), imp.getWidth());
	    int h = Math.min((int) (height / magnification), imp.getHeight());
	    int x = srcRect.x;
	    if (x + w > imp.getWidth())
		x = w - imp.getWidth();
	    int y = srcRect.y;
	    if (y + h > imp.getHeight())
		y = h - imp.getHeight();
	    srcRect.setRect(x, y, w, h);
	    repaint();
	}

	public void paint(Graphics g) {
	    Rectangle srcRect = getSrcRect();
	    double mag = getMagnification();
	    int dw = (int) (srcRect.width * mag);
	    int dh = (int) (srcRect.height * mag);
	    g.setClip(0, 0, dw, dh);

	    super.paint(g);

	    int w = getWidth();
	    int h = getHeight();
	    g.setClip(0, 0, w, h);

	    // Paint away the outside
	    g.setColor(getBackground());
	    g.fillRect(dw, 0, w - dw, h);
	    g.fillRect(0, dh, w, h - dh);
	}
    }

    /**
     * Listeners
     */
    private ActionListener listener = new ActionListener() {
	public void actionPerformed(final ActionEvent e) {
	    // listen to the buttons on separate threads not to block
	    // the event dispatch thread
	    exec.submit(new Runnable() {
		public void run() {
		    if (e.getSource() == loadImageButton) {
			try {
			    loadImage();
			} catch (Exception e) {
			    e.printStackTrace();
			}
		    } else if (e.getSource() == loadTrainingButton) {
			getTrainingShapes();
		    } else if (e.getSource() == segmentImageButton) {
			segmentImage();
		    } else {
			if (e.getSource() == templateList) {
			    deleteSelected(e);
			}
		    }

		}

		private void deleteSelected(ActionEvent e) {
		    // TODO Auto-generated method stub

		}
	    });
	}

    };

    /**
     * Compute the energy of the current segmentation
     * 
     * @param shape
     * @param lik
     * @param prior
     * @return
     */
    public double getEnergy(ImplicitShape2D shape, IntensityModel lik,
	    ShapePrior prior) {

	for(int x=0;x<this.displayImage.getWidth();x++){
	    for(int y=0;y<this.displayImage.getHeight();y++){
		
	    }
	}
	return alignmentlambda;

    }

}
