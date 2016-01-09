package ucla.chou.graphcutshapes;

import java.awt.Color;
import java.awt.Component;
import java.awt.Dimension;
import java.awt.Graphics;
import java.awt.GridBagConstraints;
import java.awt.GridBagLayout;
import java.awt.Insets;
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
import java.util.ArrayList;
import java.util.List;
import java.util.Vector;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;

import javax.swing.BorderFactory;
import javax.swing.JButton;
import javax.swing.JFileChooser;
import javax.swing.JList;
import javax.swing.JPanel;
import javax.swing.SwingUtilities;

import ucla.chou.graphcutshapes.util.ImageOverlay;

import java.io.File;

import fiji.util.gui.OverlayedImageCanvas;

import ij.IJ;
import ij.ImagePlus;
import ij.ImageStack;
import ij.WindowManager;
import ij.gui.ImageWindow;
import ij.gui.Roi;
import ij.gui.StackWindow;
import ij.io.Opener;
import ij.plugin.PlugIn;
import ij.process.AutoThresholder;
import ij.process.ByteProcessor;
import ij.process.FloatProcessor;
import ij.process.ImageProcessor;
import ij.process.ImageStatistics;
import ij.process.LUT;
import ij.process.StackStatistics;

public class GraphcutsShapeprior3D_Plugin implements PlugIn {

    /** maximum number of classes (labels) allowed on the GUI */
    private static final int MAX_NUM_TEMPLATES = 256;
    /** array of lists of Rois for each class */
    private Vector<Roi> templates = new Vector<Roi>();

    /*
     * Segmentation objects
     */

    ShapePrior3D shapep;
    GraphCutSegmenter3D graphSegmenter;
    IntensityModel likelihood;

    /*
     * Image objects
     */

    ImagePlus trainingTemplates;
    ImagePlus displayImage;

    int numTemplates = 0;

    /* GUI objects */
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

    private java.awt.List templateList[];

    private JList templateListBox;

    /** array of buttons for adding each trace class */
    private JButton[] addTemplateButton;

    final ExecutorService exec = Executors.newFixedThreadPool(1);
    ImageOverlay resultOverlay;

    private boolean imageLoaded = false;

    LUT overlayLUT;

    public GraphcutsShapeprior3D_Plugin() {
	final byte[] red = new byte[256];
	final byte[] green = new byte[256];
	final byte[] blue = new byte[256];
	overlayLUT = new LUT(red, green, blue);
	templateList = new java.awt.List[MAX_NUM_TEMPLATES];

	loadImageButton = new JButton("Load stack to segment");
	loadImageButton.setToolTipText("Opens a file selector");

	loadTrainingButton = new JButton("Load stack as a traning template");
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

	templateListBox = new JList(templates);

	this.likelihood = new GaussianMixture();

    }

    public void run(String arg0) {
	if (WindowManager.getCurrentImage() != null) {
	    displayImage = (ImagePlus) WindowManager.getCurrentImage().clone();
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

	SwingUtilities.invokeLater(new Runnable() {
	    public void run() {
		win = new GuiWindow(displayImage);
		win.pack();
	    }
	});
    }

    /**
     * TODO Fix this for 3d shapes
     */
    public void getTrainingShapes() {
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

	ArrayList<ImplicitShape3D> shapes = new ArrayList<ImplicitShape3D>(
		slices);
	AutoThresholder at = new AutoThresholder();

	for (int s = 1; s <= slices; s++) {
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
	    boolean[][][] mask = new boolean[newImage.getStackSize()][sliceip
		    .getWidth()][sliceip.getHeight()];

	    for (int slice = 0; slice < newImage.getStackSize(); slice++) {
		for (int x = 0; x < sliceip.getWidth(); x++) {
		    for (int y = 0; y < sliceip.getHeight(); y++) {

			mask[slice][x][y] = newImage.getStack()
				.getProcessor(slice).getPixelValue(x, y) > threshold ? true
				: false;
		    }

		}
	    }

	    ImplicitShape3D shape = new ImplicitShape3D(mask);
	    shapes.add(shape);
	    // IJ.log(""+shape.scale());
	    IJ.log("" + shape.orientation[0] + "," + shape.orientation[1]);
	    templates.add(shape.getRoi(0));
	    // signedDistanceStack.addSlice(new FloatProcessor
	    // shape.getLSFloatProcessor());
	    
	    templateListBox.add(fileChooser, "shape "+templates.size(),templates.size());
	    

	}

	IJ.log("constructing shape prior...");

	// createShapePrior(shapes);
	this.segmentImageButton.setEnabled(true);
	IJ.log("" + shapes.get(0).scale());
	IJ.log("" + shapes.get(0).center[0]);
	ImagePlus test = new ImagePlus();

	//

    }

    public void createShapePrior(ArrayList<ImplicitShape3D> shapes) {
	this.numTemplates = shapes.size();
	this.shapep = new ShapePrior3D(shapes);

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

	win.repaintAll();
	this.imageLoaded = true;
	this.segmentImageButton.setEnabled(true);
    }

    public void getImage() {

    }

    /**
     * Perform segmentation of currently displayed image
     */
    public void segmentImage() {
	IJ.log("Performing stack segmentation....");
	IJ.log("slices: " + displayImage.getStackSize());
	graphSegmenter = new GraphCutSegmenter3D(this.displayImage.getStack());
	likelihood = new LaplaceMixture();
	if (shapep != null)
	    graphSegmenter.setShapeKernelDensityEstimate(shapep);

	// Contains sequence of segmentations that should converge to the
	// result!!
	Vector<ImplicitShape3D> segSequence = new Vector<ImplicitShape3D>();

	Roi userSelection = this.displayImage.getRoi();

	boolean[][][] mask = new boolean[this.displayImage.getStackSize()][this.displayImage
		.getWidth()][this.displayImage.getHeight()];
	IJ.log("Initializing the segmentation algorithm using thresholding...");

	StackStatistics stat = new StackStatistics(displayImage);
	int[] histogram =  stat.histogram;

	AutoThresholder at = new AutoThresholder();
	
	int threshold = at.getThreshold("Huang", histogram);
	int slices = this.displayImage.getStackSize();

	for (int s = 0; s < slices; s++) {
	    IJ.showProgress(s, slices);
	    ImageProcessor copy = displayImage.getStack().getProcessor(s + 1)
		    .convertToByte(true);
	    for (int y = 0; y < this.displayImage.getHeight(); y++) {
		for (int x = 0; x < this.displayImage.getWidth(); x++) {
		    mask[s][x][y] = copy.getPixelValue(x, y) >= threshold ? true
			    : false;
		}

	    }
	}
	IJ.log("Inferring image statistics...");

	likelihood.Infer(displayImage.getStack(), mask,false);
	likelihood.multiplyPrecision(.5);
	IJ.log("IM: " + likelihood.getPosteriorMean(true) + " OM: "+likelihood.getPosteriorMean(false) + " IP: " +likelihood.getPosteriorPrecision(true));
	/**
	 * The thresholded shape
	 */
	ImplicitShape3D thresholdedShape = new ImplicitShape3D(mask);
	IJ.log("Finished creating 3d shape");

	graphSegmenter.setIntensityModel(likelihood);
	graphSegmenter.setLengthPenalty(displayImage.getWidth()/28);
	graphSegmenter.setNodeWeights(likelihood);
	IJ.log("Done setting t-link weights");
	graphSegmenter.setEdgeWeights();
	IJ.log("Done setting n-link weights");
	IJ.log("Finding the minimum cut...");
	float minE = graphSegmenter.relaxEnergy();
	IJ.log("Found a minimum cut of " + minE);
	IJ.log("Inner: "+graphSegmenter.getVolume());
	mask = graphSegmenter.returnMask();

	final int iters = 0;
	for (int i = 0; i < iters; i++) {
	    IJ.log("Re-inferring the image statistics...");
	    likelihood.Infer(displayImage.getStack(), mask,false);
		IJ.log("IM: " + likelihood.getPosteriorMean(true) + " OM: "+likelihood.getPosteriorMean(false) + " IP: " +likelihood.getPosteriorPrecision(true));

	    graphSegmenter = new GraphCutSegmenter3D(displayImage.getStack());
	    graphSegmenter.setIntensityModel(likelihood);
	    IJ.log("Weighting the graph..");
	    graphSegmenter.setLengthPenalty(displayImage.getWidth()/32);
	    IJ.log("Setting node weights");
	    graphSegmenter.setNodeWeights(likelihood);
	    graphSegmenter.setEdgeWeights();

	    IJ.log("Finding the minimum cut...");
	    minE = graphSegmenter.relaxEnergy();

	    IJ.log("Found minimum cut of " + graphSegmenter.Energy);
	    IJ.log("Inner: "+graphSegmenter.getVolume());
	    mask = graphSegmenter.returnMask();
	}

	
	ImagePlus im = graphSegmenter.returnMaskedImage();
	im.show();
	im.updateAndDraw();

	// train the initial image statistics

	// likelihood.Infer(this.displayImage.getStack(), mask);

	// align the templates

	if (this.shapep != null) {

	} else {
	    // No shape prior, just use the length penalty
	}

	// find the "mean template" to use as an initial value

	final int repetitions = 5;

	// redo the likelihood estimation
	likelihood.Infer(this.displayImage.getStack(), mask);

	for (int r = 0; r < repetitions; r++) {
	    // set the weights

	    // segment, rinse, repeat
	}

    }

    private int[] getStackHistogram(ImageStack stack, int breaks) {
	// TODO Auto-generated method stub
	return null;
    }

    /* GUI stuff below: */

    private class GuiWindow extends StackWindow {
	/**
	 * 
	 */
	private static final long serialVersionUID = 8865007897342896520L;

	GridBagLayout boxAnnotation = new GridBagLayout();
	JPanel buttonsPanel = new JPanel();
	Panel all = new Panel();
	JPanel trainingJPanel = new JPanel();
	JPanel templatePanel = new JPanel();
	JPanel optionsJPanel = new JPanel();
	private CustomCanvas canvas;

	GuiWindow(ImagePlus imp) {

	    super(imp, new CustomCanvas(imp));

	    final CustomCanvas canvas = (CustomCanvas) getCanvas();
	    removeAll();
	    setTitle("3D-Graph cuts with statistical shape priors");

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

	    boxAnnotation
		    .setConstraints(deleteShapeButton, templateConstraints);
	    templatePanel.add(deleteShapeButton);
	    templateConstraints.gridy++;

	    // templateConstraints.gridy++;

	    templatePanel.setBorder(BorderFactory.createTitledBorder("Shapes"));
	    templatePanel.setLayout(boxAnnotation);

	    for (int i = 0; i < numTemplates; i++) {
		templateList[i].addActionListener(listener);
		templateList[i].addItemListener(itemListener);

		templateConstraints.fill = GridBagConstraints.HORIZONTAL;
		templateConstraints.insets = new Insets(5, 5, 6, 6);

	    }

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
	    all.add(templatePanel, allConstraints);

	    allConstraints.gridx++;
	    allConstraints.weightx = 1;
	    allConstraints.weighty = 1;
	    allConstraints.gridy--;
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
     * Item listener for the trace lists
     */
    private ItemListener itemListener = new ItemListener() {
	public void itemStateChanged(final ItemEvent e) {
	    exec.submit(new Runnable() {
		public void run() {
		    for (int i = 0; i < numTemplates; i++) {
			if (e.getSource() == templateList[i])
			    listSelected(e, i);
		    }
		}
	    });
	}
    };

    /**
     * Select a list and deselect the others
     * 
     * @param e
     *            item event (originated by a list)
     * @param i
     *            list index
     */
    void listSelected(final ItemEvent e, final int i) {
	// drawExamples();
	displayImage.setColor(Color.YELLOW);

	for (int j = 0; j < numTemplates; j++) {
	    if (j == i) {
		final Roi newRoi = templates.get(templateList[i]
			.getSelectedIndex());
		// Set selected trace as current ROI
		newRoi.setImage(displayImage);
		displayImage.setRoi(newRoi);
	    } else
		templateList[j].deselect(templateList[j].getSelectedIndex());
	}

	displayImage.updateAndDraw();
    }
    


}