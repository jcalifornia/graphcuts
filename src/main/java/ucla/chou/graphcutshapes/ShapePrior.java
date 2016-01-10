package ucla.chou.graphcutshapes;

import ij.IJ;

import java.util.ArrayList;

public class ShapePrior {
	public ArrayList<ImplicitShape2D> shapes;
	public double beta;
	double multiplier = 1.0;
	public double[] priorWeights = null;
	public double inertialScale; // the universal inertialScale for the prior

	public ShapePrior(ArrayList<ImplicitShape2D> shapes) {
		this(shapes, null);

	}

	public ShapePrior(ArrayList<ImplicitShape2D> shapes, double[] priorWeights) {
		IJ.log("Constructing the shape prior using " + shapes.size()
				+ " shapes");
		this.shapes = shapes;
		if (priorWeights == null) {
			int N = shapes.size();
			this.priorWeights = new double[N];
			for (int i = 0; i < N; i++) {
				this.priorWeights[i] = 1.0 / N;
			}
		} else
			this.priorWeights = priorWeights;
		normalizeWeights();
		updateBeta();
		this.inertialScale = shapes.get(0).inertialScale;
		IJ.log("Shape prior precision " + this.beta);

	}

	public void setMultiplier(double value) {
		this.multiplier = value;
	}

	/*
     * Just add another template into the shape prior
     */
	public void addTemplate(ImplicitShape2D shape) {
		// we assume shape is already normalized
		this.shapes.add(shape);
		updateBeta();
	}

	/*

    \beta is defined in the JMIV paper
     */
	public double getBeta(boolean mult) {
		//return 1;
		return mult ? this.beta * this.multiplier : this.beta;
	}

	public double computeDistance(ImplicitShape2D Omega, ImplicitShape2D Lambda) {
		return Omega.untransformedDistance(Lambda);
	}

	public void updateBeta() {
		// We assume that the shapes are already normalized and set to the
		// correct inertialScale/inertialOrientation!
		/**
		 * Go through shapes, calculate distances
		 */
		IJ.log("Updating $\\beta$");
		double minDistances[] = new double[this.shapes.size()];
		double tempDistance;
		double sigmaS = 0;
		for (int j = 0; j < this.shapes.size(); j++) {
			IJ.showStatus("Updating $\\beta$");
			IJ.showProgress(j + 1, shapes.size());
			minDistances[j] = Double.MAX_VALUE;
			for (int k = 0; k < this.shapes.size(); k++) {
				if (k != j) {
					tempDistance = computeDistance(shapes.get(k), shapes.get(j));
					if (tempDistance < minDistances[j])
						minDistances[j] = tempDistance;
				}
			}
			sigmaS += this.priorWeights[j] * minDistances[j];
			IJ.log("Sigma " + sigmaS);
		}

		this.beta = 1.0 / sigmaS;

	}

	public ImplicitShape2D getPeakShape() {
		int maxIndex = 0;
		{
			for (int i = 0; i < this.shapes.size(); i++) {
				if (this.priorWeights[i] > this.priorWeights[maxIndex]) {
					maxIndex = i;
				}
			}
		}
		IJ.log("Template "+maxIndex+" is the most probable shape");
		return this.shapes.get(maxIndex);
	}

	public void normalizeWeights() {
		double sum = 0;
		for (int i = 0; i < priorWeights.length; i++) {
			sum += priorWeights[i];
		}
		for (int i = 0; i < priorWeights.length; i++) {
			priorWeights[i] /= sum;
		}
	}

}