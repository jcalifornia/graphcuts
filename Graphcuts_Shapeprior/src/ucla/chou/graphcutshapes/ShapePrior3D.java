package ucla.chou.graphcutshapes;

import ij.IJ;

import java.util.ArrayList;

public class ShapePrior3D {
    public ArrayList<ImplicitShape3D> shapes;
    double tauSquared;
    public double[] weights;
    public double scale; // the universal inertialScale for the prior
    public double[] orientation; // the universal inertialOrientation for the prior
    public double[] center;

    ShapePrior3D(ArrayList<ImplicitShape3D> shapes) {
	this(shapes,null);
	
    }
    
    ShapePrior3D(ArrayList<ImplicitShape3D> shapes, double[] weights){
	this.shapes = shapes;
	if(weights==null){
	    int N = shapes.size();
	    weights = new double[N];
	    for(int i = 0; i<N; i++) weights[i]=1.0/N;
	}
	else this.weights = weights;
	normalizeWeights();
	IJ.log("Constructing the shape prior using "+shapes.size()+" shapes");
	updateTauSquared();
	
    }

    /*
     * Just add another template into the shape prior
     */
    public void addTemplate(ImplicitShape3D shape) {
	// we assume shape is already normalized
	//shape.normalize();
	this.shapes.add(shape);
	updateTauSquared();
    }
    
    
    
    public double getTauSquared(){
	return this.tauSquared;
    }

    public double computeDistance(ImplicitShape3D priorShape,
	    ImplicitShape3D implicitShape3D) {
	return priorShape.distance(implicitShape3D);
    }
    
    public void updateTauSquared(){
	// We assume that the shapes are already normalized and set to the correct inertialScale/inertialOrientation!
	/**
         * Go through shapes, calculate distances
         */
	IJ.log("Updating tau^2");
        double minDistances[] = new double[this.shapes.size()];
        double tempDistance;
        double sigmaS = 0;
        for (int j = 0; j < this.shapes.size(); j++) {
            IJ.showProgress(j+1, shapes.size());
            minDistances[j] = Double.MAX_VALUE;
            for (int k = 0; k < this.shapes.size(); k++) {
                if (k == j)
                    continue;
                else {
                    tempDistance = computeDistance(shapes.get(k), shapes.get(j));
                    if (tempDistance < minDistances[j])
                        minDistances[j] = tempDistance;
                }
            }
            sigmaS += this.weights[j] * minDistances[j];
            IJ.log(""+sigmaS);
        }

        this.tauSquared = 1/sigmaS;
	
    }

    
    public void normalizeWeights(){
	double sum=0;
	for(int i =0;i<weights.length;i++){
	    sum+=weights[i];
	}
	for(int i=0;i<weights.length;i++){
	    weights[i]/=sum;
	}
    }

}
