package beast.core;

import java.util.List;
import java.util.Random;

import beast.core.parameter.RealParameter;
import beast.evolution.tree.TraitSet;
import beast.evolution.tree.Tree;
import beast.evolution.tree.coalescent.BayesianSkyline;
//import beast.core.parameter.RealParameterList;
import beast.evolution.tree.coalescent.Coalescent;
import beast.evolution.tree.coalescent.ExponentialGrowth;
import beast.evolution.tree.coalescent.IntervalList;
import beast.evolution.tree.coalescent.PopulationFunction;

//package populationPrior;

public class populationPrior extends Distribution {
	final public Input<Coalescent> coalescentPopulationModelInput = new Input<>("coalescentPopulationModel",
            "the coalescent population model");
	
	final public Input<BayesianSkyline> skylineInput = new Input<>("skyline",
            "the bayesian skyline population model");
	
	
	final public Input<RealParameter> casesInput = new Input<>("cases",
            "list of number of observed cases");
	
	final public Input<RealParameter> datesInput = new Input<>("dates",
            "list of dates for observed cases");
	
	final public Input<Tree> treeInput =  new Input<>("tree",
            "the tree for maxdate");
	
	final public Input<RealParameter> ascertainmentRateInput = new Input<>("ascertainmentRateInput",
            "ascertainment rate from real to observed cases");
	
	final public Input<RealParameter> maxDateInput = new Input<>("maxDate",
            "discard the cases later than this date");
    
	public double ascertainmentRate = 0;
	
	int sourceModel = 0;
	public BayesianSkyline srcSkyline = null;
	public PopulationFunction srcPopFunction = null;

	double sigma = 0;
int day0 = 30+25;  // January 25 (latest sequence in the dataset)
//double ascertainmentRate = 0.15;

	public double getPopSize(double t) {
		if (sourceModel == 0) {
			if (coalescentPopulationModelInput.get() != null) {
				srcPopFunction = (PopulationFunction)coalescentPopulationModelInput.get().popSizeInput;
				sourceModel = 1;
			}
			if (skylineInput.get() != null) {
				srcSkyline = skylineInput.get();
				sourceModel = 2;
			}
		}
		if (sourceModel == 1) return srcPopFunction.getPopSize(t);
		if (sourceModel == 2) return srcSkyline.getPopSize(t);
		return 0;
	}
	
	@Override 
	public double calculateLogP() {
		//Coalescent coal = coalescentPopulationModelInput.get();
		//PopulationFunction popsizeFunction = coal.popSizeInput.get();
		Tree tree = treeInput.get();
		TraitSet ts = tree.getDateTrait();
		if (sigma == 0) {
			ascertainmentRate = ascertainmentRateInput.get().getValue();
			sigma = Math.sqrt(ascertainmentRate*(1-ascertainmentRate)*(1-ascertainmentRate)+ascertainmentRate*ascertainmentRate*(1-ascertainmentRate));
			
		}
		
		//double t0 = tree.getDateTrait().maxValue;
		double t0 = ts.getDate(0);
		double logL = 0;
		double alpha = 365.0/5.0/5.0;// alphaInput.get().getArrayValue();
		
		for (int j = 0; j < casesInput.get().getDimension() ; j++) { // we consider the epidemic curve only From Dec 1 to Jan 11
			double cases = casesInput.get().getArrayValue(j);
			double date = datesInput.get().getArrayValue(j);
			
			if (maxDateInput.get().getValue() != 0 && date > maxDateInput.get().getValue()) continue;
			
			double dt = t0 - date;
			//double popSize = popsizeFunction.getPopSize(dt);
			double popSize = getPopSize(dt);
			double modeledCovariate = popSize * alpha; 
			double normalizedDiff = Math.abs(cases - modeledCovariate*ascertainmentRate) / Math.sqrt(modeledCovariate) / sigma;
			// normalizedDiff follows approximately a standard normal law
			logL -= normalizedDiff * 2 / Math.sqrt(2*Math.PI); // 1st order approximation
		}				
		return logL;
	}

	@Override
	public List<String> getArguments() {
		// TODO Auto-generated method stub
		return null;
	}

	@Override
	public List<String> getConditions() {
		// TODO Auto-generated method stub
		return null;
	}

	@Override
	public void sample(State state, Random random) {
		// TODO Auto-generated method stub
		
	}
}
