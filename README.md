The new populationPrior class is called with

```
<distribution id="popPrior" spec="util.CompoundDistribution">
				<distribution id="populationPrior.t:alignment" spec="populationPrior" 
				 cases="0 0 0 0 0 0 0 0 0 0 1 4 0 0 4 2 7 3 2 12 6 11 13 12 13 7 26 11 16 24 16 73 63 65 46 78 64 71 125 117 187 168 214 281 296 404 458 572 653 667 1254 1222 1552 1937"
				 dates="2019.9164383561645 2019.9191780821918 2019.9219178082192 2019.9246575342465 2019.927397260274 2019.9301369863015 2019.9328767123288 2019.9356164383562 2019.9383561643835 2019.941095890411 2019.9438356164383 2019.9465753424658 2019.9493150684932 2019.9520547945206 2019.954794520548 2019.9575342465753 2019.9602739726026 2019.9630136986302 2019.9657534246576 2019.968493150685 2019.9712328767123 2019.9739726027397 2019.9767123287672 2019.9794520547946 2019.982191780822 2019.9849315068493 2019.9876712328767 2019.990410958904 2019.9931506849316 2019.995890410959 2019.9986301369863 2020.0013661202186 2020.0040983606557 2020.0068306010928 2020.00956284153 2020.0122950819673 2020.0150273224044 2020.0177595628415 2020.0204918032787 2020.0232240437158 2020.025956284153 2020.02868852459 2020.0314207650274 2020.0341530054645 2020.0368852459017 2020.0396174863388 2020.042349726776 2020.045081967213 2020.0478142076502 2020.0505464480875 2020.0532786885246 2020.0560109289618 2020.058743169399 2020.061475409836 2020.0642076502731 2020.0669398907103 2020.0696721311476 2020.0724043715848 2020.0751366120219 2020.077868852459 2020.0806010928961 2020.0833333333333"
				 tree="@Tree.t:consensus"
				 ascertainmentRateInput="0.15"
				 maxDate="2020.02466"
				 coalescentPopulationModel="@CoalescentExponential.t:alignment"
				 > <!-- you can replace coalescentPopulationModel=".." by skyline="..." if you are using skyline -->
				</distribution>
</distribution>```

Cases and dates are the early SARS-CoV-2 epidemic curve


Then the relevant code computing the likelihood of the observed epidemic curve given the coalescent population curve is  


```double alpha = 365.0/5.0/5.0;// alphaInput.get().getArrayValue();
		
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
		return logL;```
