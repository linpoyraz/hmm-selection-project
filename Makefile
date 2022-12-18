simulate_scenario1:
	cd simulate && Rscript wf_selection_sim.R scenario_1.txt 10000 0.1
	$(info ---Simulated files for Scenario 1 are saved to ./simulate/simulated_files/)

simulate_scenario2:
	cd simulate && Rscript wf_selection_sim.R scenario_2.txt 10000 0.1
	$(info ---Simulated files for Scenario 1 are saved to ./simulate/simulated_files/)

analyze_scenario1:
	Rscript hmm_single_allele.R ./simulate/simulated_files/Ne_10000_init_freq_0.1_scenario_1.txt Ne_10000_init_freq_0.1_scenario_1
	$(info ---The plot for Scenario 1 is saved to ./results/single_allele/)
	$(info ---This should take around a minute, each selection coefficient set is printed)

analyze_scenario2:
	Rscript hmm_single_allele.R ./simulate/simulated_files/Ne_10000_init_freq_0.1_scenario_2.txt Ne_10000_init_freq_0.1_scenario_2
	$(info ---The plot for Scenario 2 is saved to ./results/single_allele/)
	$(info ---This should take around 4 minutes, each selection coefficient set is printed)

single_allele:
	make simulate_scenario1
	make simulate_scenario2
	make analyze_scenario1
	make analyze_scenario2

multi_allele:
	$(warning ---This takes close to 4 hours to complete)
	Rscript hmm_scan.R

clean:
	rm ./results/single_allele/*
	rm ./simulate/simulated_files/*
	$(info ---Removed simulated data and removed generated plots)
	

