############################################################################################################
###### This file is for a brief description for the structure or organization of "SignatureAnalzyer"
###### used in PCAWG7 analysis.
############################################################################################################

0. Description overall for execution and installation guide
	a. All scripts or codes are self-contained and standalone, and don't require any external dependency.
	b. All scripts or codes are writtend in R and exectued in R environments (R-3.3.3 in Mac OS X). 
	c. Some R libraries will be required - gridExtra, ggplot2, ggplots, reshape2, and grid

1. Description for directory structure 
	a. INPUT_SignatureAnalzyer - contains all inputs (lego matrix for 96 and 1536 SNV contexts, DBSs, and INDELs) 
				and the information on putative POLE, MSI, and single TMZ sample.
	b. OUTPUT_SignatureAnalzyer - all final outputs will be save in this directory.
	c. TEMPORARY_SignatureAnalzyer - all intermeidate files will be saved at this directory.
	d. OUTPUT_DEMO - all outputs from "SignatureAnalzyer.demo.R" below will be saved here.

2. Description for main scripts or codes
	a. SignatureAnalzyer.demo.R 
		a1 Demonstrating how SignatureAnalzyer extrats signatures for 35 PCAWG Biliary samples using 96 contexts.
		a2 About 20 minutes for 10 independent BayesNMF runs with tol = 1.e-07 and K = 25 (maximum signatures) 
		    on the processon 2.6 GHz Intel Cori7 in MacBook (OS X El Captian).
	b. SignatureAnalyzer.PCAWG.COMPOSITE.R - COMPOSITE signature extraction and activity assignment for 2780 PCAWG samples.
	c. SignatureAnalyzer.PCAWG.DNP.R - DBS (double-base substitution) signature extraction and activity assignment for 2780 PCAWG samples.
	d. SignatureAnalyzer.PCAWG.INDEL.R - INDEL signature extraction and activity assignment for 2780 PCAWG samples.
	e. SignatureAnalyzer.PCAWG.function.R - contains all necessary functions.

3. Description of key functions 
	a. BayesNMF.L1W.L2H 
		a1 Bayesian non-negative matrix factorization algorithm with an exponential prior for W and a half-normal prior for H
		a2 This function is used in all signature extraction in PCAWG activity.
	b. BayesNMF.L1.KL.fixed_W.Z
		a1 Adopted from the Bayesian non-ngative matrix factorization algorithm with an exponential prior for both W and H.
		a2 This function is used in the activity attribution step to select an optimal set of signatures.
	c. BayesNMF.L1.KL.fixed_W.Z.sample
		a1 Adopted from the Bayesian non-ngative matrix factorization algorithm with an exponential prior for both W and H.
		a2 This function is used for the activity attribution with a set of selected signatures.
	d. Detailed descriptions for all other functions are contained in each script or code.

4. Relevant references 
	a. Tan, VY, Fevotte, C. Automatic relevance determination in nonnegative matrix factorization with the beta-divergence. IEEE Trans Pattern Anal Mach Intell 2013;35:1592-605.
	b. Kim J, et al. Somatic ERCC2 mutations are associated with a distinct genomic signature in urothelial tumors (2016). Nat Genet 48, 600-606.
	c. Kasar, S, Kim J et al. Whole-genome sequencing reveals activation-induced cytidine deaminase signatures during indolent chronic lymphocytic leukaemia evolution. 
		Nat Commun. 6:8866 doi: 10.1038/ncomms9866 (2015).
	d. P. Polak, Kim J, L. Brounstein  et al, A mutational signature reveals alterations underlying deficient homologous recombination repair in breast cancer. 
		Nature Genetics (2017), doi:10.1038/ng.3934
	e. Haradhvala NJ, Kim J, Maruvka YE et al, Distinct mutational signatures characterize concurrent loss of polymerase proofreading and mismatch repair, 
		Nature Commun. (2018), PMID: 29717118
