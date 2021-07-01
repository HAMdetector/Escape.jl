var documenterSearchIndex = {"docs":
[{"location":"advanced_usage/#Advanced-usage","page":"Advanced usage","title":"Advanced usage","text":"","category":"section"},{"location":"advanced_usage/#Running-sub-models","page":"Advanced usage","title":"Running sub-models","text":"","category":"section"},{"location":"advanced_usage/","page":"Advanced usage","title":"Advanced usage","text":"By default, the run_model function uses the full HAMdetector model that includes the horseshoe prior, phylogeny information and epitope prediction. Reduced versions of the full model can be run by supplying an additional positional argument:","category":"page"},{"location":"advanced_usage/","page":"Advanced usage","title":"Advanced usage","text":"using Escape\n\ndata = HLAData(\n    alignment_file = \"/home/user/Desktop/alignment.fasta\",\n    tree_file = \"/home/user/Desktop/phylogeny.tree\"\n)\n\nresult = Escape.run_model(HLAModel{4}(), data)","category":"page"},{"location":"advanced_usage/","page":"Advanced usage","title":"Advanced usage","text":"The different versions are:","category":"page"},{"location":"advanced_usage/","page":"Advanced usage","title":"Advanced usage","text":"1: A simple logistic regression model\n2: Like model 1, but additionally with horseshoe prior\n3: Like model 2, but additionally with phylogeny information included\n4: Full HAMdetector model, like model 3 but with epitope prediction included.","category":"page"},{"location":"advanced_usage/","page":"Advanced usage","title":"Advanced usage","text":"Note that the full HAMdetector model includes all the previous models as a special case.","category":"page"},{"location":"advanced_usage/#Modifying-sampler-parameters","page":"Advanced usage","title":"Modifying sampler parameters","text":"","category":"section"},{"location":"advanced_usage/","page":"Advanced usage","title":"Advanced usage","text":"By default, 4 chains with 1000 iterations each and a warmup of 300 samples are drawn. This should be plenty for all use-cases, but can be changed as shown below:","category":"page"},{"location":"advanced_usage/","page":"Advanced usage","title":"Advanced usage","text":"using Escape\n\ndata = HLAData(\n    alignment_file = \"/home/user/Desktop/alignment.fasta\",\n    tree_file = \"/home/user/Desktop/phylogeny.tree\"\n)\n\nresult = Escape.run_model(\n    data, chains = 4, iter = 1000, warmup = 300, stan_args = \"adapt_delta=0.95\"\n)","category":"page"},{"location":"advanced_usage/","page":"Advanced usage","title":"Advanced usage","text":"For a full list of the CmdStan options, please refer to the CmdStan sampling parameters. It should not be necessary to modify the stan_args argument, it is just listed here for completeness. ","category":"page"},{"location":"advanced_usage/#Using-4-digit-HLA-alleles","page":"Advanced usage","title":"Using 4-digit HLA alleles","text":"","category":"section"},{"location":"advanced_usage/","page":"Advanced usage","title":"Advanced usage","text":"By default, HAMdetector uses HLA alleles with 2-digit accuracy. To use HLA alleles with 4-digit accuracy (if supplied), set the allele_depth keyword argument to 2 (default is 1).","category":"page"},{"location":"advanced_usage/","page":"Advanced usage","title":"Advanced usage","text":"using Escape\n\ndata = HLAData(\n    alignment_file = \"/home/user/Desktop/alignment_with_4_digit_alleles.fasta\",\n    tree_file = \"/home/user/Desktop/phylogeny.tree\",\n    allele_depth = 2\n)","category":"page"},{"location":"quick_start/#Quick-start","page":"Quick start","title":"Quick start","text":"","category":"section"},{"location":"quick_start/#Running-HAMdetector","page":"Quick start","title":"Running HAMdetector","text":"","category":"section"},{"location":"quick_start/","page":"Quick start","title":"Quick start","text":"HAMdetector requires two things: A set of aligned amino acid sequences that are annotated  with HLA alleles and a phylogenetic tree. The sequence alignment has to be provided in FASTA format, the phylogenetic tree in Newick format.","category":"page"},{"location":"quick_start/","page":"Quick start","title":"Quick start","text":"The easiest way of running HAMdetector is by providing the input files in a standardized way, that is, the sequences in the .fasta file are annotated by an underscore separated list of HLA alleles and the tips of the phylogenetic tree are labeled by their position in the multiple sequence alignment.","category":"page"},{"location":"quick_start/","page":"Quick start","title":"Quick start","text":"An example of valid input files is shown below:","category":"page"},{"location":"quick_start/","page":"Quick start","title":"Quick start","text":"alignment.fasta","category":"page"},{"location":"quick_start/","page":"Quick start","title":"Quick start","text":">HLA-A*32_HLA-A*03_HLA-B*52_HLA-B*14_HLA-C*16_HLA-C*16\nMG-RASVMG-RA*V\n>HLA-A*68_HLA-A*31_HLA-B*44_HLA-B*08_HLA-C*07_HLA-C*07\nMGARASVMG-RT*V\n>HLA-A*24_HLA-A*24_HLA-B*15_HLA-B*46_HLA-C*12_HLA-C*12\nMSARASVMGSRSSV\n>HLA-A*02_HLA-A*03_HLA-B*58_HLA-B*44_HLA-C*14_HLA-C*06\nMSARASVMG-RRSV\n>HLA-A*03_HLA-A*31_HLA-B*55_HLA-B*40_HLA-C*16_HLA-C*03\nMSARASVMG-RMSV","category":"page"},{"location":"quick_start/","page":"Quick start","title":"Quick start","text":"phylogeny.tree","category":"page"},{"location":"quick_start/","page":"Quick start","title":"Quick start","text":"((2:0.033310,1:0.118635):33.615353,(3:0.107192,5:0.857526):0.102147,4:0.247706);","category":"page"},{"location":"quick_start/","page":"Quick start","title":"Quick start","text":"Assuming these files were saved as /home/user/Desktop/alignment.fasta and /home/user/Desktop/phylogeny.tree, HAMdetector can be run with the following commands:","category":"page"},{"location":"quick_start/","page":"Quick start","title":"Quick start","text":"using Escape\n\ndata = HLAData(\n    alignment_file = \"/home/user/Desktop/alignment.fasta\",\n    tree_file = \"/home/user/Desktop/phylogeny.tree\"\n)\n\nresult = run_model(data)\n\nsummary = replacement_summary(result)","category":"page"},{"location":"quick_start/","page":"Quick start","title":"Quick start","text":"In a first step, the input data are converted into a suitable format to run the HAMdetector model. This step includes calculating phylogeny information and running MHCFlurry for epitope prediction.","category":"page"},{"location":"quick_start/","page":"Quick start","title":"Quick start","text":"In a second step, the HAMdetector model is run on the input data.","category":"page"},{"location":"quick_start/","page":"Quick start","title":"Quick start","text":"The replacement_summary function then computes a table of all observed replacements, sorted by decreasing estimated probability of being HLA allele associated. The column posterior_p denotes the posterior probability of the HLA regression coefficient being larger than 0, conditioned on the model and the observed data. A posterior probability of 1 denotes that the model is certain that a given replacement is an HLA-associated mutation, whereas a posterior probability of 0.5 means that there is no evidence of the replacement being HLA-associated (0.5 because a regression coefficient centered around 0 allocates equal probability to negative and positive values). ","category":"page"},{"location":"quick_start/","page":"Quick start","title":"Quick start","text":"The returned structure is a DataFrame object and elements can be accessed using standard indexing, e.g. summary[1:10, :posterior_p] returns the posterior probabilities of the first ten entries.","category":"page"},{"location":"quick_start/","page":"Quick start","title":"Quick start","text":"DataFrames can be saved as .csv files using the CSV package:","category":"page"},{"location":"quick_start/","page":"Quick start","title":"Quick start","text":"using CSV\n\nCSV.write(\"/home/user/Desktop/summary.csv\", summary)","category":"page"},{"location":"quick_start/","page":"Quick start","title":"Quick start","text":"The HLAData and run_model functions have an optional keyword argument save_file that can be used to save the converted input data and the result object.","category":"page"},{"location":"quick_start/","page":"Quick start","title":"Quick start","text":"data = HLAData(\n    alignment_file = \"/home/user/Desktop/alignment.fasta\",\n    tree_file = \"/home/user/Desktop/phylogeny.tree\",\n    save_file = \"/home/user/Desktop/input_data.jld2\"\n)\n\nresult = run_model(data, save_file = \"/home/user/Desktop/model_result.jld2\")\n\nsummary = replacement_summary(result)","category":"page"},{"location":"quick_start/","page":"Quick start","title":"Quick start","text":"The file extension \".jld2\" is used because HAMdetector uses the JLD2 package for saving and loading data.","category":"page"},{"location":"quick_start/","page":"Quick start","title":"Quick start","text":"To restore saved results from disk, use the load_data and load_result functions:","category":"page"},{"location":"quick_start/","page":"Quick start","title":"Quick start","text":"using Escape\n\nresult = load_result(\"/home/user/Desktop/model_result.jld2\")\ndata = load_data(\"/home/user/Desktop/input_data.jld2\")","category":"page"},{"location":"quick_start/#Providing-an-HLA-annotation-file","page":"Quick start","title":"Providing an HLA annotation file","text":"","category":"section"},{"location":"quick_start/","page":"Quick start","title":"Quick start","text":"For some datasets, sequences are not annotated with HLA information directly, but with a unique sequence identifier. For these cases, the HLAData function accepts an additional argument to provide a file that matches these sequence identifiers with the HLA annotation. The HLA annotation file has to be a comma-, tab- or semicolon-separated list of values with the columns \"identifier\", \"A1\", \"A2\", \"B1\", \"B2\", \"C1\", \"C2\". An example of a valid annotation file is shown below:","category":"page"},{"location":"quick_start/","page":"Quick start","title":"Quick start","text":"hla_annotations.csv","category":"page"},{"location":"quick_start/","page":"Quick start","title":"Quick start","text":"idenfier,A1,A2,B1,B2,C1,C2\nseq1,HLA-A*74,HLA-A*24,HLA-B*44,HLA-B*40,HLA-C*04,HLA-C*02\nseq2,HLA-A*11,HLA-A*11,HLA-B*51,HLA-B*07,HLA-C*01,HLA-C*01\nseq3,HLA-A*11,HLA-A*11,HLA-B*15,HLA-B*35,HLA-C*05,HLA-C*07\nseq4,HLA-A*24,HLA-A*23,HLA-B*54,HLA-B*41,HLA-C*06,HLA-C*05\nseq5,HLA-A*31,HLA-A*01,HLA-B*57,HLA-B*35,HLA-C*01,HLA-C*06","category":"page"},{"location":"quick_start/","page":"Quick start","title":"Quick start","text":"Note that, as with providing HLA annotation directly in the .fasta file, the order of the individual HLA alleles does not matter, e.g. a sequence annotated with HLA-B*18 in the B1 column and HLA-B*51 in the B2 column is treated the same as an otherwise identical sequence annotated with HLA-B*51 in the B1 column and HLA-B*18 in the B2 column.","category":"page"},{"location":"quick_start/","page":"Quick start","title":"Quick start","text":"HLA alleles do not have to be provided according to the standard HLA nomenclature, as HAMdetector is able to parse a broad range of commonly used HLA formats, also see section Parsing HLA section.","category":"page"},{"location":"quick_start/","page":"Quick start","title":"Quick start","text":"using Escape\n\ndata = HLAData(\n    alignment_file = \"/home/user/Desktop/alignment_with_ids.fasta\",\n    tree_file = \"/home/user/Desktop/phylogeny.tree\",\n    hla_annotation_file = \"/home/user/Desktop/hla_annotations.csv\"\n)","category":"page"},{"location":"quick_start/","page":"Quick start","title":"Quick start","text":"If an HLA annotation file is provided, the tips of the phylogenetic tree may either be labeled by the position of the corresponding sequence in the multiple sequence alignment - as explained previously- or by the unique sequence identifier.","category":"page"},{"location":"quick_start/#Convergence-diagnostics","page":"Quick start","title":"Convergence diagnostics","text":"","category":"section"},{"location":"quick_start/","page":"Quick start","title":"Quick start","text":"HAMdetector uses CmdStan for model fitting, which in turn uses a variant of Hamiltonian Monte Carlo sampling to generate samples from the typical set.  To check that these samples faithfully reflect the posterior distribution, HAMdetector uses CmdStan's diagnostic tool. The output can be retrieved with:","category":"page"},{"location":"quick_start/","page":"Quick start","title":"Quick start","text":"using Escape\n\nEscape.diagnostics(result)","category":"page"},{"location":"quick_start/","page":"Quick start","title":"Quick start","text":"The output is not nicely formated yet, but a successful model run should end in \"Processing complete, no problems detected\".","category":"page"},{"location":"quick_start/#Posterior-predictive-checks","page":"Quick start","title":"Posterior predictive checks","text":"","category":"section"},{"location":"quick_start/","page":"Quick start","title":"Quick start","text":"Bayesian statistics allows for a form of model checking called posterior predictive checks. A simple form of posterior predictive checks are calibration plots, which can be generated with:","category":"page"},{"location":"quick_start/","page":"Quick start","title":"Quick start","text":"using Escape, Plots\n\nEscape.calibration_plot(result)","category":"page"},{"location":"quick_start/","page":"Quick start","title":"Quick start","text":"For a well calibrated model, the vertical bars should span the diagonal line. Strong deviations from this suggest that something went wrong. A similar plot that specifically checks calibration of the internally calculated phylogeny can be generated with Escape.phylogeny_calibration(result).","category":"page"},{"location":"quick_start/","page":"Quick start","title":"Quick start","text":"Based on the quality of the provided phylogenetic tree, the calibration plot might show stronger deviations. This is usually not an issue because HAMdetector is able to handle noisy phylogenetic inputs.","category":"page"},{"location":"quick_start/#Parsing-HLA-alleles","page":"Quick start","title":"Parsing HLA alleles","text":"","category":"section"},{"location":"quick_start/","page":"Quick start","title":"Quick start","text":"For many datasets, HLA annotations do not follow the official HLA nomenclature. HAMdetector allows parsing of a wide range of HLA allele notations and tries to be as flexible as possible. To check if HAMdetector is able to parse HLA allele annotations in your dataset, use the parse_allele function:","category":"page"},{"location":"quick_start/","page":"Quick start","title":"Quick start","text":"using Escape\n\nEscape.parse_allele(\"B51\")\nEscape.parse_allele(\"B5101\")\nEscape.parse_allele(\"A*0103\")\nEscape.parse_allele(\"Cw*0302\")\nEscape.parse_allele(\"HLA-B51\")","category":"page"},{"location":"quick_start/","page":"Quick start","title":"Quick start","text":"Please note that HAMdetector right now supports either 2-digit or 4-digit HLA alleles, but no mixture of both. This limitation can in principle be solved by applying a technique called partial pooling,  which is going to be implemented in future release of HAMdetector.","category":"page"},{"location":"quick_start/","page":"Quick start","title":"Quick start","text":"Missing HLA alleles are supported by providing a blank HLA entry. Internally, missing HLA alleles are not imputed, but treated as 0 entries in the HLA matrix. This means that any effect that cannot be explained by one of the present HLA alleles is accounted for by the error term (the intercept), which is a reasonable default if HLA alleles are missing at random.","category":"page"},{"location":"#HAMdetector:-A-Bayesian-regression-model-that-integrates-information-to-detect-HLA-associated-mutations","page":"HAMdetector: A Bayesian regression model that integrates information to detect HLA-associated mutations","title":"HAMdetector: A Bayesian regression model that integrates information to detect HLA-associated mutations","text":"","category":"section"},{"location":"","page":"HAMdetector: A Bayesian regression model that integrates information to detect HLA-associated mutations","title":"HAMdetector: A Bayesian regression model that integrates information to detect HLA-associated mutations","text":"(Image: ) (Image: )  (Image: License: MIT)","category":"page"},{"location":"","page":"HAMdetector: A Bayesian regression model that integrates information to detect HLA-associated mutations","title":"HAMdetector: A Bayesian regression model that integrates information to detect HLA-associated mutations","text":"Escape.jl (HAMdetector) is a software package for identifying HLA associated substitutions based on annotated sequence data (aligned viral sequences paired to host HLA class I data).","category":"page"},{"location":"","page":"HAMdetector: A Bayesian regression model that integrates information to detect HLA-associated mutations","title":"HAMdetector: A Bayesian regression model that integrates information to detect HLA-associated mutations","text":"HAMdetector makes use of information from epitope prediction via MHCflurry 2.0 and phylogeny (based on RAxML-NG). The model is fit using Stan.","category":"page"},{"location":"","page":"HAMdetector: A Bayesian regression model that integrates information to detect HLA-associated mutations","title":"HAMdetector: A Bayesian regression model that integrates information to detect HLA-associated mutations","text":"See the documentation:","category":"page"},{"location":"","page":"HAMdetector: A Bayesian regression model that integrates information to detect HLA-associated mutations","title":"HAMdetector: A Bayesian regression model that integrates information to detect HLA-associated mutations","text":"(Image: )","category":"page"},{"location":"#Software-requirements","page":"HAMdetector: A Bayesian regression model that integrates information to detect HLA-associated mutations","title":"Software requirements","text":"","category":"section"},{"location":"","page":"HAMdetector: A Bayesian regression model that integrates information to detect HLA-associated mutations","title":"HAMdetector: A Bayesian regression model that integrates information to detect HLA-associated mutations","text":"Linux 64-bit\nJulia 1.6 or newer","category":"page"},{"location":"#Installation","page":"HAMdetector: A Bayesian regression model that integrates information to detect HLA-associated mutations","title":"Installation","text":"","category":"section"},{"location":"","page":"HAMdetector: A Bayesian regression model that integrates information to detect HLA-associated mutations","title":"HAMdetector: A Bayesian regression model that integrates information to detect HLA-associated mutations","text":"In the Julia REPL, switch to Pkg-mode by typing the ] character.  Then add HAMdetector_registry to your registry and install the Escape package:","category":"page"},{"location":"","page":"HAMdetector: A Bayesian regression model that integrates information to detect HLA-associated mutations","title":"HAMdetector: A Bayesian regression model that integrates information to detect HLA-associated mutations","text":"pkg> registry add https://github.com/HAMdetector/HAMdetector_registry.git\npkg> add Escape","category":"page"},{"location":"","page":"HAMdetector: A Bayesian regression model that integrates information to detect HLA-associated mutations","title":"HAMdetector: A Bayesian regression model that integrates information to detect HLA-associated mutations","text":"To make sure that HAMdetector was installed correctly, a large suite of unit tests can be run with ] test Escape.","category":"page"},{"location":"#Reporting-Issues-and-Contributions","page":"HAMdetector: A Bayesian regression model that integrates information to detect HLA-associated mutations","title":"Reporting Issues and Contributions","text":"","category":"section"},{"location":"","page":"HAMdetector: A Bayesian regression model that integrates information to detect HLA-associated mutations","title":"HAMdetector: A Bayesian regression model that integrates information to detect HLA-associated mutations","text":"If you experience issues during installation or usage HAMdetector, please open an issue in this repository. Contributions in the form of suggestions, improvements to the documentation or pull requests are extremely welcome.","category":"page"},{"location":"#License","page":"HAMdetector: A Bayesian regression model that integrates information to detect HLA-associated mutations","title":"License","text":"","category":"section"},{"location":"","page":"HAMdetector: A Bayesian regression model that integrates information to detect HLA-associated mutations","title":"HAMdetector: A Bayesian regression model that integrates information to detect HLA-associated mutations","text":"The HAMdetector source code is distributed under a permissive MIT license (see file LICENSE).","category":"page"}]
}
