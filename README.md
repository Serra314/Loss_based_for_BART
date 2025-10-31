# Loss-based prior for CART and BART
The repository contains the code used to produce the results and the figure illustrated in the article "Loss-based prior for CART and BART models"

The code folder contains the functions for the MCMC algorithm, analyse posterior results, visualise trees for the CART model.
The code for the BART model implementation is in the `diabetes_application` folder

The file `Prior_figures.R` produces Figures 1 to 5

The file `Synthetic_data_example.R` contains the code to compute the results of the synthetic experiment in Section 5. The code produces Figures 6 to 8, the results in Table 1, and Figures from 1 to 6 of the supplementary material

The file `Breast_cancer_example.R` contains the code to compute the result of the real data example on breast cancer data reported in Section 6.1. The code produces Figures 9 to 13, the results in Table 2, and Figures from 7 to 10 of the supplementary material

The folder `diabetes_application` contains data and code to reproduce the result of the real data example on diabetes data reported in Section 6.2. The script `code_to_run.R` contains code to produce Figures 14, and 15, and Figures 11, and 12 of the supplementary material.
