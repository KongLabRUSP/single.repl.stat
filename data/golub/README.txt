Source: 
https://www.kaggle.com/crawford/gene-expression/data

Description from Kaggle:
https://www.kaggle.com/crawford/gene-expression

Context
This dataset comes from a proof-of-concept study published in 1999 by Golub et al. It showed how new cases of cancer could be classified by gene expression monitoring (via DNA microarray) and thereby provided a general approach for identifying new cancer classes and assigning tumors to known classes. These data were used to classify patients with acute myeloid leukemia (AML) and acute lymphoblastic leukemia (ALL).

Content
Golub et al "Molecular Classification of Cancer: Class Discovery and Class Prediction by Gene Expression Monitoring"

There are two datasets containing the initial (training, 38 samples) and independent (test, 34 samples) datasets used in the paper. These datasets contain measurements corresponding to ALL and AML samples from Bone Marrow and Peripheral Blood. Intensity values have been re-scaled such that overall intensities for each chip are equivalent.

Acknowledgements
Molecular Classification of Cancer: Class Discovery and Class Prediction by Gene Expression

Science 286:531-537. (1999). Published: 1999.10.14

T.R. Golub, D.K. Slonim, P. Tamayo, C. Huard, M. Gaasenbeek, J.P. Mesirov, H. Coller, M. Loh, J.R. Downing, M.A. Caligiuri, C.D. Bloomfield, and E.S. Lander

These datasets have been converted to a comma separated value files (CSV).

Inspiration
These datasets are great for classification problems. The original authors used the data to classify the type of cancer in each patient by their gene expressions.

Description from Dhammika & Javier's book (Chapter 1.4.1 The Golub Data):
The Golub data first appeared in ...
27 patients with acute lympoblastic leukemia (ALL) 
11 with acute myeloid leukemia (AML)
3051 genes measured

NOTE: Istead of Kaggle data, used R dataset from multtest package: data(golub)