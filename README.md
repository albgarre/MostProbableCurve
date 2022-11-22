# Code for The Most Probable Curve method - A robust approach to estimate kinetic models from low plate count data resulting in reduced uncertainty

This repository includes the replication code for the article The Most Probable Curve method - A robust approach to estimate kinetic models from low plate count data resulting in reduced uncertainty (Garre et al., 2022; International Journal of Food Microbiology, 380, 109871). The article is available in Open Access here: 

https://doi.org/10.1016/j.ijfoodmicro.2022.109871

# Abstract 

A novel method is proposed for fitting microbial inactivation models to data on liquid media: the Most Probable Curve (MPC) method. It is a multilevel model that makes a separation between the “true” microbial concentration according to the model, the “actual” concentration in the media considering chance, and the actual counts on the plate. It is based on the assumptions that stress resistance is homogeneous within a microbial population, and that there is no aggregation of microbial cells. Under these assumptions, the number of colonies in/on a plate follows a Poisson distribution with expected value depending on the proposed kinetic model, the number of dilutions and the plated volume.

The novel method is compared against (non)linear regression based on a normal likelihood distribution (traditional method), Poisson regression and gamma-Poisson regression using data on the inactivation of Listeria monocytogenes. The conclusion is that the traditional method has limitations when the data includes plates with low (or zero) cell counts, which can be mitigated using more complex (discrete) likelihoods. However, Poisson regression uses an unrealistic likelihood function, making it unsuitable for survivor curves with several log-reductions. Gamma-Poisson regression uses a more realistic likelihood function, even though it is based mostly on empirical hypotheses. We conclude that the MPC method can be used reliably, especially when the data includes plates with low or zero counts. Furthermore, it generates a more realistic description of uncertainty, integrating the contribution of the plating error and reducing the uncertainty of the primary model parameters. Consequently, although it increases modelling complexity, the MPC method can be of great interest in predictive microbiology, especially in studies focused on variability analysis.



