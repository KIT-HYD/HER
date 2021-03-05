# HER: an information theoretic alternative for geostatistics

In HER method, we propose a stochastic, geostatistical estimator which combines information theory with probability aggregation methods for minimizing predictive uncertainty, and predicting distributions directly based on empirical probability. Histogram via entropy reduction (HER) relaxes parametrizations, avoiding the risk of adding information not present in data (or losing available information). It provides a framework for uncertainty estimation that takes into account both spatial configuration and data values, while allowing to infer (or introduce) continuous or discontinuous characteristics of the field. 
We investigate the framework utility using synthetically generated datasets from Gaussian Processes with different sample sizes and data properties (different spatial correlation distances and addition of noise). 
HER method brings a new perspective of spatial interpolation and uncertainty analysis to geostatistics and statistical learning, using the lens of information theory.

The code and datasets are complementary parts of the study proposed by Thiesen, Vieira, Mälicke, Loritz, Wellmann and Ehret (2020):

>_Thiesen, S.; Vieira, D.; Mälicke, M.; Loritz, R.; Wellmann, J. F.; Ehret, U. Histogram via entropy reduction (HER): an information-theoretic alternative for geostatistics, Hydrol. Earth Syst. Sci., https://doi.org/10.5194/hess-24-4523-2020, 24(9), 4523–4540, 2020._ 

## License agreement

The HER method comes with ABSOLUTELY NO WARRANTY. You are welcome to modify and redistribute it within the license agreement. The HER method is published under the CreativeCommons "CC-BY-4.0" license together with a ready-to-use sample data set. To view a full version of the license agreement please visit [CC-BY-4.0](https://creativecommons.org/licenses/by/4.0/).

## Requisites

* MATLAB (tested on 2018b).

## Usage

See HER.m

## File structure

* HER script ... .m
* functions/ ... .m
* datasets/ ... .mat

### HER

The script is divided in 7 sections:

__1. Load dataset__
	Loads the dataset.
	
__2. Define infogram and Geo3 properties__
	Definition of the infogram properties, aggregation method, z threshold (optional).
	
__3. Geo1: Spatial characterization__
	Extracts spatial correlation patterns. ```f_her_infogram.m```

__4. Geo2: Weight optimization__
	Optimizes weights for the aggregation method based on entropy minimization. ```f_her_weight.m```

__5. Geo3: z PMF prediction__
	Applies spatial characterization and optimal weights for PMF prediction. ```f_her_predict.m```
	
__6. Extract PMF statistics__
	Obtains mean, median, mode and probability of a z threshold (optional) of the predicted z PMFs and plots the results.

__7. Calculate performance metrics__
	Calculates Root Mean Square Error (RMSE), Mean Error (ME), Mean Absolute Error (MAE), Nash-Sutcliffe model efficiency and scoring rule (DKL) of the validation set.  
	
__8. Clear__
	Clears intermediate variables.

### Functions

The functions are detailed in their own source code body. Examples of how to use them are available in the `HER.m` script. 

```
f_DKL_w_AND.m
f_DKL_w_OR.m
f_diff.m
f_entropy.m
f_euclidean_dist.m
f_linear_aggregation.m
f_loglinear_aggregation.m
f_performance_det.m
f_performance_prob.m
f_her_infogram.m
f_her_weight.m
f_her_predict.m
f_extract_pmf_statistics.m
f_plot_infogram.m
f_plot_weights.m
f_plot_prediction.m
f_plot_probabilitymap.m
```

```f_plot``` functions were specifically built for the dataset of the study.

### Dataset of the study

The folder contains synthetic observations used in the paper case study. Four synthetic 2D spatial datasets with grid size 100x100 were generated from known Gaussian processes. We use rational quadratic kernel as the covariance function, with correlation lengths of 6 and 18 units. For both, short- and long-range fields, a white noise was introduced given by Gaussian distribution with mean 0 and standard deviation equal to 0.5.

The generated sets comprise:
	* SR0: short-range field without noise 
	* SR1: short-range field with noise
	* LR0: long-range field without noise 
	* LR1: long-range field with noise
We randomly shuffled the data, and then divided it in three mutually exclusive sets: one to generate the calibration subsets (sizes of 200, 400, 600, 800, 1000, 1500, and 2000), one for validation (containing 2000 data points), and another 2000 data points as test set.

Each dataset file contains:

* __idx_rand_full:__ index of the randomly shuffled data (same for all files)
* __sample_size:__ all calibration sizes available of the dataset (same for all files)
* __data:__ matrix with z values of the full generated dataset
* __txt:__ dataset type (SR0, SR1, LR0, LR1)
* __idx_cal:__ index of the calibration set
* __idx_val:__ index of the validation set
* __idx_test:__ index of the test set
* __x:__ matrix with x coordinates of the full dataset
* __x_cal:__ vector with x coordinates of the calibration set (x_cal=x(idx_cal))
* __x_val:__ vector with x coordinates of the validation set (x_val=x(idx_val))
* __x_test:__ vector with x coordinates of the test set (x_test=x(idx_test))
* __y:__ matrix with y coordinates of the full dataset
* __y_cal:__ vector with y coordinates of the calibration set (y_cal=y(idx_cal))
* __y_val:__ vector with y coordinates of the validation set (y_val=y(idx_val))
* __y_test:__ vector with y coordinates of the test set (y_test=y(idx_test))
* __z:__ matrix with z values of the full generated dataset (z=data)
* __z_cal:__ vector with z values of the calibration dataset (z_cal=z(idx_cal))
* __z_val:__ vector with z values of the validation dataset (z_val=z(idx_val))
* __z_test:__ vector with z values of the test dataset (z_test=z(idx_test))
* __dim_cal:__ size of the calibration set (dim_cal=length(idx_cal))
* __dim_val:__ size of the validation set (dim_val=length(idx_val))
* __dim_test:__ size of the test set (dim_test=length(idx_test))

The synthetic field generator, using Gaussian processes, is available in scikit-learn (Pedregosa et al., 2011), while the code producing the fields can be found at
https://github.com/mmaelicke/random_fields.

## Contact

Stephanie Thiesen | stephanie.thiesen@kit.edu
Uwe Ehret | uwe.ehret@kit.edu

