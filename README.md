# fplib 

### Purpose and Advantages
fplib is a data plotting and curve-fitting library. It

- automates the redundant and time-consuming raw data plots, curve-fitting, residual plots and Latex-figure-Creation tasks
into ~10 lines which can simply be copy-pasted every time
- prints full latex code for curve-fitting figure + results for direct copy-paste into Overleaf
- is specifically tailored to comply with the formatting regulations of the "Grundpraktikum (GP)", 
and by extension, the "Fortgeschrittenenpraktikum (FP)"
- has lots of pre-defined models, but also enables easy custom model definitions

From the FP lectures, it copies exactly the color scheme, font size, figure format and size, labels, and 
lists optimized parameter results and $\chi^2 /N_{dof}$ in a structured way. 

As a result, the figures created with this library leave no room for point deduction due to formatting errors.

### Content
fplib is structured into two modules: 

| | |
|---|---|
| `fpplot` | plotting (e.g., raw data) |
| `fpfit`  | curve fitting |

### Installation
Simply download the repository by clicking on Code -> Download Zip 

Or clone the repository via 
```bash 
git clone git@github.com:aristidgrossmann/FPlib.git
```

fplib is designed to be easily customized and extended, thats why no build instructions are described here.


### Repository Structure

<pre>
FPlib/
├── fplib/                       # contains the library
│   ├── fpplot.py                # plotting functions
│   ├── fpfit.py                 # curve fitting functions
│   ├── ModelTemplate.py         # abstract model class (inherit from this when defining custom models)
│   ├── utils/                   # utility functions used by fpplot and fpfit
│   └── models/                  # lots of pre-defined models here 
|
├── examples/                    # examples 
│   └── plots/                   # plotting examples
│   └── fits/                    # curve fit examples
|
├── originalVersion/
│   └── fp_library.py            #old and messy version (dont look at it)
|
├── .gitignore
├── README.md
└── requirements.txt

</pre>
---

## Plotting

In your script, import `fpplot` via

```python 
from fplib import fpplot
```

From `fpplot` you have access to 4 methods: 
| Method | Description |
|------|------------|
| `plot_raw_data_no_uncertainty` | plots x-y-data without errorbars |
| `plot_raw_data_no_uncertainty_zoom_window` | same as above. On top of that, it plots a zoom window on a specified region |
| `plot_raw_data_with_uncertainty` | plots x-y-data with x-y-errorbars |
| `plot_raw_data_with_uncertainty_zoom_window` | same as above. On top of that, it plots a zoom window on a specified region |


The 4 methods have the following **required arguments**: 


| Argument | Required by | Type |  Description |
|------|------------|----------------|-----------|
| `xdata` | all | `np.ndarray` | x data |
| `ydata` | all | `np.ndarray` | y data |
| `xerr` | `plot_raw_data_with_uncertainty`, `plot_raw_data_with_uncertainty_zoom_window` | `np.ndarray` | uncertainty on xdata |
| `yerr` | `plot_raw_data_with_uncertainty`,  `plot_raw_data_with_uncertainty_zoom_window` | `np.ndarray` | uncertainty on ydata |
| `title` | all | `str` | plot title |
| `xlabel` | all | `str` | x axis label |
| `ylabel` | all | `str` | y axis label |
| `legend_loc` | all | `str` | location of the legend, f.e.  `upper center`  (see matplotlib for options) |
| `file_name` | all | `str` | name under which the plot is saved |
| `zoom_xlim` | `plot_raw_data_no_uncertainty_zoom_window`, `plot_raw_data_with_uncertainty_zoom_window` | `str` | x range of the zoom region |
| `zoom_ylim` | `plot_raw_data_no_uncertainty_zoom_window`, `plot_raw_data_with_uncertainty_zoom_window` | `str` | y range of the zoom region |
| `scaling_factor` | `plot_raw_data_no_uncertainty_zoom_window`, `plot_raw_data_with_uncertainty_zoom_window` | `float` | scales the size of the zoom window |
| `zoom_window_position` | `plot_raw_data_no_uncertainty_zoom_window`, `plot_raw_data_with_uncertainty_zoom_window` | `str` | location of the zoom window, f.e.  `upper center`  (see matplotlib for options) |


Additionally, all 4 methods have the following **optional arguments**: 


| Argument | Optional for | Type |  Description |
|------|------------|----------------|-----------|
| `xlims` | all | `tuple` | x axis limits. Default: `None`. Can be set to for example (2, 10) |
| `ylims` | all | `tuple` | y axis limits. Default: `None`. Can be set to for example  (2, 10) |
| `custom_plot` | all | `boolean` | Default: `False`. Setting to `True` returns the axes object, allowing further modifications |


---


## Curve Fitting


### Methodology

In analogy to the `Praktikumsbibliothek` provided for the GP (see https://grundpraktikum.physik.rwth-aachen.de/software/), curve fits of a function $f(x|\boldsymbol{a})$ with parameters $\boldsymbol{a}$ on data with uncertainty only on the y-value $\sigma_y$ by minimizing 
$\chi^2 = \sum_i \frac{(y_i - f(x_i|\boldsymbol{a}))^2}{\sigma_{y,i}^2}$.
The function $\mathrm{scipy.optimize.curve \\_fit}$ is used to solve the minimization problem. 

When presented with uncertainty in the x and y values simultaneously, the Orthogonal Distance Regression (ODR) method  is used, which minimizes a similar objective function. To solve the minimization problem numerically, the algorithm $\mathrm{scipy.odr}$ in analogy to the `Praktikumsbibliothek` is utilized.

In accordance with the GP protocol guidelines, the fit parameters and their standard deviations are reported. After consulting the FP advisors, correlations are not reported.

To quantify the quality of our fit, the value for normalized
$\chi^2/N_{\mathrm{dof}} = \frac{1}{N_{\mathrm{dof}}}\sum_i \frac{(y_i - f(x_i|\textbf{a}))^2}{\sigma_{y,i}^2 + (f'(x_i|\textbf{a})\sigma_{x,i})^2}$ is reported.


Additionally, the residuals $y_i - f(x_i|\textbf{a})$ with the uncertainty $\sqrt{\sigma_{y,i}^2 + (f'(x_i|\textbf{a})\sigma_{x,i})^2}$ are plotted, in analogy to the GP guidelines.


### fpfit

In your script, import `fpfit` via

```python 
from fplib import fpfit
```
Inside `fpfit`, the function  `general_curve_fit` supports all-purpose curve fitting + redidual plots + LaTex figure code generation. 
The code for LaTex figures is written to a .txt file with the same name as `file_name`.

Its returns are:
| Return variable | Type |  Description |
|------|------------|----------------|
| `popt` |  `np.ndarray` | optimized model parameters in the order that they were specified in |
| `popt_std` |  `np.ndarray` | uncertainties on optimized model parameters in the same order |
| `pcorr` |  `np.ndarray` | correlation matrix |

They can be captured via 
```python 
popt, popt_std, pcorr = fpfit.general_curve_fit(....)
```

`general_curve_fit` has the following **required and optional arguments**: 

| Argument | Required | Type |  Description |
|------|------------|----------------|-----------|
| `xdata` | yes | `np.ndarray` | x data |
| `ydata` | yes | `np.ndarray` | y data |
| `xerr` | yes | `np.ndarray` | uncertainty on xdata |
| `yerr` | yes | `np.ndarray` | uncertainty on ydata |
| `model` | yes | `Model` | model class that is fit to the data (see further down) |
| `p0` | yes | `np.ndarray` | initial guess. Must have the same length as the number of model parameters |
| `fit1_label` | yes | `str` | label of the fit line (in the legend) |
| `plot_uncertainties` | yes | `boolean` | `True`: plots uncertainties. `False`: only data points without error bars    |
| `xlabels` | yes | `str` | x axis label |
| `curveFitPlot_ylabel` | yes | `str` | y axis label of the fit plot |
| `residualPlot_ylabel` | yes | `str` | y axis label of the Residual plot |
| `plot_title` | yes | `str` | Title of the plot |
| `legend_loc` | yes | `str` | location of the legend in the curve-fit-plot, f.e.  `upper center`  (see matplotlib for options) |
| `file_name` | yes | `str` | file name under which the plot and the Latex figure code is saved |
| `xlims` | optional | `tuple` | x axis limits. Default: `None`. Can be set to for example (2, 10) |
| `ylims` | optional | `tiple` | y axis limits of the top plot. Default: `None`. Can be set to for example (2, 10) |
| `compressed_Latex_output` | optional | `bool` | Default: `True`. Setting to `False` prints more curve fit info (f.e. correlation) to the Latex txt output  |
| `custom_plot` | optional | `boolean` | Default: `False`. Setting to `True` returns the axes objects of the curve fit and the residual plots, allowing further modifications |
| `exclude_indices` | yes | `np.ndarray` | Indices of the data points that are excluded when fitting. Data points are still plotted though |
| `exclude_zero_count_data_points` | optional | 'boolean | Default: `False`. Set to `True` when working with count data (since 0-counts have zero uncertainty)  |
| `peak_index` | optional | `int` | Only used when fitting `DoubleGaussian`. Can be set to `0` or `1` to specify which peak is the background|
| `peak1_label` | optional | `str` | Only used when fitting `DoubleGaussian`. Label of the first peak|
| `peak2_label` | optional | `str` | Only used when fitting `DoubleGaussian`. Label of the first peak|


---


## Models

There already are a lot of pre-defined models in the `models` folder. However, custom models can also be defined easily. 
Here an example, which defines the custom model $f(x) = A \cos^2(\omega x + \phi) + a_0 + a_1 x + a_2 x^2$ with the model parameters 
$A, \omega, \phi, a_0, a_1, a_2$

```python 
from fplib.ModelTemplate import ModelTemplate
import numpy as np

class CustomModel(ModelTemplate)  #inherit from ModelTemplate

    # Model description (must me specified)
    MODEL_NAME = 'Custom Model'   
    MODEL_EQUATION_LATEX = "f(x) = A \cos^2(\omega x + \phi) + a_0 + a_1 x + a_2 x^2"
    MODEL_PARAMETER_LABELS = ["$A$", "$\\omega$", "$\\phi$", "$a_0$", "$a_1$", "$a_2$"]

    @staticmethod
    def modelFunction(data, A, omega, phi, a_0, a_1, a_2):   #must be in the order data, *params
        return A * np.cos(omega *data + phi)**2 + a_0 + a_1 *x + a_2 * x**2
    
    @staticmethod
    def modelFunctionDerivative(data, C, mu, sigma):   #must be in the order data, *params
        return - 2*A*omega* np.sin(omega *data + phi)*np.cos(omega *data + phi)  + a_1 + 2*a_2
```
After that, you can simply import your model and pass it to the all-purpose curve fit function `fpfit.general_curve_fit`. For example, if your model is inside a file named `CustomModel.py` inside the `models` folder, simply import it via 
```python 
from fplib.models.CustomModel import CustomModel
```

The LaTex figure code generation is automated for all models inheriting from `ModelTemplate`, but if you wish to implement a custom Result formatting, simply overwrite the `@classmethod` `print_curve_fit_popt_general` (for example as in `DoubleGaussian`). 

Also, if you define a model whose number of input parameters is not fixed, you have to overwrite the `@classmethod` `getModelParameterLabels` (for example as in `NGaussians`). 








