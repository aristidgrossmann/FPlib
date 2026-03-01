# fplib 

### Purpose and Advantages
fplib is a data plotting and curve-fitting library. It

- automates the redundant and time-consuming plotting, curve-fitting and Latex-figure-Creation tasks
- prints full latex code for curve-fitting figure + results for direct copy-paste into Overleaf
- specifically tailored to comply with the formatting regulations of the "Grundpraktikum (GP)", 
and by extension, the "Fortgeschrittenenpraktikum (FP)"
- lots of pre-defined models, but also enables easy custom model definitions

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


| Argument | Required by | type |  Description |
|------|------------|----------------|-----------|
| `xdata` | all | `np.ndarray` | x data |
| `ydata` | all | `np.ndarray` | y data |
| `xerr` | `plot_raw_data_with_uncertainty`, `plot_raw_data_with_uncertainty_zoom_window` | `np.ndarray` | x errors |
| `yerr` | `plot_raw_data_with_uncertainty`,  `plot_raw_data_with_uncertainty_zoom_window` | `np.ndarray` | y errors |
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


| Argument | Optional for | type |  Description |
|------|------------|----------------|-----------|
| `xlims` | all | `tuple` | x axis limits. Default: `None`. set to f.e. (2, 10) |
| `ylims` | all | `tuple` | y axis limits. Default: `None`. set to f.e. (2, 10) |
| `custom_plot` | all | `boolean` | Default: `False`. Setting to `True` returns the axes object, allowing further modifications |


