---
title: 'EPP-TB: The ERP Post-Processing Tool Box  <img src="doc\logo.png" align="right" />'
output:
  html_document:
    keep_md: true
---



Once completing pre-processing in eeglab/erplab/Net Station, you're ready to get to the fun stuff: plotting and measuring with EPP-TB!
This package is aimed at simple (readable), concise (one function per action) and reproducible code writing.
Additionally, a major component is the ability to export data and plots to be used and further manipulated elsewhere (say, `ggplot2`?).

All user-end functions start with `epp_*`.

## Getting Started

### Prerequisites

To use EPP-TB you will need:

- Matlab (2015a+).
- [eeglab](https://sccn.ucsd.edu/eeglab/index.php) (14.X.X) for importing `.set` files, and plotting topos.
- R + [`tidyverse`](https://github.com/tidyverse/tidyverse) + [`scales`](https://github.com/hadley/scales) for plotting in R.

### Installing

Install the package [downloading](https://github.com/mattansb/EPP-TB/releases) and adding the EPP-TB folder (and sub-folders) to your Matlab paths.

## Importing Data into an EPP structure

Three import methods are currently supported:

- `epp_loadeeglab` - Import multiple `.set` files from eeglab (supports wavelet analysis based on [Mike X. Cohen's great book and code](http://www.mikexcohen.com/)).
- `epp_loaderplab` - Import from [erplab](http://www.erpinfo.org/erplab.html).
- `epp_loadegimat` - Import multiple `.mat` files exported from Net Station (EGI).

The resulting structure array has `length(struct)` equal to the number of conditions, and contains the following fields, per condition:

- **Condition**: the name of the conditions.
- **timeLine**: a vector of time points.
- **IDs**: a table with two variables - _**ID**_ and _**nTrials**_ (the number of trials a ERP has been averaged across).
- **Data**: a channels$\times$`length(timeLine)`$\times$`length(IDs)` matrix for ERP data.

If a wavelet analysis has been preformed, the Data field is replaced with:

- **Freqs**: a vector of frequencies used.
- **ersp** and **itc**: channels$\times$`length(Freqs)`$\times$`length(timeLine)`$\times$`length(IDs)` matrices.

### Data Reduction

These functions can be used to compute new conditions or manipulate existing data:

- Merge 2 or more conditions with `epp_combineconds`.
- Compute differences between 2 conditions with `epp_diffwave`.
- Compute LRP with `epp_LRP`
- Compute global field potentials with `epp_GFP`.



## Plotting

### ERP Plots

#### Grand averages
Grand averages can be plotted by specifying the conditions and channel indices (averaged across) to plot:

```Matlab
conds         = {'Control_Crr_plus','Control_L5_plus'};
channel_inds  = [5 6 11 12];
epp_plotgrands(study,conds,channel_inds)
```

![](doc/grand_plot-1.png)<!-- -->


You can also plot error envelopes:

```Matlab
epp_plotgrands(study,conds,channel_inds,'errorType','SE')
% can also be set to 'SD' or 'CIXX' 
% (with XX replaced with any percent: 'CI95','CI80', etc...).
```

![](doc/grand_plot_se-1.png)<!-- -->


#### Topo Plots

Topo plotting is dependent on eeglab functions. Additionally, you will need to provide a eeglab-like chanlocs structure to these functions, with `length(chanlocs)==size(study.Data,1)`.

```Matlab
times = [185];
epp_plottopo(study,chanlocs,conds,channel_inds,times)
```
![](doc/topo_plot-1.png)<!-- -->

#### Butterfly and Trace Plots

Butterfly plots can be used to plot the mean ERP for each subject individually.

```Matlab
epp_plotbutterfly(study,conds,channel_inds)
```

![](doc/butterfly_plot-1.png)<!-- -->


Trace plots are similar to butterfly plots, but the mean activation (across subjects) is plotted for each channel separately.

```Matlab
epp_plotbutterfly(study,conds,channel_inds,'all',true)
% if channels are indicated, they are ignored.
```
![](doc/trace_plot-1.png)<!-- -->

### TF Plots

Time Frequency plots plot both ersp and itc:

```Matlab
epp_plotTF(study,conds,channel_inds)
```

![](doc/TFplot-1.png)<!-- -->

Topo Plots-  is a wrapper around `epp_plottopo`.
Need to specify bands:

```Matlab
bands = [4, 8; 8, 12];
% a matrix of frequencies, with each row containing a range of frequencies
% to plot (1st column is lower limit, 2nd column is upper limit of band).
epp_plottopoTF(study,chanlocs,conds,channel_inds,times,bands)
```
![](doc/TFtopo-1.png)<!-- -->

### Exporting to R

All plots can be exported to R and plotted with [`ggplot2`](https://github.com/tidyverse/ggplot2) by setting `'R',true` in any of the plotting function. This produces two time-stamped files: 

- A data file (`*_data.csv`)
- A code file (`*_code.R`), to plot said data using `ggplot2`.


## Measuring

When saving, all are saved in a 2-sheet xlsx file, with the second sheet containing the parameters used in measuring - this is to insure reproducibility of results. Thus, if you have the output `.xlsx` file, you'll never find yourself asking "what was the time window we used last time?".

### ERP Amplitude

The function `epp_getamplitude` implements the methods for measuring amplitudes described in chapter 9 of [Steven J. Luck's intro to ERP book](https://mitpress.mit.edu/books/introduction-event-related-potential-technique-0).

 - (Local) Peak amplitude
 - Mean amplitude
 - Integral
 - Area

### ERP Latency

The function `epp_getlatency` implements the methods for measuring latencies described in chapter 9 of [Steven J. Luck's intro to ERP book](https://mitpress.mit.edu/books/introduction-event-related-potential-technique-0), as well as those described in [Kiesel et al.'s Jakknife paper](http://onlinelibrary.wiley.com/doi/10.1111/j.1469-8986.2007.00618.x/abstract).

- (Local) Peak latency
- Relative criterion
- Fractional area
- Baseline deviation
- Absolute criterion

### TF Power

The function `epp_getTF` measure the *mean* ersp / itc from selected channels, within a specified time-window,separately for specified frequency bands (as shown above, for `epp_plottopoTF`).

## Authors

- **Mattan S. Ben-Shachar** [aut, cre].
- **Rachel Rac** [ctb].

## Acknowledgments

- **Mike X. Cohen** - who's [code](http://www.mikexcohen.com/) was implemented in `suppWaveletConv3.m`.
- **Matt Craddock** - who's [code](https://craddm.github.io/blog/2017/02/25/EEG-topography) was implemented in `epp_plottopo.R`.
- **Winston Chang** - who's [code](http://www.cookbook-r.com/Graphs/Plotting_means_and_error_bars_(ggplot2)/) for computing within-subject (and mixed model) CI's was implemented in `epp_plotgrands.m`.
