""""
This python script contains the find_hvgs() function, which identifies highly-variable genes according to the method of Klein et al.
It's been adapted from the code of Sam Wollock (helper_funtions.py, get_vscores() and filter_genes(), to run on modern scanpy `adata` objects
(Original github link:  https://github.com/AllonKleinLab/klunctions/blob/77d2571e5083faed53b3950ffffc9dd07f60a9e5/sam/Analysis/scBasics/helper_functions.py)

Usage: find_hvgs(adata)
Output: indices/locations of genes in adata.var which have been identified as highly variable genes
    if plot_result = True: a plot of mean expression versus fano factor, indicating which genes have been delineated as highly variable
"""
import numpy as np
import pandas as pd
import scipy.optimize
import matplotlib.pyplot as plt


def find_hvgs(adata, min_mean = 0, min_cells = 3, min_counts = 3, min_vscore_percentile = 85, plot_result = True):
    """
    Identifies highly-variable genes according to the method of Klein et al. 2015
    Adapted by Laura Bagamery from code written by Sam Wollock (helper_fuctions.py, get_vscores() and filter_genes())
    github: https://github.com/AllonKleinLab/klunctions/blob/77d2571e5083faed53b3950ffffc9dd07f60a9e5/sam/Analysis/scBasics/helper_functions.py

    This function assumes normalized counts but NOT log-transformed data

    Parameters:
    min_mean: HVGs must have, at minimum, a mean expression level of `min_mean`
    min_cells: HVGs must be expressed in, at minimum, `min_cells` number of cells
    min_vscore_percentile: minimum v-score (above-poisson noise) for calling HVGs
    """

    #E = adata.copy()
    E = adata.X.copy()
    E = E.todense()
    E = np.asarray(E)
    # find mean expression for each gene
    mu_genes = np.mean(E, axis = 0)
    # find indices for genes expressed above minimum threshold
    i_genes = np.where(mu_genes > min_mean)[0]
    # select only genes above minimum threshold
    mu_genes = mu_genes[i_genes]
    # calculate variance
    var_genes = np.var(E[:, i_genes], axis = 0)

    # calculate fano factor
    fano = var_genes / mu_genes

    x = np.log(mu_genes)
    y = np.log(fano / mu_genes)

    # bin into `n_bins` based on values of x
    # find percentile `p_bin` of y within each [x-defined] bin
    n_bins = 50
    p_bin = 0.1
    x_by_bin, bins = pd.cut(x, n_bins, labels = False, retbins = True)
    # convert values from edges of bins (n_bins + 1 values)
    # to centers of binds (n_bins values) by moving average
    x_bins = np.array(pd.Series(bins).rolling(window = 2).mean())[1:]
    # iterate over bins
    y_bins = np.zeros(x_bins.shape)
    for i_bin in range(n_bins):
        # find y-values, quantile within bin
        subset = np.where(x_by_bin == i_bin)[0]
        if len(subset) > 0:
            y_bins[i_bin] = np.percentile(y[subset], p_bin)
        # if nothing in bin: replace with previous value
        # or NaN, if necessary (the first value)
        else:
            if i_bin == 0:
                y_bins[i_bin] = np.nan
            else:
                y_bins[i_bin] = y_bins[i_bin - 1]

    # drop NaNs (if applicable) generated in y; slice x accordingly
    y_reals = np.where(y_bins == y_bins)[0]
    x = x_bins[y_reals]
    y = y_bins[y_reals]

    def glog(x_, b_, c_):
        return np.log((c_*np.exp(-x_)) + b_)

    # histogram of fano factor
    hist, hbins = np.histogram(np.log(fano), bins = 200)
    # set bin values to the center of the bin ranges
    hbins = np.array(pd.Series(hbins).rolling(window = 2).mean()[1:])
    i_max = np.argmax(hist)
    # calculate c from bin with highest frequency of genes (or default to 1, if < 1)
    c = np.max((np.exp(hbins[i_max]), 1))

    def errorfunc(b, x_, c_, y_):
        return np.sum(abs(glog(x_, b, c_) - y_))

    # minimize error function
    # calculate parameters as in Klein dropseq paper
    b = scipy.optimize.fmin(func = errorfunc, x0 = [0.1], args = (x, c, y,), disp = False)
    a = (c / (1 + b)) - 1

    # calculate v scores, cvs (Klein 2015 theory supplement)
    v_scores = fano / ((1 + a) * (1 + b) + b * mu_genes)
    cv_eff = np.sqrt((1 + a) * (1 + b) -1)
    cv_input = np.sqrt(b)

    # find v scores > 0; slice other statistics accordingly
    iv = np.where(v_scores > 0)[0]
    i_genes = i_genes[iv]
    mu_genes = mu_genes[iv]
    fano = fano[iv]

    # threshold based on determined percentile rank of v score
    min_vscore = np.percentile(v_scores, min_vscore_percentile)
    # slicing: expression data for i_genes (above min_mean), above min_count
    # which are expressed in at least min_cells
    # and which have a v-score above the minimum threshold
    i_hvgs = (((E[:, i_genes] >= min_counts).sum(axis = 0).squeeze()  >= min_cells) & (v_scores >= min_vscore))

    # plotting the results
    if plot_result == True:
        # set x range to plot--0.5x below and 2x above the observed extrema
        x_min = 0.5 * np.min(mu_genes)
        x_max = 2 * np.max(mu_genes)
        # x-values, x-values along threshold line to plot
        x_threshold = x_min * np.exp(np.log(x_max / x_min) * np.linspace(0, 1, 100))
        y_threshold = (1 + a) * (1 + b) + (b * x_threshold)
        f, ax = plt.subplots()
        plt.scatter(np.log10(mu_genes), np.log10(fano), color = 'gray', alpha = 0.3, s = 3);
        plt.scatter(np.log10(mu_genes)[i_hvgs], np.log10(fano)[i_hvgs], color = 'black', alpha = 0.3, s = 3);
        plt.plot(np.log10(x_threshold), np.log10(y_threshold))
        ax.set_xlabel('log10(mean)')
        ax.set_ylabel('log10(fano factor)')

    return i_genes[i_hvgs]
