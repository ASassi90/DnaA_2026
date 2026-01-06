from scipy.optimize import curve_fit
import math
import numpy as np
from matplotlib import pyplot as plt

def hill_function(x, H, K):
    return x**H / (K**H + x**H)

def get_Hill(x_data, y_data, p0):
    """
        computes Hill coefficient and activation threshold from the fit of x_data and 
        y_data with a hill function.
    """
    popt, _ = curve_fit(hill_function, x_data, y_data, p0=p0)
    return popt[0], popt[1]

def get_c(dnaa, K, c_tot):
    """
        compute the concentration of occupied binding sites.
    """
    c=((dnaa+K+c_tot)-math.sqrt((dnaa+K+c_tot)**2.-4.*dnaa*c_tot))/2.
    return c

def create_figure(layout='single', figsize=(4,4), xlabel=None, ylabel=None, n_stacked=2, sharex=True, sharey=False, 
                  xlim=None, ylim=None, xticks=None, yticks=None, labelsize=14):
    """
    Creates a figure with a consistent design for scientific papers.
    
    Parameters:
    - layout: 'single', 'stacked', or 'grid' (2x2)
    - figsize: size of the square figure in inches
    - xlabel: label(s) for the x-axis (must be a list if layout='grid')
    - ylabel: label(s) for the y-axis (must be a list if layout='stacked' or 'grid')
    - n_stacked: number of stacked plots (only applies to 'stacked' layout)
    - sharex: whether to share the x-axis in stacked plots
    - sharey: whether to share the y-axis in stacked plots
    - xlim: range(s) for the x-axis (must be a list of tuples if layout='grid', single tuple if layout='stacked')
    - ylim: range(s) for the y-axis (must be a list of tuples matching panel count)
    - xticks: tick values for the x-axis (must be a list of lists if layout='grid', single list if layout='stacked')
    - yticks: tick values for the y-axis (must be a list of lists matching panel count)
    
    Returns:
    - fig, ax: Matplotlib figure and axes objects
    """
    plt.rcParams.update({
        'font.size': 12,       # Global font size
        'axes.labelsize': labelsize,  # Label font size
        'xtick.labelsize': 12, # X tick font size
        'ytick.labelsize': 12, # Y tick font size
        'axes.linewidth': 1.5,   # Adjusted Axes border thickness
        'xtick.major.width': 1.5, # Adjusted Tick thickness
        'ytick.major.width': 1.5,
        "text.usetex": False,  # Do not use external LaTeX
        "font.family": "serif",
        "mathtext.fontset": "stix",  # STIX fonts look similar to LaTeX
    })
    
    fig, ax = None, None
    
    if layout == 'single':
        fig, ax = plt.subplots(figsize=figsize)
        if isinstance(xlabel, str):
            ax.set_xlabel(xlabel)
        if isinstance(ylabel, str):
            ax.set_ylabel(ylabel)
        if xlim:
            ax.set_xlim(xlim)
        if ylim:
            ax.set_ylim(ylim)
        if xticks:
            ax.set_xticks(xticks)
        if yticks:
            ax.set_yticks(yticks)
    
    elif layout == 'stacked':
        fig, axes = plt.subplots(n_stacked, 1, figsize=figsize, sharex=sharex, sharey=sharey)
        #fig.subplots_adjust(hspace=0.3)  # Adjust vertical spacing to prevent overlap
        if ylabel is None or not isinstance(ylabel, list) or len(ylabel) != n_stacked:
            ylabel = [''] * n_stacked
        
        for i, ax in enumerate(axes):
            ax.tick_params(axis='both', width=1.5, labelsize=12)
            ax.set_ylabel(ylabel[i])
            if xlim:
                ax.set_xlim(xlim)  # Same x-axis range for all stacked panels
            if ylim and isinstance(ylim, list) and len(ylim) == n_stacked:
                ax.set_ylim(ylim[i])
            if xticks:
                ax.set_xticks(xticks)  # Same xticks for all stacked panels
            if yticks and isinstance(yticks, list) and len(yticks) == n_stacked:
                ax.set_yticks(yticks[i])
        
        axes[-1].set_xlabel(xlabel)
        ax = axes

    elif layout == 'grid':
        fig = plt.figure(figsize=figsize)
        
        # Define fixed panel positions (x0, y0, x1, y1) in figure coordinates
        panel_positions = [
            [0.13, 0.58, 0.32, 0.32],  # Top-left
            [0.58, 0.58, 0.32, 0.32],  # Top-right
            [0.13, 0.13, 0.32, 0.32],  # Bottom-left
            [0.58, 0.13, 0.32, 0.32],  # Bottom-right
        ]
        
        axes = []
        for pos in panel_positions:
            ax = fig.add_axes(pos)
            #ax.set_box_aspect(1)  # Preserve square aspect ratio
            ax.tick_params(axis='both', width=1.5, labelsize=12)
            axes.append(ax)
        
        if xlabel:
            axes[2].set_xlabel(xlabel[2])
            axes[3].set_xlabel(xlabel[3])
            axes[0].set_xlabel(xlabel[0])
            axes[1].set_xlabel(xlabel[1])
        if ylabel:
            axes[0].set_ylabel(ylabel[0])
            axes[1].set_ylabel(ylabel[1])
            axes[2].set_ylabel(ylabel[2])
            axes[3].set_ylabel(ylabel[3])
        
        ax=np.array(axes).reshape(2, 2)
    
    else:
        raise ValueError("Invalid layout. Choose 'single', 'stacked', or 'grid'.")
    
    return fig, ax