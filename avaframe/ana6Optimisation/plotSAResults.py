import numpy as np
import logging
import pathlib
from scipy.stats import qmc
from SALib.sample import morris
from SALib.analyze import morris as morris_analyze
import pickle

import pandas as pd
import configparser
import os
import matplotlib.pyplot as plt
import seaborn as sns
import rasterio
from rasterio.mask import mask
import geopandas as gpd
import json
import re
from adjustText import adjust_text

import avaframe.out3Plot.plotUtils as pU
from avaframe.in3Utils import cfgUtils
from avaframe.in3Utils import logUtils
from avaframe.in3Utils import cfgHandling
from avaframe.in3Utils import fileHandlerUtils as fU
import avaframe.in2Trans.rasterUtils as IOf
import avaframe.in1Data.computeFromDistribution as cP
import avaframe.com1DFA.deriveParameterSet as dP
from avaframe.in3Utils import geoTrans as gT
from avaframe.out3Plot import statsPlots as sP
from avaframe.in1Data import getInput as gI

import avaframe.out1Peak.outPlotAllPeakDiffs as oPD
import avaframe.in2Trans.shpConversion as shpConv
import avaframe.com1DFA.DFAtools as DFAtls


from avaframe.com1DFA import com1DFA
from avaframe.com8MoTPSA import com8MoTPSA
from avaframe.ana4Stats import probAna
from datetime import datetime




def barplot(SiDF, avaName, outDir):
    """
    Create a bar plot of Morris sensitivity results.

    Bars show μ* as percentage of the total sensitivity, with bar width
    proportional to σ. A vertical dashed line is drawn where the cumulative μ*
    first reaches 80%.

    Parameters
    ----------
    SiDF : pandas.DataFrame
        DataFrame with at least the following columns:
        - ``Parameter`` : str, parameter names
        - ``mu_star`` : float, mean absolute elementary effect
        - ``mu_star_conf`` : float, optional, confidence interval of μ*
        - ``sigma`` : float, standard deviation of elementary effects
    avaName : str
        Avalanche name, used in the saved filename.
    outDir : str or pathlib.Path
        Directory where the plot is saved.
    """

    # 1) Sort by mu_star (descending)
    df = SiDF.sort_values("mu_star", ascending=False).reset_index(drop=True)

    # 2) Normalize μ* so the percentages sum to 100
    total_mu = df["mu_star"].sum()
    mu_pct = 100 * df["mu_star"] / total_mu

    # (optional) scale the confidence interval to the same percentage units
    if "mu_star_conf" in df:
        mu_pct_conf = 100 * df["mu_star_conf"] / total_mu
    else:
        mu_pct_conf = None

    # 3) Map σ to bar widths
    wmin, wmax = 0.3, 0.9
    rng = df["sigma"].max() - df["sigma"].min()
    sigma_norm = (df["sigma"] - df["sigma"].min()) / (rng if rng else 1)
    bar_widths = wmin + (wmax - wmin) * sigma_norm

    # 4) Plot
    x = np.arange(len(df))
    fig, ax = plt.subplots(figsize=(12, 6))
    bars = ax.bar(
        x, mu_pct, width=bar_widths, edgecolor="black",
        #yerr=mu_pct_conf if mu_pct_conf is not None else None,
        capsize=3 if mu_pct_conf is not None else 0
    )

    ax.set_xticks(x, df["Parameter"], rotation=60, ha="right")
    ax.set_ylabel("μ* (% of total)")
    ax.set_title("Sensitivity: μ* (percent of total) with σ mapped to bar width")
    #ax.set_ylim(0, max(105, mu_pct.max() * 1.15))  # headroom

    # Show value labels
    for xi, yi in zip(x, mu_pct):
        ax.text(xi, yi, f"{yi:.1f}%", ha="center", va="bottom", fontsize=9)


    # 5) Find where cumulative μ* first reaches 80%
    cum = mu_pct.cumsum()                    # cumulative percent (will end at 100)
    idx80 = int(np.searchsorted(cum, 80.0, side="left"))
    idx80 = min(idx80, len(df) - 1)          # safety: clamp to last if needed

    # Position of the right edge of the bar that brings us to >=80%
    boundary_x = x[idx80] + bar_widths[idx80] / 2

    # 6) Draw the vertical line and annotate
    ax.axvline(boundary_x, linestyle="--")
    ax.text(
        boundary_x, ax.get_ylim()[1] * 0.95,
        f"cum = {cum.iloc[idx80]:.1f}%)",
        rotation=90, va="top", ha="right"
    )

    # 7) Save figure
    # include date, format: YYYYMMDD
    date = datetime.now().strftime("%Y%m%d")
    figName = f"{outDir}/{avaName}_paramRanking_{date}.png"
    plt.savefig(figName, dpi=300, bbox_inches="tight")



def scatterplot(SiDF, avaName, outDir):
    """
    Create a scatter plot of Morris sensitivity results (μ* vs σ).

    Parameters
    ----------
    SiDF : pandas.DataFrame
        DataFrame with at least the following columns:
        - ``Parameter`` : str, parameter names
        - ``mu_star`` : float, mean absolute elementary effect
        - ``sigma`` : float, standard deviation of elementary effects
    avaName : str
        Avalanche name, used in the saved filename.
    outDir : str or pathlib.Path
        Directory where the plot is saved.
    """

    # Scatter Plot erstellen
    plt.figure(figsize=(12, 7))
    plt.scatter(SiDF['mu_star'], SiDF['sigma'], color='blue')

    # Punkte mit Namen beschriften
    for i, txt in enumerate(SiDF['Parameter']):
        plt.annotate(txt, (SiDF['mu_star'][i], SiDF['sigma'][i]), fontsize=9, xytext=(5, 5), textcoords='offset points')

    # Achsen und Titel
    plt.xlabel('mu_star (Einflussstärke)', fontsize=12)
    plt.ylabel('sigma (Nichtlinearität / Interaktionen)', fontsize=12)
    plt.title('Morris Sensitivitätsanalyse: mu_star vs sigma', fontsize=14)
    plt.grid(True)

    # Save figure
    # include date, format: YYYYMMDD
    date = datetime.now().strftime("%Y%m%d")
    figName = f"{outDir}/{avaName}_Scatterplot_{date}.png"
    plt.savefig(figName, dpi=300, bbox_inches="tight")



def scatterplotUncertainty(SiDF, avaName, outDir):
    """
    Create a scatter plot of Morris sensitivity results with uncertainty.
    Plots μ* vs σ with horizontal error bars given by ``mu_star_conf``.

    Parameters
    ----------
    SiDF : pandas.DataFrame
        DataFrame with at least the following columns:
        - ``Parameter`` : str, parameter names
        - ``mu_star`` : float, mean absolute elementary effect
        - ``mu_star_conf`` : float, confidence interval of μ*
        - ``sigma`` : float, standard deviation of elementary effects
    avaName : str
        Avalanche name, used in the saved filename.
    outDir : str or pathlib.Path
        Directory where the plot is saved.
    """
    # Plot with error bars
    plt.figure(figsize=(12, 7))
    plt.errorbar(
        SiDF['mu_star'], SiDF['sigma'],
        xerr=SiDF['mu_star_conf'],
        fmt='o', color='blue', ecolor='gray', elinewidth=1.5, capsize=4
    )

    # Annotations with adjustText
    texts = [plt.text(SiDF['mu_star'][i], SiDF['sigma'][i], SiDF['Parameter'][i], fontsize=9) for i in range(len(SiDF))]
    adjust_text(texts, arrowprops=dict(arrowstyle='->', color='gray', lw=0.5))

    # Axes and layout
    plt.xlabel('mu_star (Einflussstärke)', fontsize=12)
    plt.ylabel('sigma (Nichtlinearität / Interaktionen)', fontsize=12)
    plt.title('Morris Sensitivitätsanalyse: mu_star vs sigma (mit Unsicherheit)', fontsize=14)

    # Save figure
    # include date, format: YYYYMMDD
    date = datetime.now().strftime("%Y%m%d")
    figName = f"{outDir}/{avaName}_ScatterplotUncertainty_{date}.png"
    plt.savefig(figName, dpi=300, bbox_inches="tight")











