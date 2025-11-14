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
import time

import avaframe.out3Plot.plotUtils as pU
from avaframe.in3Utils import cfgUtils, initializeProject, logUtils
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


from sklearn.model_selection import KFold, cross_validate
from sklearn.preprocessing import StandardScaler
from sklearn.pipeline import Pipeline
from sklearn.metrics import make_scorer, mean_squared_error, r2_score, mean_absolute_error

from sklearn.gaussian_process import GaussianProcessRegressor
from sklearn.gaussian_process.kernels import Matern, WhiteKernel, ConstantKernel

from sklearn.ensemble import ExtraTreesRegressor  # robust tree surrogate (TPE-like behavior)

from scipy.stats import norm, qmc


def readParamSetDF(inDir, varParList, paramSelected=None):
    """
    Read parameter sets from .ini files in a directory and build a DataFrame.

    Parameters
    ----------
    inDir : str or pathlib.Path
        Path to directory containing .ini files.
    varParList : list of str
        List of parameter names to extract values from each .ini file.

    Returns
    -------
    paramSetDF : pandas.DataFrame
        DataFrame with simName, parameterSet and order as columns
    """

    # List to hold all parameters sets
    paramSet = []
    order = []
    sampleMethods = []
    filenames = []

    # Loop over all files in the folder
    for filename in os.listdir(inDir):
        # Select only .ini files that contain 'AnbruchDrone' in the name
        if filename.endswith('.ini') and 'sourceConfiguration' not in filename:
            filepath = os.path.join(inDir, filename)

            # Load the .ini file
            config = configparser.ConfigParser()
            config.read(filepath)

            if 'VISUALISATION' in config.sections():
                # config is inifile
                index = config['VISUALISATION']['scenario']

                if 'VISUALISATION' in config.sections():
                    # config is inifile
                    index = config['VISUALISATION']['scenario']
                    if 'sampleMethod' in config['VISUALISATION']:
                        sampleMethod = config['VISUALISATION']['sampleMethod']
                    else:
                        sampleMethod = np.nan

            row = []  # row contains 1 row
            for param in varParList:
                section = probAna.fetchParameterSection(config, param)
                value = config[section][param]
                value = float(value)
                row.append(value)

            order.append(index)
            sampleMethods.append(sampleMethod)
            paramSet.append(row)  # rows contains all rows
            filenames.append(os.path.splitext(filename)[0])

    # convert to pandas DF
    paramSetDF = pd.DataFrame({
        'simName': filenames,
        'parameterSet': paramSet,  # [row for row in paramSet], # Wrap each row as a list
        'order': pd.to_numeric(order),  # convert to int
        'sampleMethod': sampleMethods
    })

    # add selected parameters as columns if not None
    if paramSelected is not None:
        cols_to_add = paramSelected
    else:
        cols_to_add = varParList  # use all available parameter names

    dummy = pd.DataFrame(paramSetDF['parameterSet'].tolist(), columns=varParList)[cols_to_add]
    paramSetDF = pd.concat([paramSetDF, dummy], axis=1)


    return paramSetDF


def readArealIndicators(inDir):
    """
    Read areal indicator results from a pickle file and convert to a DataFrame.

    Parameters
    ----------
    inDir : str or pathlib.Path
        Path to pickle file containing indicator results.

    Returns
    -------
    indicatorsDF : pandas.DataFrame
        DataFrame with simName, and and areal indicators,
    """

    with open(inDir, "rb") as f:
        all_results = pickle.load(f)

    indicatorsDF = pd.DataFrame(all_results)
    return indicatorsDF


def readMaxVelocity(avaDir):
    '''

    Parameters
    ----------
    avaDir

    Returns
    -------
    df: pandas.Dataframe
    with u2_max and simName
    '''

    # Point to the folder that contains the text files
    folder = pathlib.Path(avaDir, "Outputs", "com8MoTPSA")
    results = []

    # Loop through all matching files
    for file_path in folder.glob("*_DataTime.txt"):
        # Read each file
        df = pd.read_csv(file_path, sep=r"\s+", header=0)
        df.columns = df.columns.str.lstrip('#')
        # Compute the max value of 's2_max(m)'
        s2_max_value = df['s2_max(m)'].max()
        # Extract the base name without the suffix
        simName = file_path.stem.replace("_DataTime", "")
        results.append({"simName": simName, "u2_max": s2_max_value})

    # Create a summary DataFrame
    return pd.DataFrame(results)


def addLossMetrics(df):
    """
    Nicht kpomplett, Compute evaluation metrics (recall, precision, F1, Tversky score) and an
    optimisation variable from a given DataFrame.

    The metrics are based on area (number of pixel would also be possible). Invalid values
    (division by zero) are replaced with 0.

    Parameters
    ----------
    df : pandas.DataFrame
        Input DataFrame with at least the columns
        ``TP_SimRef_area``, ``FP_SimRef_area``, ``FN_SimRef_area``.

    Returns
    -------
    df : pandas.DataFrame
        Same DataFrame with additional columns:
        - ``recall`` : float
        - ``precision`` : float
        - ``f1_score`` : float
        - ``tversky_score`` : float
        - ``optimisationVariable`` : float (1 - tversky_score)
    """

    # decide if loss function is based on ncells or area
    basedOn = '_area'

    TP = df[f"TP_SimRef{basedOn}"]
    FP = df[f"FP_SimRef{basedOn}"]
    FN = df[f"FN_SimRef{basedOn}"]

    # recall = TP / (TP + FN)
    denomRecall = TP + FN
    df["recall"] = np.where(denomRecall != 0, TP / denomRecall, 0.0)

    # precision = TP / (TP + FP)
    denomPrecision = TP + FP
    df["precision"] = np.where(denomPrecision != 0, TP / denomPrecision, 0.0)

    # f1 Score
    denomF1 = df['precision'] + df['recall']
    df["f1_score"] = np.where(denomF1 != 0, 2.0 * df['precision'] * df['recall'] / denomF1, 0.0)

    # Tversky score = TP / (TP + alpha * FP + beta * FN), gives penalty to overshoot --> alpha
    alpha = 2
    beta = 1
    denomTversky = TP + alpha * FP + beta * FN
    df["tversky_score"] = 1 - np.where(denomTversky != 0, TP / denomTversky,
                                       0.0)  # Subtract 1 to ensure that 0 are good values and 1 bad

    df['optimisationVariable'] = df['tversky_score']
    '''
    # u2_max
    u2_max_ref = 247.931  # from ref sim: "d1867f98e6"
    u2_diff = u2_max_ref - df['u2_max']

    # Min-Max-Normalisierung für u2_diff
    u2_diff_min = u2_diff.min()
    u2_diff_max = u2_diff.max()
    u2_diff_norm = (u2_diff - u2_diff_min) / (u2_diff_max - u2_diff_min)

    df['optimisationVariable'] = 0.5 * df['tversky_score'] + 0.5 * u2_diff_norm
    '''
    return df


def buildFinalDF(avalancheDir, cfgProb, paramSelected=None):
    """
    Build the final merged DataFrame for a given avalanche.

    Combines parameter sets, AIMEC results, and areal indicators into one DataFrame,
    then computes evaluation metrics 'via addLossMetrics'.

    Parameters
    ----------
    avalancheDir : str
        Path of avalanche directory
    cfgProb : configparser.ConfigParser
        Config parser
    Returns
    -------
    finalDF : pandas.DataFrame
        Final DataFrame containing:
        - ``simName``
        - ``parameterSet``
        - ``order``
        - Areal indicator columns
        - Evaluation metrics (recall, precision, f1_score, tversky_score, optimisationVariable)
    """
    # Load variables that a varied
    varParList = cfgProb['PROBRUN']['varParList'].split('|')

    # Folder where ini files from simulations are
    inDir = pathlib.Path(avalancheDir, 'Outputs/com8MoTPSA/configurationFiles')
    # read parameterSetDF
    paramSetDF = readParamSetDF(inDir, varParList, paramSelected)

    # dataframe from AIMEC
    # df_aimec = pd.read_csv(
    #    '../' + avalancheDir + '/Outputs/ana3AIMEC/com8MoTPSA/Results_' + avaName + '_ppr_lim_1_w_600resAnalysisDF.csv')
    # merge aimec with parameterSet
    # df_merged = pd.merge(paramSetDF, df_aimec, on='simName', how='inner')

    # Folder where pickle file is saved
    arealIndicatorDir = pathlib.Path(avalancheDir, 'Outputs', 'out1Peak', 'arealIndicators.pkl')
    # read areal indicators
    indicatorsDF = readArealIndicators(arealIndicatorDir)

    # read u2_max df
    # u2_maxDF = readMaxVelocity(avalancheDir)

    # merge df
    df_merged = pd.merge(paramSetDF, indicatorsDF, on='simName', how='inner')
    # df_merged = df_merged.merge(u2_maxDF, on='simName')

    # add optimisation variables
    finalDF = addLossMetrics(df_merged)

    return finalDF


def createDFParameterLoss(df, paramAll, paramSelected):
    """
    Create DataFrames linking selected parameters with the loss function.

    Expands the ``parameterSet`` column from the input DataFrame into a
    parameter table, selects the relevant parameters, and attaches the
    optimisation variable ("Loss"). Also returns a version where the
    parameters are min–max scaled to [0, 1].

    Parameters
    ----------
    df : pandas.DataFrame
    paramAll : list of str
        Names of all parameters in the same order as stored in ``parameterSet``.
    paramSelected : list of str
        Subset of parameters to include in the output DataFrames.

    Returns
    -------
    paramLossDF : pandas.DataFrame
        DataFrame with one column per selected parameter and an additional
        ``Loss`` column with the raw values of ``optimisationVariable``.
    paramLossDFScaled : pandas.DataFrame
        Same as ``paramLossDF`` but with the selected parameters normalised
        to the range [0, 1] using min–max scaling.
    """
    paramLossDF = df[paramSelected].copy()
    paramLossDFScaled = (paramLossDF - paramLossDF.min()) / (paramLossDF.max() - paramLossDF.min())  # normalise
    paramLossDF['Loss'] = df[
        'optimisationVariable']  # merging works with different index, but it is the right value
    paramLossDFScaled['Loss'] = df[
        'optimisationVariable']  # merging works with different index, but it is the right value
    return paramLossDF, paramLossDFScaled


def fitSurrogate(df):

    # Prepare X, y
    y_col = 'Loss'

    X = df.drop(columns=[y_col]).to_numpy(dtype=float)
    y = df[y_col].to_numpy(dtype=float).reshape(-1)

    n_features = X.shape[1]
    print(f"Data shape: X={X.shape}, y={y.shape}, features={list(df.drop(columns=[y_col]).columns)}")

    # ---- GP kernel (Matern-Kovarianz für glattere Funktionen)
    kernel = (
            ConstantKernel(1.0, (1e-3, 1e3))  # Output varianz, informationen wie stark Y variiert
            * Matern(length_scale=np.ones(n_features),
                     length_scale_bounds=(1e-2, 1e2),  # in z score Raum, also wie oft mal Standardabweichen
                     nu=2.5)
    )

    gp = GaussianProcessRegressor(
        kernel=kernel,
        alpha=1e-8,
        normalize_y=False,
        n_restarts_optimizer=10,
        random_state=0,
    )

    etr = ExtraTreesRegressor(
        n_estimators=500,
        min_samples_leaf=1,
        bootstrap=False,
        random_state=0,
        n_jobs=-1,
    )

    # ---- Pipelines (feature scaling + model)
    gp_pipe = Pipeline([("x_scaler", StandardScaler()), ("model", gp)])
    etr_pipe = Pipeline([("x_scaler", StandardScaler()), ("model", etr)])

    return X, y, gp_pipe, etr_pipe


def KFoldCV(X,y, pipe, pipeName):
    # For losses, sklearn uses "neg_*" because higher-is-better internally.
    rmse_scorer = "neg_root_mean_squared_error"
    mae_scorer = "neg_mean_absolute_error"
    r2_scorer = "r2"

    cv = KFold(n_splits=5, shuffle=True, random_state=0)

    scores = cross_validate(
        pipe, X, y, cv=cv,
        scoring={"rmse": rmse_scorer, "mae": mae_scorer, "r2": r2_scorer},
        return_train_score=True,
        error_score="raise"  # fail fast if something else is wrong
    )
    # NOTE: rmse/mae were returned as NEGATIVE numbers because the higher is better
    test_rmse = -scores['test_rmse']
    test_mae = -scores['test_mae']
    train_rmse = -scores['train_rmse']
    train_mae = -scores['train_mae']

    print(f"\n{pipeName} 5-fold CV:")
    print(f"  RMSE  (mean ± sd): {test_rmse.mean():.4g} ± {test_rmse.std():.4g}")
    print(f"  MAE   (mean ± sd): {test_mae.mean():.4g} ± {test_mae.std():.4g}")
    print(f"  Test R²    (mean ± sd): {scores['test_r2'].mean():.4g} ± {scores['test_r2'].std():.4g}")
    print(f"  Train RMSE        : {train_rmse.mean():.4g}")
    print(f"  Train R²          : {scores['train_r2'].mean():.4g}")

    return scores


def optimiseNonSeqV1(pipe, paramBounds):
    """
    Creates a number of samples, predicts loss function with surrogate. Return a dictionary with stats of the best N
    surrogate samples.
    Parameters
    ----------
    pipe :
    paramBounds:

    Returns
    -------
    topNStats : pandas.DataFrame
        DataFrame with simName, and and areal indicators,
    """

    paramSelected = list(paramBounds.keys())
    bounds = np.array(list(paramBounds.values()), dtype=float)  # shape (d,2)
    d = bounds.shape[0]

    # get LH canditates
    # sampler = qmc.LatinHypercube(d=d)
    sampler = qmc.LatinHypercube(d=d, seed=12345)
    sample = sampler.random(n=1000000)
    X0 = qmc.scale(sample, bounds[:, 0], bounds[:, 1])

    # Vorhersage des Verlusts mit GP-Modell
    mu, sigma = pipe.predict(X0, return_std=True)

    df_candidates = pd.DataFrame(X0, columns=paramSelected)

    topNStat, _ = analyzeTopCandidates(df_candidates, mu, sigma, paramSelected, N=100)
    return topNStat


def analyzeTopCandidates(df_candidates, mu, sigma, param_cols, N=5):
    """
    Analysis of top N candidates and return statistic with topN and best surrogate values.

    """

    # --- Top N ---
    idx_topN = np.argsort(mu)[:N]
    topNData = df_candidates.iloc[idx_topN].copy()
    topNData["mu"] = mu[idx_topN]
    topNData["sigma"] = sigma[idx_topN]

    mean_params = topNData[param_cols].mean()
    std_params = topNData[param_cols].std()
    mean_mu = topNData["mu"].mean()
    std_mu = topNData["mu"].std()
    mean_sigma = topNData["sigma"].mean()
    std_sigma = topNData["sigma"].std()

    print(f"\n🔍 Mittelwerte ± Std (Top {N}):")
    for p in param_cols:
        m, s = mean_params[p], std_params[p]
        perc = (s / m * 100) if m != 0 else np.nan
        print(f"  {p:30s}: {m:.6f} ± {s:.6f} ({perc:.1f}%)")
    perc_mu = (std_mu / mean_mu * 100) if mean_mu != 0 else np.nan
    perc_sigma = (std_sigma / mean_sigma * 100) if mean_sigma != 0 else np.nan
    print(f"📉 mu:    {mean_mu:.4f} ± {std_mu:.4f} ({perc_mu:.1f}%)")
    print(f"📊 sigma: {mean_sigma:.4f} ± {std_sigma:.4f} ({perc_sigma:.1f}%)")

    # --- Bester einzelner Punkt ---
    idx_best = np.argmin(mu)
    best_params = df_candidates.iloc[idx_best].copy()
    best_loss = mu[idx_best]
    best_sigma = sigma[idx_best]

    print("\n🔍 Beste einzelne Parameterkombination laut GP:")
    for p in param_cols:
        print(f"  {p:30s}: {best_params[p]:.4f}")
    print(f"📉 mu:    {best_loss:.4f}")
    print(f"📊 sigma: {best_sigma:.4f}")

    return {
        f"TopNBest": {
            "mean_params": mean_params,
            "std_params": std_params,
            "mean_mu": mean_mu,
            "std_mu": std_mu,
            "mean_sigma": mean_sigma,
            "std_sigma": std_sigma,
        },
        "Best": {
            "params": best_params,
            "mu": best_loss,
            "sigma": best_sigma,
        }
    }, topNData


def saveTopCandidates(results_dict, finalDF, paramSelected, out_path="analysisTable.png", title=None, simName=None):
    """
    Create two result tables (Surrogate and Model), save them as PNG and CSV.
    Both tables are plotted into a single PNG image (stacked) and exported into
    one CSV file with section headers.

    Parameters
    ----------
    results_dict : dict
        Dictionary containing surrogate optimisation results with keys:
          - "TopNBest": dict with entries "mean_params", "std_params",
            "mean_mu", "std_mu", "mean_sigma", "std_sigma".
          - "Best": dict with entries "params", "mu", "sigma".
    finalDF : pandas.DataFrame
        DataFrame with all model results. Must contain:
          - Column "optimisationVariable"
          - Column "simName"
          - Parameter columns referenced by paramSelected.
    paramSelected : list of str
        List of parameter column names to analyse.
    out_path : str or pathlib.Path, optional
        Path to save the PNG image. Default = "analysisTable.png".
        The CSV file will be saved in the same folder with suffix "_tables.csv".
    title : str, optional
        Global title for the PNG figure. Default = None.
    simName : str, optional
        Simulation name to highlight in Table 1 (Surrogate). If provided,
        its optimisation variable is added as a row.

    Returns
    -------
    pathlib.Path
        Path to the saved PNG figure.

    Notes
    -----
    - This code was written with AI.
    """
    # --- Table 1: TopNBest (Surrogate) + Best ---
    top = results_dict["TopNBest"]
    mean_params = pd.Series(top["mean_params"])
    std_params = pd.Series(top["std_params"])
    rel_params = std_params / mean_params * 100.0

    df_top = pd.DataFrame({"mean": mean_params, "std": std_params, "relStd [%]": rel_params})
    df_top = pd.concat([
        df_top,
        pd.DataFrame(
            {"mean": [top["mean_mu"], top["mean_sigma"]],
             "std": [top["std_mu"], top["std_sigma"]],
             "relStd [%]": [top["std_mu"] / top["mean_mu"] * 100.0,
                            top["std_sigma"] / top["mean_sigma"] * 100.0]},
            index=["optimisationVariableSurrogate", "sigma"]  # formerly "mu"
        )
    ])

    if simName is not None:
        val = finalDF.loc[finalDF["simName"] == simName, "optimisationVariable"].iloc[0]
        df_top.loc["optimisationVariable", ["mean", "std", "relStd [%]"]] = [val, None, None]

    best = results_dict["Best"]
    df_best = pd.Series(best["params"], name="best")
    df_best.loc["optimisationVariableSurrogate"] = best["mu"]
    df_best.loc["sigma"] = best["sigma"]
    df_top["best"] = df_best

    # --- Table 2: Top10 (Model) + Best + optimisationVariable ---
    best_idx = finalDF["optimisationVariable"].idxmin()
    top10 = finalDF.nsmallest(10, "optimisationVariable")

    top10_params_num = top10[paramSelected].apply(pd.to_numeric, errors="coerce")
    df_top10 = pd.DataFrame({
        "mean": top10_params_num.mean(),
        "std": top10_params_num.std(),
        "relStd [%]": top10_params_num.std() / top10_params_num.mean() * 100.0,
        "best": pd.to_numeric(finalDF.loc[best_idx, paramSelected], errors="coerce"),
    })

    opt_mean = pd.to_numeric(top10["optimisationVariable"], errors="coerce").mean()
    opt_std = pd.to_numeric(top10["optimisationVariable"], errors="coerce").std()
    df_top10.loc["optimisationVariable", ["mean", "std", "relStd [%]", "best"]] = [
        opt_mean, opt_std, opt_std / opt_mean * 100.0, finalDF.loc[best_idx, "optimisationVariable"]
    ]

    # --- Titles ---
    title_1 = f"Surrogate: Top N Best + Single Best {f'({simName})' if simName else ''}"
    best_sim = finalDF.at[best_idx, "simName"]
    title_2 = f"Model: Top 10 Best + Single Best {f'({best_sim})' if best_sim else ''}"

    # --- Number formatting ---
    def formatSig(x):
        try:
            x = float(x)
        except (TypeError, ValueError):
            return ""
        if pd.isna(x):
            return ""
        if x == 0.0:
            return "0"
        if x < 1e-3:
            return f"{x:.2e}"  # small numbers in scientific notation
        elif x < 1000:
            s = f"{x:.2g}"  # 2 significant digits
            if "e" in s or "E" in s:
                return f"{x:.6f}".rstrip("0").rstrip(".")
            return s
        else:
            return str(int(x))  # large numbers as integers

    # --- Apply formatting to tables ---
    def fmt_df(df):
        g = df.copy()
        for c in ["mean", "std", "relStd [%]", "best"]:
            if c in g.columns:
                g[c] = g[c].map(formatSig)
        return g

    df_top_disp = fmt_df(df_top)
    df_top10_disp = fmt_df(df_top10)

    # --- Plot with two tables ---
    fig_h = 1.2 + 0.38 * (len(df_top_disp) + len(df_top10_disp))
    fig, axes = plt.subplots(2, 1, figsize=(10, fig_h))
    if title:
        fig.suptitle(title, fontsize=14, y=0.99)

    for ax, df_disp, t in [(axes[0], df_top_disp, title_1), (axes[1], df_top10_disp, title_2)]:
        ax.axis("off")
        tbl = ax.table(
            cellText=df_disp.values,
            rowLabels=df_disp.index.tolist(),
            colLabels=df_disp.columns.tolist(),
            loc="center"
        )
        tbl.auto_set_font_size(False)
        tbl.set_fontsize(9)
        tbl.scale(1, 1.2)
        ax.set_title(t, fontsize=12, pad=10)

    plt.tight_layout(rect=(0, 0, 1, 0.98))

    # --- Save PNG ---
    out_path = pathlib.Path(out_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=200, bbox_inches="tight")
    plt.close(fig)

    # --- Save CSV (both tables stacked) ---
    csv_path = out_path.with_name(f"{out_path.stem}_tables.csv")
    cols = ["row", "mean", "std", "relStd [%]", "best"]

    t1 = fmt_df(df_top).reset_index(names=["row"])[cols]
    t2 = fmt_df(df_top10).reset_index(names=["row"])[cols]

    def header(txt):
        return pd.DataFrame([{"row": txt, "mean": "", "std": "", "relStd [%]": "", "best": ""}])

    csv_both = pd.concat([header(title_1), t1, header(title_2), t2], ignore_index=True)
    csv_both.to_csv(csv_path, index=False)

    return out_path



def expectedImprovement(mu, sigma, f_best, xi=0.02):
    sigma = np.maximum(sigma, 1e-12)  # numeric safety
    imp = f_best - mu - xi  # minimization, thats why the sign is different, xi for finetunig exploitation

    Z = imp / sigma
    ei = imp * norm.cdf(Z) + sigma * norm.pdf(Z)
    ei[sigma <= 1e-12] = 0.0  # set EI to zero where sigma is 1e-12
    return ei

    # https://ekamperi.github.io/machine%20learning/2021/06/11/acquisition-functions.html


def lowerConfidenceBound(mu, sigma, k=2.0):
    return -mu + k * sigma  # then do: x_next = X0[np.argmax(lcb)]


def EINextPoint(pipe, y, paramBounds):
    paramSelected = list(paramBounds.keys())
    bounds = np.array(list(paramBounds.values()), dtype=float)  # shape (d,2)
    d = bounds.shape[0]
    f_best = np.nanmin(y)

    # get LH canditates
    sampler = qmc.LatinHypercube(d=d)
    sample = sampler.random(n=100000)
    X0 = qmc.scale(sample, bounds[:, 0], bounds[:, 1])

    # predict with pipe
    mu, sigma = pipe.predict(X0, return_std=True)

    print(mu.mean(), mu.max(), sigma.mean(), sigma.max())

    # EI or LCB for minimization
    ei = expectedImprovement(mu, sigma, f_best)
    lcb = lowerConfidenceBound(mu, sigma)

    xBest = X0[np.argmax(ei)].copy()
    # xBest = X0[np.argmax(lcb)].copy()
    xBestDict = {feat: float(val) for feat, val in zip(paramSelected, xBest)}

    return xBest, xBestDict, np.max(ei), np.max(lcb)


def runCom8MoTPSA(avalancheDir, xBestDict, cfgMain, i=0, optimisationType=None):
    # Time the whole routine
    startTime = time.time()

    # log file name; leave empty to use default runLog.log
    logName = 'runCom8MoTPSA'

    # Start logging
    log = logUtils.initiateLogger(avalancheDir, logName)
    log.info('MAIN SCRIPT')
    log.info('Current avalanche: %s', avalancheDir)
    # ----------------
    # Clean input directory(ies) of old work and output files
    # If you just created the ``avalancheDir`` this one should be clean but if you
    # already did some calculations you might want to clean it::
    initializeProject.cleanSingleAvaDir(avalancheDir, deleteOutput=False)
    # Get module config
    cfgCom8MoTPSA = cfgUtils.getModuleConfig(com8MoTPSA, toPrint=False)

    # overwrite cfgcom8 with xBest values
    for param, val in xBestDict.items():
        # print(param, val)
        section = probAna.fetchParameterSection(cfgCom8MoTPSA, param)
        cfgCom8MoTPSA[section][param] = str(val)
    # give visualisation unique scenario for identifying later
    cfgCom8MoTPSA['VISUALISATION']['scenario'] = str(i)
    if optimisationType == 'nonSeq':
        cfgCom8MoTPSA["VISUALISATION"]["sampleMethod"] = 'nonSeq'
    else:
        cfgCom8MoTPSA["VISUALISATION"]["sampleMethod"] = 'EI/LCB'

    # ----------------
    # Run psa
    simName = com8MoTPSA.com8MoTPSAMain(cfgMain, cfgInfo=cfgCom8MoTPSA, returnSimName=True)
    # Print time needed
    endTime = time.time()
    log.info('Took %6.1f seconds to calculate.' % (endTime - startTime))

    return simName


def saveBestRow(df, y, ei=None, lcb=None, simName=None, csv_path='dummy.csv'):
    '''
    If simName = none, dann wird beste bestehende sim gespeichert
    if simName = not None, dann wird zuletzt ausgeführte sim gespeichert

    '''

    if simName is not None:

        row = df.loc[df['simName'] == simName]
    else:
        idx = df[y].idxmin()
        row = df.loc[[idx]].copy()

    # ensure the optional columns always exist
    if ei is not None:
        row["ei"] = ei
    if lcb is not None:
        row["lcb"] = lcb

    path = pathlib.Path(csv_path)
    row.to_csv(path, mode="a", index=False, header=not path.exists())
