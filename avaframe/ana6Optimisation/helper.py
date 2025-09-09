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

    # add selected parameters as columns
    dummy = pd.DataFrame(paramSetDF['parameterSet'].tolist(), columns=varParList)[paramSelected]
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


def addLossMetrics(df):
    """
    Compute evaluation metrics (recall, precision, F1, Tversky score) and an
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
    df["tversky_score"] = np.where(denomTversky != 0, TP / denomTversky, 0.0)

    df['optimisationVariable'] = 1 - df['tversky_score']  # Subtract 1 to ensure that 0 are good values and 1 bad
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

    # merge df
    df_merged = pd.merge(paramSetDF, indicatorsDF, on='simName', how='inner')

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
    paramLossDF = df[paramSelected]
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
    # n_features must be defined elsewhere to match your X
    kernel = (
            ConstantKernel(1.0, (1e-6, 1e6))
            * Matern(length_scale=np.ones(n_features),
                     length_scale_bounds=(1e-3, 1e6),
                     nu=1.5)
            + WhiteKernel(noise_level=1e-3, noise_level_bounds=(1e-8, 1e2))
    )

    gp = GaussianProcessRegressor(
        kernel=kernel,
        alpha=0.0,
        normalize_y=True,
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
    sampler = qmc.LatinHypercube(d=d, seed=12345)
    sample = sampler.random(n=1000000)
    X0 = qmc.scale(sample, bounds[:, 0], bounds[:, 1])

    df_candidates = pd.DataFrame(X0, columns=paramSelected)

    # Vorhersage des Verlusts mit GP-Modell
    mu, sigma = pipe.predict(df_candidates, return_std=True)

    topNStat, _ = analyzeTopCandidates(df_candidates, mu, sigma, paramSelected, N=100)

    return topNStat


def analyzeTopCandidates(df_candidates, mu, sigma, param_cols, N=100):
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
    Speichert drei Tabellen als ein Bild:
      1) TopNBest (Mittelwert/Std/Rel.-Std der Parameter + mu/sigma + optimisationVariable der gewählten Simulation)
      2) Best (beste Parameterkombination + mu/sigma)
      3) FinalDF-Eintrag (beste Zeile nach optimisationVariable) transponiert: ausgewählte Parameter + optimisationVariable
    """
    # --- TopNBest-DataFrame bauen ---
    top = results_dict["TopNBest"]
    mean_params = pd.Series(top["mean_params"]).astype(float)
    std_params = pd.Series(top["std_params"]).astype(float)
    rel_params = (std_params / mean_params.replace(0, np.nan) * 100.0)

    df_top = pd.DataFrame({
        "mean": mean_params,
        "std": std_params,
        "rel_std_%": rel_params
    })

    # mu/sigma Zeilen anhängen
    mean_mu = float(top["mean_mu"])
    std_mu = float(top["std_mu"])
    mean_sigma = float(top["mean_sigma"])
    std_sigma = float(top["std_sigma"])

    rel_mu = (std_mu / mean_mu * 100.0) if mean_mu != 0 else np.nan
    rel_sigma = (std_sigma / mean_sigma * 100.0) if mean_sigma != 0 else np.nan

    df_top = pd.concat([
        df_top,
        pd.DataFrame({
            "mean": [mean_mu, mean_sigma],
            "std": [std_mu, std_sigma],
            "rel_std_%": [rel_mu, rel_sigma]
        }, index=["mu", "sigma"])
    ])

    # zusätzliche Zeile: optimisationVariable der gewünschten Simulation (aus finalDF gefiltert per simName)
    if simName is not None:
        match = finalDF.loc[finalDF["simName"] == simName, "optimisationVariable"]
        if not match.empty:
            df_top.loc[f"optimisationVariable ({simName})", ["mean", "std", "rel_std_%"]] = [float(match.iloc[0]),
                                                                                             np.nan, np.nan]

    # --- Best-DataFrame bauen ---
    best = results_dict["Best"]
    best_params = pd.Series(best["params"]).astype(float)
    df_best = pd.DataFrame({"value": best_params})
    df_best.loc["mu", "value"] = float(best["mu"])
    df_best.loc["sigma", "value"] = float(best["sigma"])

    # --- FinalDF (beste Zeile nach optimisationVariable), transponiert ---
    if isinstance(paramSelected, (list, tuple, pd.Index, np.ndarray)):
        param_cols = list(paramSelected)
    else:
        param_cols = [paramSelected]
    need = param_cols + ["optimisationVariable"]

    best_idx = finalDF["optimisationVariable"].idxmin()
    row = finalDF.loc[best_idx, need]
    df_final = pd.DataFrame({"value": row})

    # --- Formatierhilfe ---
    def fmt_df(df, cols=None, nd=6):
        df = df.copy()
        cols = df.columns if cols is None else cols
        for c in cols:
            orig = df[c]
            nums = pd.to_numeric(orig, errors="coerce")
            fmt_vals = nums.map(lambda x: f"{x:.{nd}g}" if pd.notnull(x) else None)
            df[c] = np.where(nums.notnull(), fmt_vals, orig.astype(str))
        return df

    df_top_disp = fmt_df(df_top, ["mean", "std", "rel_std_%"], nd=6)
    df_best_disp = fmt_df(df_best, ["value"], nd=6)
    df_final_disp = fmt_df(df_final, ["value"], nd=6)

    # --- Plot: drei Tabellen ---
    n_rows_top = len(df_top_disp)
    n_rows_best = len(df_best_disp)
    n_rows_final = len(df_final_disp)
    fig_h = 1.2 + 0.35 * (n_rows_top + n_rows_best + n_rows_final)

    fig, axes = plt.subplots(3, 1, figsize=(10, fig_h))
    if title:
        fig.suptitle(title, fontsize=14, y=0.99)

    axes[0].axis("off")
    tbl1 = axes[0].table(
        cellText=df_top_disp.values,
        rowLabels=df_top_disp.index.tolist(),
        colLabels=df_top_disp.columns.tolist(),
        loc="center"
    )
    tbl1.auto_set_font_size(False);
    tbl1.set_fontsize(9);
    tbl1.scale(1, 1.2)
    axes[0].set_title("Surrogate TopNBest – Parameterstatistik & μ/σ", fontsize=12, pad=10)

    axes[1].axis("off")
    tbl2 = axes[1].table(
        cellText=df_best_disp.values,
        rowLabels=df_best_disp.index.tolist(),
        colLabels=df_best_disp.columns.tolist(),
        loc="center"
    )
    tbl2.auto_set_font_size(False);
    tbl2.set_fontsize(10);
    tbl2.scale(1, 1.2)
    axes[1].set_title("Surrogate Best – Best Parameterkombination ", fontsize=12, pad=10)

    axes[2].axis("off")
    tbl3 = axes[2].table(
        cellText=df_final_disp.values,
        rowLabels=df_final_disp.index.tolist(),
        colLabels=df_final_disp.columns.tolist(),
        loc="center"
    )
    tbl3.auto_set_font_size(False);
    tbl3.set_fontsize(9);
    tbl3.scale(1, 1.2)
    axes[2].set_title("FinalDF – Best Model Run", fontsize=12, pad=10)

    plt.tight_layout(rect=(0, 0, 1, 0.98))
    out_path = pathlib.Path(out_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=200, bbox_inches="tight")
    plt.close(fig)
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

    # ensure the optional columns always exist (prevents CSV schema drift)
    if ei is not None:
        row["ei"] = ei
    if lcb is not None:
        row["lcb"] = lcb

    path = pathlib.Path(csv_path)
    row.to_csv(path, mode="a", index=False, header=not path.exists())
