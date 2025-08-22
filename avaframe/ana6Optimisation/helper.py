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


from sklearn.model_selection import KFold, cross_validate
from sklearn.preprocessing import StandardScaler
from sklearn.pipeline import Pipeline
from sklearn.metrics import make_scorer, mean_squared_error, r2_score, mean_absolute_error

from sklearn.gaussian_process import GaussianProcessRegressor
from sklearn.gaussian_process.kernels import Matern, WhiteKernel, ConstantKernel

from sklearn.ensemble import ExtraTreesRegressor  # robust tree surrogate (TPE-like behavior)

from pyDOE import lhs




def readParamSetDF(inDir, varParList):
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
    # List to hold the corresponding filenames
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
            row = []  # row contains 1 row
            for param in varParList:
                section = probAna.fetchParameterSection(config, param)
                value = config[section][param]
                value = float(value)
                row.append(value)

            order.append(index)
            paramSet.append(row)  # rows contains all rows
            filenames.append(os.path.splitext(filename)[0])

    # convert to pandas DF
    paramSetDF = pd.DataFrame({
        'simName': filenames,
        'parameterSet': paramSet,  # [row for row in paramSet], # Wrap each row as a list
        'order': pd.to_numeric(order)  # convert to int
    })
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

    rows = []
    for entry in all_results:
        # Remove _ppr at the end if it exists
        clean_name = re.sub(r"_ppr$", "", entry["sim_name"])
        indicators = entry["indicator_dict"]
        row = {"simName": clean_name}

        for key, subdict in indicators.items():
            # Store short names: e.g., TP_cells, TP_area
            short_key = key.replace("truePositive", "TP_SimRef") \
                .replace("falsePositive", "FP_SimRef") \
                .replace("falseNegative", "FN_SimRef")
            row[f"{short_key}_cells"] = subdict.get("nCells", None)
            row[f"{short_key}_area"] = subdict.get("areaSum", None)
        rows.append(row)

    indicatorsDF = pd.DataFrame(rows)
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


def buildFinalDF(avaName, cfgProb):
    """
    Build the final merged DataFrame for a given avalanche.

    Combines parameter sets, AIMEC results, and areal indicators into one DataFrame,
    then computes evaluation metrics 'via addLossMetrics'.

    Parameters
    ----------
    avaName : str
        ava + Avalanche name
    cfgProb : configparser.ConfigParser
        Config parser
    Returns
    -------
    finalDF : pandas.DataFrame
        Final DataFrame containing:
        - ``simName``
        - ``parameterSet``
        - ``order``
        - AIMEC results columns
        - Areal indicator columns
        - Evaluation metrics (recall, precision, f1_score, tversky_score, optimisationVariable)
    """

    avalancheDir = 'data/' + avaName

    # Load variables that a varied
    varParList = cfgProb['PROBRUN']['varParList'].split('|')

    # Folder where ini files from simulations are
    inDir = pathlib.Path('../' + avalancheDir + '/Outputs/com8MoTPSA/configurationFiles')
    # read parameterSetDF
    paramSetDF = readParamSetDF(inDir, varParList)

    # dataframe from AIMEC
    df_aimec = pd.read_csv(
        '../' + avalancheDir + '/Outputs/ana3AIMEC/com8MoTPSA/Results_' + avaName + '_ppr_lim_1_w_600resAnalysisDF.csv')

    # merge aimec with parameterSet
    df_merged = pd.merge(paramSetDF, df_aimec, on='simName', how='inner')

    # Folder where pickle file is saved
    arealIndicatorDir = pathlib.Path('../' + avalancheDir + '/Outputs' + '/out1Peak' + '/arealIndicators.pkl')
    # read areal indicators
    indicatorsDF = readArealIndicators(arealIndicatorDir)

    # merge df
    df_merged = df_merged.merge(indicatorsDF, on="simName", how="left")
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
    paramLossDF = pd.DataFrame(df['parameterSet'].tolist(), columns=paramAll)
    paramLossDF = paramLossDF[paramSelected]
    paramLossDFScaled = (paramLossDF - paramLossDF.min()) / (paramLossDF.max() - paramLossDF.min())  # normalise
    paramLossDF['Loss'] = df['optimisationVariable']  # merging works with different index, but it is the right value
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
                     nu=2.5)
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


def optimiseNonSeqV1(pipe, paramSelected):

    # ToDo create LH samples with avaframe functions

    param_bounds = {
        'Dry-friction coefficient (-)': (0.15, 0.4),
        'Density (kg/m^3)': (50, 300),
        'Deposition rate 21  (m/s)': (0.15, 0.5),
        'Basal drag coeff. 1-2 (-)': (0.02, 0.05),
        'Top drag coeff. (-)': (1e-5, 1e-3),
        'Avalanche shear strength (Pa)': (0.0, 5.0),
        'Turbulent drag coefficient (-)': (0.0015, 0.005),
    }

    # 10.000 zufällige Punkte
    N = 10000
    lhs_unit = lhs(len(paramSelected), samples=N)  # neue LHS samples scaled

    X_candidates = np.array([
        lhs_unit[:, i] * (param_bounds[p][1] - param_bounds[p][0]) + param_bounds[p][0]
        for i, p in enumerate(paramSelected)
    ]).T

    df_candidates = pd.DataFrame(X_candidates, columns=paramSelected)

    # Vorhersage des Verlusts mit GP-Modell
    mu, sigma = pipe.predict(df_candidates, return_std=True)

    '''
    # Besten Punkt finden
    idx_best = np.argmin(mu)
    best_params = df_candidates.iloc[idx_best]
    best_loss = mu[idx_best]

    print("🔍 Beste Parameterkombination laut GP:")
    print(best_params)
    print(f"📉 Erwarteter Verlust: {best_loss:.4f}")
    '''
    # Top 100 Punkte mit kleinstem mu
    idx_top = np.argsort(mu)[:100]

    top_params = df_candidates.iloc[idx_top].copy()
    top_params["mu"] = mu[idx_top]
    top_params["sigma"] = sigma[idx_top]

    # Statistiken
    mean_params = top_params[paramSelected].mean()
    std_params = top_params[paramSelected].std()
    mean_mu = top_params["mu"].mean()
    std_mu = top_params["mu"].std()
    mean_sigma = top_params["sigma"].mean()
    std_sigma = top_params["sigma"].std()

    print("🔍 Mittelwerte ± Std der besten 100 Parameterkombinationen:")

    for p in paramSelected:
        m = mean_params[p]
        s = std_params[p]
        print(f"  {p:30s}: {m:.4f} ± {s:.4f}")

    print(f"\n📉 mu:    {mean_mu:.4f} ± {std_mu:.4f}")
    print(f"📊 sigma: {mean_sigma:.4f} ± {std_sigma:.4f}")
























