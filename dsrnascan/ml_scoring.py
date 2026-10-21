"""
ML scoring module for dsRNAscan predictions.

Scores predicted dsRNA structures using pre-trained machine learning models:
  - Stability model: YDF RandomForest trained to predict RNA editing status
    from structural features alone (no annotation required).
  - Probing model: sklearn RandomForest trained to predict agreement with
    in vivo RNA structure probing data.

Models are bundled with the package and loaded on demand.
"""

import functools
import os
import logging
import numpy as np
import pandas as pd

logger = logging.getLogger(__name__)

MODELS_DIR = os.path.join(os.path.dirname(__file__), 'models')

# Youden thresholds (optimized on human editing data)
# These may not be directly applicable to non-human organisms
# Stability: pure structural RF (1000 trees, depth 5, no genic feature)
#   Test set (chr1, chr19, chrX): AUC=0.875, Sens=94.7%, Spec=72.3%
#   Full dataset: AUC=0.888, Sens=94.8%, Spec=74.9%
#   (0.2471 belonged to the published model, which also used a genic feature)
STABILITY_YOUDEN_THRESHOLD = 0.5526
# Probing: the model predicts a power-weighted advantage (about -0.07 to 0.16). It is reported on 0-1 by
# min-max scaling over the 5.13M human dsRNAs, as on the dsRNA browser; the Youden cutoff against
# editing, raw 0.0315, is 0.440 on that scale (the published resource's structure_probing_high_conf).
PROBING_RAW_MIN = -0.06649859532060612
PROBING_RAW_MAX = 0.1561966513188369
PROBING_YOUDEN_THRESHOLD_RAW = 0.0315
PROBING_YOUDEN_THRESHOLD = (PROBING_YOUDEN_THRESHOLD_RAW - PROBING_RAW_MIN) / (PROBING_RAW_MAX - PROBING_RAW_MIN)


def _check_ydf():
    """Check if ydf is available."""
    try:
        import ydf
        return True
    except ImportError:
        return False


def _check_sklearn():
    """Check if sklearn is available."""
    try:
        import joblib
        import sklearn
        return True
    except ImportError:
        return False


def check_ml_dependencies():
    """Check which ML scoring models are available based on installed packages.

    Returns dict with model availability and messages.
    """
    status = {
        'stability': _check_ydf(),
        'probing': _check_sklearn(),
    }
    missing = []
    if not status['stability']:
        missing.append('ydf (for stability model)')
    if not status['probing']:
        missing.append('scikit-learn (for probing model)')
    status['missing'] = missing
    return status


def derive_ml_features(df):
    """Derive ML model input features from dsRNAscan output columns.

    Parameters
    ----------
    df : pd.DataFrame
        dsRNAscan merged results with standard column names.

    Returns
    -------
    pd.DataFrame
        DataFrame with columns: bp_count, er_energy, percent_paired,
        loop_length, length.
    """
    # bp_count = longest contiguous helix
    bp_count = df['longest_helix'].astype(float)

    # er_energy = free energy of hybridization
    er_energy = df['dG(kcal/mol)'].astype(float)

    # percent_paired = percentage of paired bases
    percent_paired = df['percent_paired'].astype(float)

    # loop_length = distance between the two arms
    loop_length = (df['j_start'] - df['i_end']).abs().astype(float)

    # length = length of the longer arm (end - start), the training data's `length`: it reproduces the
    # published probing scores (99% exact on chr19; the shorter arm gives 84%)
    arm_i = (df['i_end'] - df['i_start']).abs().astype(float)
    arm_j = (df['j_end'] - df['j_start']).abs().astype(float)
    length = pd.concat([arm_i, arm_j], axis=1).max(axis=1)

    return pd.DataFrame({
        'bp_count': bp_count,
        'er_energy': er_energy,
        'percent_paired': percent_paired,
        'loop_length': loop_length,
        'length': length,
    })


def score_stability(features_df, threads=None):
    """Score predictions with the stability model (YDF RandomForest).

    Parameters
    ----------
    features_df : pd.DataFrame
        Output of derive_ml_features().

    Returns
    -------
    np.ndarray
        Probability scores (0-1) indicating likelihood of being an
        edited dsRNA based on structural properties.
    """
    import ydf

    model_path = os.path.join(MODELS_DIR, 'stability_rf_depth5')
    if not os.path.exists(model_path):
        raise FileNotFoundError(f"Stability model not found at {model_path}")

    model = _stability_model(model_path)
    try:
        scores = model.predict(features_df, num_threads=threads)
    except TypeError:  # older ydf
        scores = model.predict(features_df)
    return scores


@functools.lru_cache(maxsize=None)
def _stability_model(model_path):
    """Loaded once per run (results are scored batch by batch)."""
    import ydf
    previous = ydf.verbose(0)  # no "Loading model" lines on the console
    try:
        return ydf.load_model(model_path)
    finally:
        ydf.verbose(previous)


def score_probing(features_df, threads=None):
    """Score predictions with the probing model (sklearn RandomForest).

    Parameters
    ----------
    features_df : pd.DataFrame
        Output of derive_ml_features().

    Returns
    -------
    np.ndarray
        Predicted power-weighted advantage scores indicating likelihood
        of agreeing with in vivo RNA structure probing data.
    """
    import joblib

    model_dir = os.path.join(MODELS_DIR, 'probing_rf')
    model_path = os.path.join(model_dir, 'power_weighted_advantage_model.joblib')
    scaler_path = os.path.join(model_dir, 'power_weighted_advantage_scaler.joblib')

    if not os.path.exists(model_path):
        raise FileNotFoundError(f"Probing model not found at {model_path}")

    model, scaler = _probing_model(model_path, scaler_path)

    # Scaler expects features in order: length, percent_paired, loop_length, bp_count, er_energy
    probing_features = features_df[['length', 'percent_paired', 'loop_length', 'bp_count', 'er_energy']].copy()
    X_scaled = scaler.transform(probing_features)
    model.n_jobs = threads or -1  # the saved model uses every core on the machine
    scores = model.predict(X_scaled)
    return scores


@functools.lru_cache(maxsize=None)
def _probing_model(model_path, scaler_path):
    """Loaded once per run (results are scored batch by batch)."""
    import joblib
    import warnings
    with warnings.catch_warnings():
        warnings.filterwarnings('ignore', category=UserWarning, message='.*InconsistentVersionWarning.*')
        warnings.filterwarnings('ignore', message='.*Trying to unpickle.*')
        return joblib.load(model_path), joblib.load(scaler_path)


def score_results(results_df, models=None, threads=None, log=True):
    """Score dsRNAscan results with available ML models.

    Parameters
    ----------
    results_df : pd.DataFrame
        dsRNAscan merged results DataFrame.
    models : list of str, optional
        Which models to run. Default: all available.
        Options: 'stability', 'probing'.

    Returns
    -------
    pd.DataFrame
        Input DataFrame with appended score columns:
        - stability_model_score (if stability model available)
        - probing_model_score (if probing model available)
    """
    status = check_ml_dependencies()

    if models is None:
        models = []
        if status['stability']:
            models.append('stability')
        if status['probing']:
            models.append('probing')

    if not models:
        logger.warning("No ML scoring models available. Install ydf and/or scikit-learn.")
        return results_df

    features_df = derive_ml_features(results_df)
    n = len(results_df)

    if 'stability' in models:
        if not status['stability']:
            logger.warning("ydf not installed, skipping stability model")
        else:
            try:
                scores = score_stability(features_df, threads)
                results_df = results_df.copy()
                results_df['stability_model_score'] = np.round(scores, 3)
                results_df['likely_edited'] = np.where(scores >= STABILITY_YOUDEN_THRESHOLD, 'Yes', 'No')
                n_high = (results_df['likely_edited'] == 'Yes').sum()
                if log:
                    logger.info(f"  likely_edited (stability model): {n_high:,} of {n:,} ({n_high/n*100:.1f}%)")
            except Exception as e:
                logger.error(f"Stability model scoring failed: {e}")

    if 'probing' in models:
        if not status['probing']:
            logger.warning("scikit-learn not installed, skipping probing model")
        else:
            try:
                raw = score_probing(features_df, threads)
                scores = np.clip((raw - PROBING_RAW_MIN) / (PROBING_RAW_MAX - PROBING_RAW_MIN), 0, 1)
                if 'stability_model_score' not in results_df.columns:
                    results_df = results_df.copy()
                results_df['probing_model_score'] = np.round(scores, 3)
                results_df['likely_forms'] = np.where(raw >= PROBING_YOUDEN_THRESHOLD_RAW, 'Yes', 'No')
                n_high = (results_df['likely_forms'] == 'Yes').sum()
                if log:
                    logger.info(f"  likely_forms (probing model): {n_high:,} of {n:,} ({n_high/n*100:.1f}%)")
            except Exception as e:
                logger.error(f"Probing model scoring failed: {e}")

    return results_df


def score_results_file(input_path, output_path=None):
    """Score an existing dsRNAscan results file (*.dsrnas.txt).

    Parameters
    ----------
    input_path : str
        Path to a *.dsrnas.txt file.
    output_path : str, optional
        Path to write scored results. Default: overwrite input.

    Returns
    -------
    pd.DataFrame
        Scored results.
    """
    df = pd.read_csv(input_path, sep='\t')
    if len(df) == 0:
        logger.info("No predictions to score.")
        return df

    scored = score_results(df)

    if output_path is None:
        output_path = input_path
    scored.to_csv(output_path, sep='\t', index=False)
    return scored
