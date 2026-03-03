"""
Load and validate mouse clustering log files.

The mouse log is a CSV mapping filenames to Mouse_IDs:
    Filename,Mouse_ID
    cell1.abf,Mouse1
    cell2.abf,Mouse1
    cell3.abf,Mouse2
"""

import os
import logging
import pandas as pd
from typing import Dict

logger = logging.getLogger(__name__)


def load_mouse_log(path: str) -> Dict[str, str]:
    """Load a mouse log CSV and return a filename -> Mouse_ID mapping.
    
    Parameters
    ----------
    path : str
        Path to the mouse log CSV. Must have columns 'Filename' and 'Mouse_ID'.
    
    Returns
    -------
    dict
        Mapping of filename (without path, e.g. 'cell1.abf') to Mouse_ID string.
    
    Raises
    ------
    FileNotFoundError
        If the file does not exist.
    ValueError
        If required columns are missing or the file is empty.
    """
    if not os.path.exists(path):
        raise FileNotFoundError(f"Mouse log not found: {path}")
    
    df = pd.read_csv(path)
    
    # Normalize column names (strip whitespace, case-insensitive match)
    df.columns = df.columns.str.strip()
    col_map = {c.lower(): c for c in df.columns}
    
    if 'filename' not in col_map:
        raise ValueError(f"Mouse log must have a 'Filename' column. Found: {list(df.columns)}")
    if 'mouse_id' not in col_map:
        raise ValueError(f"Mouse log must have a 'Mouse_ID' column. Found: {list(df.columns)}")
    
    filename_col = col_map['filename']
    mouse_id_col = col_map['mouse_id']
    
    # Drop rows with missing values
    df = df.dropna(subset=[filename_col, mouse_id_col])
    
    if df.empty:
        raise ValueError("Mouse log is empty after dropping missing values")
    
    # Build mapping (use basename only). Data always has .abf in filename;
    # if the user omitted it in the log, add that key so lookup succeeds.
    mapping = {}
    for _, row in df.iterrows():
        fname = os.path.basename(str(row[filename_col]).strip())
        mouse_id = str(row[mouse_id_col]).strip()
        if fname in mapping and mapping[fname] != mouse_id:
            logger.warning(f"Duplicate filename '{fname}' with different Mouse_IDs: "
                         f"'{mapping[fname]}' vs '{mouse_id}'. Using last value.")
        mapping[fname] = mouse_id
        if not fname.lower().endswith('.abf'):
            mapping[fname + '.abf'] = mouse_id

    logger.info(f"Loaded mouse log: {len(df)} files across "
                f"{len(set(mapping.values()))} mice")
    return mapping


def validate_mouse_log(mouse_log: Dict[str, str], filenames: list) -> list:
    """Check that all filenames in the data have a Mouse_ID entry.
    
    Parameters
    ----------
    mouse_log : dict
        Filename -> Mouse_ID mapping.
    filenames : list
        List of filenames from the loaded data.
    
    Returns
    -------
    list
        List of filenames missing from the mouse log (empty if all present).
    """
    missing = []
    for fname in filenames:
        basename = os.path.basename(str(fname).strip())
        if basename not in mouse_log:
            missing.append(basename)
    
    if missing:
        logger.warning(f"{len(missing)} files not found in mouse log: {missing[:5]}"
                      + ("..." if len(missing) > 5 else ""))
    return missing
