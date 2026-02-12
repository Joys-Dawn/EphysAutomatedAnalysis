# -*- mode: python ; coding: utf-8 -*-

import os
import importlib

block_cipher = None

# Collect submodules that PyInstaller often misses
hiddenimports = [
    # pingouin and its dependencies
    'pingouin',
    'pingouin.bayesian',
    'pingouin.circular',
    'pingouin.correlation',
    'pingouin.distribution',
    'pingouin.effsize',
    'pingouin.multicomp',
    'pingouin.nonparametric',
    'pingouin.parametric',
    'pingouin.pairwise',
    'pingouin.plotting',
    'pingouin.power',
    'pingouin.regression',
    'pingouin.reliability',
    # scikit-posthocs
    'scikit_posthocs',
    # statsmodels submodules
    'statsmodels',
    'statsmodels.stats',
    'statsmodels.stats.multitest',
    'statsmodels.stats.anova',
    'statsmodels.stats.oneway',
    'statsmodels.formula.api',
    'statsmodels.regression.mixed_linear_model',
    'statsmodels.tsa',
    'statsmodels.tsa.statespace',
    # pyabf
    'pyabf',
    # python-pptx
    'pptx',
    'pptx.util',
    'pptx.enum.text',
    'pptx.oxml.xmlchemy',
    # scipy submodules
    'scipy.optimize',
    'scipy.signal',
    'scipy.ndimage',
    'scipy.stats',
    # patsy
    'patsy',
    # seaborn
    'seaborn',
    # matplotlib backends
    'matplotlib',
    'matplotlib.backends.backend_tkagg',
    # pandas
    'pandas',
    # numpy
    'numpy',
    # application packages
    'analysis_code',
    'analysis_code.analyze_abf',
    'analysis_code.brief_current',
    'analysis_code.current_steps',
    'analysis_code.IC_gap_free',
    'analysis_code.vc_test',
    'analysis_code.trace_analysis',
    'analysis_code.subthresh_features',
    'analysis_code.time_series_utils',
    'modular_analysis',
    'modular_analysis.data_extraction',
    'modular_analysis.data_extraction.extractor',
    'modular_analysis.statistical_analysis',
    'modular_analysis.statistical_analysis.analyzer',
    'modular_analysis.statistical_analysis.designs',
    'modular_analysis.statistical_analysis.plotting',
    'modular_analysis.statistical_analysis.tests',
    'modular_analysis.statistical_analysis.tests.unpaired_ttest',
    'modular_analysis.statistical_analysis.tests.paired_ttest',
    'modular_analysis.statistical_analysis.tests.oneway_anova',
    'modular_analysis.statistical_analysis.tests.two_way_anova',
    'modular_analysis.statistical_analysis.tests.repeated_measures_anova',
    'modular_analysis.statistical_analysis.tests.mixed_anova',
    'modular_analysis.statistical_analysis.tests.frequency_analysis',
    'modular_analysis.statistical_analysis.tests.frequency_analysis_dependent',
    'modular_analysis.statistical_analysis.tests.attenuation_analysis',
    'modular_analysis.statistical_analysis.tests.attenuation_analysis_dependent',
    'modular_analysis.statistical_analysis.tests.posthoc_utils',
    'modular_analysis.shared',
    'modular_analysis.shared.config',
    'modular_analysis.shared.data_models',
    'modular_analysis.shared.utils',
]

a = Analysis(
    ['modular_analysis_app.py'],
    pathex=[],
    binaries=[],
    datas=[],
    hiddenimports=hiddenimports,
    hookspath=[],
    hooksconfig={},
    runtime_hooks=[],
    excludes=[],
    win_no_prefer_redirects=False,
    win_private_assemblies=False,
    cipher=block_cipher,
    noarchive=False,
)

pyz = PYZ(a.pure, a.zipped_data, cipher=block_cipher)

exe = EXE(
    pyz,
    a.scripts,
    [],
    exclude_binaries=True,
    name='ModularAnalysis',
    debug=False,
    bootloader_ignore_signals=False,
    strip=False,
    upx=True,
    console=False,  # No console window -- tkinter app
)

coll = COLLECT(
    exe,
    a.binaries,
    a.zipfiles,
    a.datas,
    strip=False,
    upx=True,
    upx_exclude=[],
    name='ModularAnalysis',
)
