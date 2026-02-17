"""
pathSQE: Automated analysis of single-crystal inelastic neutron scattering data
"""

__version__ = "0.1.0"
__author__ = "Aiden Sable"

# Note: This package requires mantid to be installed.
# For full functionality, use: pixi install --environment mantid
# or: conda install -c mantid mantid

try:
    # Import all modules - they require mantid
    from . import core
    from . import helper
    from . import filter_functions
    from . import simulations
    from . import plotting_and_reports
    from . import Resolution
    from . import search_and_plot

    __all__ = [
        "core",
        "helper",
        "filter_functions",
        "simulations",
        "plotting_and_reports",
        "Resolution",
        "search_and_plot",
    ]
except ImportError as e:
    if "mantid" in str(e):
        raise ImportError(
            "The pathSQE package requires mantid to be installed. "
            "Please install mantid using one of the following:\n"
            "  - pixi install --environment mantid\n"
            "  - conda install -c mantid mantid\n"
            f"\nOriginal error: {e}"
        ) from e
    raise
