"""HET Code Plot analysis module for activity database analysis.

This module provides functionality to analyze and visualize the distribution
of variants per HET code in the Activity Finder database.
"""

from .het_code_variants import (
    execute_het_code_analysis,
    plot_het_code_variants,
    generate_latex_table
)

__all__ = [
    'execute_het_code_analysis',
    'plot_het_code_variants',
    'generate_latex_table'
]
