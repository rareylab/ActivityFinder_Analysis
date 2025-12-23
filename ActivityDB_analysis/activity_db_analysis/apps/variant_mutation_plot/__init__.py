"""Variant Mutation Plot analysis module for ChEMBL database analysis.

This module provides functionality to analyze and visualize the distribution
of mutation counts in ChEMBL variant sequences.
"""

from .variant_mutation_analysis import (
    execute_mutation_analysis,
    plot_mutation_distribution,
    generate_latex_table
)

__all__ = [
    'execute_mutation_analysis',
    'plot_mutation_distribution',
    'generate_latex_table'
]
