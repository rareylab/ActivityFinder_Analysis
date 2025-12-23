"""
SQL Artifacts Management Package
===============================

This package provides comprehensive management tools for SQL artifacts used in
ligand mapping analysis, including functions, materialized views, views, and indexes.

The main module `sql_artifacts_management.py` contains a complete CLI application
built with Typer that allows fine-grained control over database objects.

Usage:
    python -m activity_db_analysis.apps.sql_artifacts.sql_artifacts_management --help
"""

from .sql_artifacts_management import app, ARTIFACTS, ArtifactType, DependencyManager, DatabaseExecutor

__all__ = [
    'app',
    'ARTIFACTS',
    'ArtifactType',
    'DependencyManager',
    'DatabaseExecutor'
]
