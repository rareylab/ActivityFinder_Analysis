"""
LaTeX Table Generation Utilities
=================================

Adds:
- Automatic numeric column S spec inference
- Forced textual columns for categorical fields
- Header alias & shortening support to reduce table width (user can pass a mapping)
- Rich caption override support (caller responsibility)

Use `build_latex_table(df, caption, label, header_aliases={'long_column_name':'short'})`
to shorten verbose headers while keeping the underlying CSV columns untouched.
"""

import re
import pandas as pd
from typing import Optional, Dict, Any, List
from loguru import logger


def _escape_latex(text: str) -> str:
    if pd.isna(text):
        return ""
    text = str(text)
    # If the string is already a math-mode wrapped token like \( ... \),
    # leave it untouched so that backslashes are not escaped into \textbackslash.
    if text.startswith(r"\(") and text.endswith(r"\)"):
        return text
    replacements = [
        ("\\", "\\textbackslash{}"),
        ("{", "\\{"),
        ("}", "\\}"),
        ("$", "\\$"),
        ("&", "\\&"),
        ("%", "\\%"),
        ("#", "\\#"),
        ("_", "\\_"),
        ("~", "\\textasciitilde{}"),
        ("^", "\\textasciicircum{}"),
    ]
    for old, new in replacements:
        text = text.replace(old, new)
    return text


def _analyze_numeric_column(series: pd.Series) -> Dict[str, Any]:
    valid_values = series.dropna()
    if len(valid_values) == 0:
        return {
            'integer_digits': 1,
            'decimal_digits': 0,
            'format_spec': 'S[table-format=1.0]'
        }
    all_integers = all(float(v).is_integer() for v in valid_values if not pd.isna(v))
    if all_integers:
        max_val = abs(valid_values).max()
        integer_digits = len(str(int(max_val)))
        return {
            'integer_digits': integer_digits,
            'decimal_digits': 0,
            'format_spec': f'S[table-format={integer_digits}.0]'
        }
    max_abs = abs(valid_values).max()
    integer_digits = len(str(int(max_abs))) if max_abs >= 1 else 1
    max_decimals = 0
    for v in valid_values:
        if not pd.isna(v):
            str_val = f"{float(v):.10f}".rstrip('0').rstrip('.')
            if '.' in str_val:
                decimals = len(str_val.split('.')[1])
                max_decimals = max(max_decimals, decimals)
    decimal_digits = min(max_decimals, 4)
    return {
        'integer_digits': integer_digits,
        'decimal_digits': decimal_digits,
        'format_spec': f'S[table-format={integer_digits}.{decimal_digits}]'
    }


def _detect_column_types(df: pd.DataFrame) -> Dict[str, str]:
    column_specs = {}
    force_text_substrings = (
        "identity_bucket",
        "chembl_type",
        "variant",
        "analysis_type",
        "category",
    )
    for col in df.columns:
        lname = col.lower()
        if any(s in lname for s in force_text_substrings):
            column_specs[col] = 'l'
            continue
        try:
            numeric_series = pd.to_numeric(df[col], errors='coerce')
            if numeric_series.notna().any():
                spec_info = _analyze_numeric_column(numeric_series)
                column_specs[col] = spec_info['format_spec']
            else:
                column_specs[col] = 'l'
        except Exception:
            column_specs[col] = 'l'
    return column_specs


def render_latex_table(
    df: pd.DataFrame,
    caption: str,
    label: str,
    index_name_behavior: str = 'preserve',
    auto_format: bool = True,
    numeric_column_specs: Optional[Dict[str, str]] = None,
    footnotesize: bool = True,
    landscape: bool = False,
    adjustbox: bool = False,
    header_aliases: Optional[Dict[str, str]] = None,
    squeeze: bool = True,
) -> str:
    df_copy = df.copy()
    if auto_format:
        column_specs = _detect_column_types(df_copy)
        if numeric_column_specs:
            column_specs.update(numeric_column_specs)
    else:
        column_specs = numeric_column_specs or {col: 'l' for col in df_copy.columns}
    if header_aliases is None:
        header_aliases = {}
    if squeeze:
        # Try to collapse wide float specs with trailing .0 to minimal width
        squeezed = {}
        for k, spec in column_specs.items():
            if spec.startswith("S[table-format=") and spec.endswith(".0]"):
                squeezed[k] = spec.replace(".0]", ".0]")
            else:
                squeezed[k] = spec
        column_specs = squeezed
    col_format = ' '.join(column_specs.values())
    for col in df_copy.columns:
        if column_specs.get(col, 'l') == 'l':
            df_copy[col] = df_copy[col].apply(lambda x: _escape_latex(str(x)) if not pd.isna(x) else '')
    lines: List[str] = []
    if landscape:
        lines.append("\\begin{landscape}")
    lines.append("\\begin{table}[!htbp]")
    lines.append("\\centering")
    if footnotesize:
        lines.append("\\footnotesize")
    lines.append(f"\\caption{{{caption}}}")
    lines.append(f"\\label{{{label}}}")
    if adjustbox:
        lines.append("\\adjustbox{width=\\textwidth}{%")
    lines.append(f"\\begin{{tabular}}{{{col_format}}}")
    lines.append("\\toprule")
    header_parts = []
    for col, spec in column_specs.items():
        raw_header = header_aliases.get(col, col)
        escaped = _escape_latex(raw_header)
        if spec.startswith('S'):
            header_parts.append(f"{{{escaped}}}")
        else:
            header_parts.append(escaped)
    lines.append(" & ".join(header_parts) + " \\\\")
    lines.append("\\midrule")
    lines.append("\\rowcolors{2}{gray!15}{white}")
    for _, row in df_copy.iterrows():
        row_parts = []
        for col in df_copy.columns:
            val = row[col]
            if pd.isna(val):
                row_parts.append("")
            elif column_specs.get(col, 'l').startswith('S'):
                row_parts.append(str(val))
            else:
                row_parts.append(str(val))
        lines.append(" & ".join(row_parts) + " \\\\")
    lines.append("\\bottomrule")
    lines.append("\\end{tabular}")
    if adjustbox:
        lines.append("}")
    lines.append("\\end{table}")
    if landscape:
        lines.append("\\end{landscape}")
    return "\n".join(lines)


def build_latex_table(
    df: pd.DataFrame,
    caption: str,
    label: str,
    numeric_column_specs: Optional[Dict[str, str]] = None,
    auto_format: bool = True,
    header_aliases: Optional[Dict[str, str]] = None,
    squeeze: bool = True,
    **kwargs
) -> str:
    return render_latex_table(
        df=df,
        caption=caption,
        label=label,
        numeric_column_specs=numeric_column_specs,
        auto_format=auto_format,
        header_aliases=header_aliases,
        squeeze=squeeze,
        **kwargs
    )
