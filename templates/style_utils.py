# Figure style helpers shared by all templates; keep this file next to fair_utils.py
# fair_utils decides WHERE a figure is saved; this file decides WHAT it looks like.
# Three things are kept stable across every figure in a project:
#   1. colour binding  - one cell type keeps one colour in every figure (UMAP, heatmap band, dot plot)
#   2. colour scale    - diverging blue-white-red, centred on the semantic zero, fixed limits
#   3. verification    - no overlapping text is allowed to reach a saved file
# Compatible with Python 3.7

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

# Diverging map for z-scored expression. Centre (#f7f7f7) marks z = 0, not the data midpoint.
DIVERGING = mpl.colors.LinearSegmentedColormap.from_list(
    "style_utils_div", ["#2166ac", "#f7f7f7", "#b2182b"], N=256
)

# Default colour binding. Override per project with resolve_palette(..., preferred=...).
# Any label not listed here is served from EXTENSION_PALETTE instead.
DEFAULT_BINDING = {
    "Enterocyte": "#1f77b4", "Stem": "#ff7f0e", "TA": "#279e68",
    "T": "#d62728", "B": "#aa40fc", "Myeloid": "#8c564b",
    "Fibroblast": "#e377c2", "Pericyte/SMC": "#b5bd61", "Other stromal": "#17becf",
    "Endothelial": "#aec7e8", "Goblet": "#ffbb78", "EEC": "#98df8a", "Paneth": "#ff9896",
}

# Reserved hues for labels absent from the binding. Chosen not to collide with the values above.
EXTENSION_PALETTE = [
    "#8c6d31", "#d6616b", "#7b4173", "#843c39", "#5254a3", "#637939", "#a55194",
    "#6b6ecf", "#b5cf6b", "#ce6dbd", "#9c9ede", "#e7969c", "#de9ed6", "#ad494a",
]

# Common spellings of the same population, so a renamed dataset keeps its colours.
ALIASES = {
    "enterocytes": "Enterocyte", "absorptive": "Enterocyte",
    "stem cell": "Stem", "stem cells": "Stem", "isc": "Stem",
    "transit amplifying": "TA", "transit-amplifying": "TA",
    "t cell": "T", "t cells": "T", "t-cell": "T", "cd4 t": "T", "cd8 t": "T",
    "b cell": "B", "b cells": "B", "b-cell": "B", "plasma": "B", "plasma cell": "B",
    "macrophage": "Myeloid", "macrophages": "Myeloid", "monocyte": "Myeloid",
    "dc": "Myeloid", "dendritic": "Myeloid", "neutrophil": "Myeloid", "mono/mac": "Myeloid",
    "fibroblasts": "Fibroblast", "stromal": "Fibroblast", "fibro": "Fibroblast",
    "pericyte": "Pericyte/SMC", "smc": "Pericyte/SMC", "smooth muscle": "Pericyte/SMC",
    "endothelium": "Endothelial", "endothelial cell": "Endothelial", "ec": "Endothelial",
    "goblet cell": "Goblet", "enteroendocrine": "EEC", "paneth cell": "Paneth",
}


def resolve_palette(labels, preferred=None, extension=None, aliases=None):
    """Bind every label to a colour: exact match > alias match > reserved extension hue.

    Keeps colours stable when a dataset renames or adds cell types, so figures from
    different datasets stay visually comparable. Returns {label: hex}.
    """
    preferred = DEFAULT_BINDING if preferred is None else preferred
    extension = EXTENSION_PALETTE if extension is None else extension
    aliases = ALIASES if aliases is None else aliases

    taken = set()
    binding = {}
    for label in labels:
        key = str(label).strip()
        if key in preferred:
            binding[key] = preferred[key]
        elif key.lower() in aliases and aliases[key.lower()] in preferred:
            binding[key] = preferred[aliases[key.lower()]]
        else:
            binding[key] = None
        if binding[key] is not None:
            taken.add(binding[key])

    spare = [c for c in extension if c not in taken]
    cursor = 0
    for key in binding:
        if binding[key] is None:
            binding[key] = spare[cursor % len(spare)]
            cursor += 1
    return binding


def group_matrix(adata, genes, group_key, sub_key=None, layer="raw", min_cells=10):
    """Mean expression per group (optionally per group x sub-group), genes as columns.

    layer: "raw" reads adata.raw, None reads adata.X, otherwise adata.layers[layer].
    Combinations with fewer than min_cells cells are dropped, never averaged in.
    Returns (matrix, dropped) where dropped lists the excluded combinations and their n.
    """
    import scipy.sparse as sp

    if layer == "raw":
        source, var_names = adata.raw.X, adata.raw.var_names
    elif layer is None:
        source, var_names = adata.X, adata.var_names
    else:
        source, var_names = adata.layers[layer], adata.var_names

    missing = [g for g in genes if g not in var_names]
    if missing:
        raise KeyError("genes absent from the matrix: " + ", ".join(missing))

    columns = source[:, [var_names.get_loc(g) for g in genes]]
    if sp.issparse(columns):
        columns = columns.toarray()
    columns = np.asarray(columns)

    primary = adata.obs[group_key].astype(str).values
    secondary = adata.obs[sub_key].astype(str).values if sub_key else None

    rows, index, dropped = [], [], []
    levels = list(pd.unique(primary)) if secondary is None else [
        (a, b) for a in pd.unique(primary) for b in pd.unique(secondary)
    ]
    for level in levels:
        if secondary is None:
            mask = primary == level
        else:
            mask = (primary == level[0]) & (secondary == level[1])
        n = int(mask.sum())
        if n < min_cells:
            dropped.append((level, n))
            continue
        rows.append(columns[mask].mean(axis=0))
        index.append(level)

    matrix = pd.DataFrame(rows, index=index, columns=list(genes))
    if secondary is not None and len(matrix):
        matrix.index = pd.MultiIndex.from_tuples(matrix.index)
    return matrix, dropped


def zscore(matrix, limit=2.0):
    """Z-score each gene across the rows of the matrix, then clip to +/- limit."""
    scaled = (matrix - matrix.mean()) / matrix.std()
    return scaled.clip(-limit, limit)


def heatmap(matrix, group_order, sub_order=None, palette=None, limit=2.0,
            gene_families=None, figsize=None, title=None):
    """Genes as rows, group (x sub-group) as columns, with a colour band naming each group.

    matrix: rows = group or (group, sub-group), columns = genes, already z-scored.
    gene_families: {family name: [genes]} draws grouped row bands in the given order.
    Returns the Figure; pass it to check_overlap() and then to fair.save_fig().
    """
    if sub_order is None:
        columns, block = list(group_order), 1
    else:
        columns = [(g, s) for g in group_order for s in sub_order]
        block = len(sub_order)
    present = [c for c in columns if c in matrix.index]
    if len(present) < len(columns):
        raise ValueError("matrix is missing %d of the requested columns" % (len(columns) - len(present)))

    if gene_families:
        genes, bounds, cursor = [], [], 0
        for family, members in gene_families.items():
            members = [g for g in members if g in matrix.columns]
            bounds.append((family, cursor, cursor + len(members)))
            genes += members
            cursor += len(members)
    else:
        genes, bounds = list(matrix.columns), None

    values = matrix.loc[present, genes].values.T.astype(float)
    palette = resolve_palette(group_order) if palette is None else palette
    if figsize is None:
        figsize = (max(6.0, 0.24 * len(present) + 2.5), max(3.0, 0.22 * len(genes) + 1.8))

    fig, ax = plt.subplots(figsize=figsize)
    ax.imshow(values, aspect="auto", cmap=DIVERGING, vmin=-limit, vmax=limit, interpolation="nearest")

    if bounds:
        family_hues = resolve_palette([f for f, _, _ in bounds], preferred={})
        for family, start, stop in bounds:
            ax.axhspan(start - 0.5, stop - 0.5, color=family_hues[family], alpha=0.10, zorder=0)
        for _, _, stop in bounds[:-1]:
            ax.axhline(stop - 0.5, color="white", linewidth=1.8, zorder=2)

    ax.set_yticks(range(len(genes)))
    ax.set_yticklabels(genes, fontsize=8, fontstyle="italic")
    ax.tick_params(axis="y", length=0)

    if sub_order is None:
        ax.set_xticks(range(len(present)))
        ax.set_xticklabels([str(c) for c in present], fontsize=7, rotation=90)
    else:
        ax.set_xticks(range(len(present)))
        ax.set_xticklabels([str(c[1]) for c in present], fontsize=6, rotation=90)
        for edge in range(block, len(present), block):
            ax.axvline(edge - 0.5, color="black", linewidth=1.3)
    ax.tick_params(axis="x", length=0, pad=2)

    ax.set_xticks(np.arange(-0.5, len(present), 1), minor=True)
    ax.set_yticks(np.arange(-0.5, len(genes), 1), minor=True)
    ax.grid(which="minor", color="white", linewidth=0.7)
    ax.tick_params(which="minor", length=0)
    for spine in ax.spines.values():
        spine.set_visible(True)
        spine.set_color("black")
        spine.set_linewidth(1.0)

    bar = fig.colorbar(ax.images[0], ax=ax, fraction=0.02, pad=0.014, shrink=0.5)
    bar.set_label("Mean expression\n(z-score across columns)", fontsize=8)
    bar.ax.tick_params(labelsize=7)
    bar.outline.set_linewidth(0.6)

    band = ax.inset_axes([0, 1.006, 1, 0.045], transform=ax.transAxes)
    for i, group in enumerate(group_order):
        left = i * block
        band.axvspan(left - 0.5, left + block - 0.5, facecolor=palette[group],
                     edgecolor="white", linewidth=1.0)
        band.text(left + block / 2.0 - 0.5, 0.5, str(group), ha="center", va="center",
                  fontsize=7.8 if len(str(group)) <= 10 else 6.6,
                  fontweight="bold", color="white")
    band.set_xlim(-0.5, len(present) - 0.5)
    band.set_ylim(0, 1)
    band.axis("off")

    if title:
        ax.set_title(title, fontsize=10, fontweight="bold", pad=22)
    fig.tight_layout()
    return fig


def check_overlap(fig, scope="axes"):
    """Return the list of overlapping visible text pairs; [] means the figure is clean.

    scope="axes" compares tick labels within the same Axes, which is what dense
    multi-panel figures need: comparing across Axes reports false positives.
    scope="figure" compares every text object in the figure.
    """
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    found = []
    if scope == "figure":
        items = [(t, t.get_window_extent(renderer)) for t in fig.findobj(mpl.text.Text)
                 if t.get_text().strip() and t.get_visible()]
        groups = [items]
    else:
        groups = []
        for ax in fig.get_axes():
            labels = ax.get_xticklabels(which="both") + ax.get_yticklabels(which="both")
            groups.append([(t, t.get_window_extent(renderer)) for t in labels if t.get_text().strip()])
    for items in groups:
        for i, (first, box_a) in enumerate(items):
            for second, box_b in items[i + 1:]:
                if box_a.overlaps(box_b):
                    found.append((first.get_text(), second.get_text()))
    return found


def save_checked(dirs, name, fig, fair, scope="axes"):
    """Verify no text overlaps, then hand the figure to fair.save_fig (PDF + PNG).

    Raises before writing anything if the figure has overlapping labels, so a
    figure with unreadable text never reaches results/.
    """
    clashes = check_overlap(fig, scope=scope)
    if clashes:
        raise ValueError("overlapping labels, figure not saved: %s" % clashes[:5])
    fair.save_fig(dirs, name, fig=fig)
