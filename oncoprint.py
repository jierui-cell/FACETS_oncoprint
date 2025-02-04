# Author: Jierui Xu
# Date: 2025-02-03

import matplotlib.pyplot as plt
import matplotlib.patches as patches
import matplotlib.collections as mc
import pandas as pd
import numpy as np

from matplotlib.patches import Rectangle
from matplotlib.colors import ListedColormap
from matplotlib.collections import PatchCollection
import matplotlib.patches as mpatches

def plot_oncoprint_arm_level_cna_on_axis_patches(
    df_arm_level, 
    samples, 
    tcn_lcn_targets, 
    ax,
    color_all_cn_states=False, 
    return_order=False, 
    use_order=None,
    add_line=None,
    sort_arm='17p'  # <--- New argument with default value '17p'
):
    """
    Plots the CNA matrix using patches, assigning a DISTINCT color 
    to each (tcn, lcn) combination. A legend at the bottom shows
    the color <-> (tcn, lcn) mapping.

    Returns the final sorted sample order if return_order=True.
    """

    # Create a color dictionary for each (tcn, lcn). 
    if not color_all_cn_states: 
        color_palette = ['#ca0020','#ededed','#92c5de','#0571b0'][::-1] 
        tcn_lcn_order = [(1,0), (2,0), (2,1), (3,1)]
    else: 
        color_palette = ['#67001f','#b2182b','#d6604d','#f4a582','#fddbc7','#d1e5f0','#92c5de','#ededed','#2166ac','#053061'][::-1]
        tcn_lcn_order = [(1, 0), (2, 0), (2, 1), (3, 0), (3, 1), (4, 1), (4, 2), (5, 1), (5, 2), (6, 2)]  

    color_for_tuple = {}
    for idx, (tcn_val, lcn_val) in enumerate(tcn_lcn_order):
        color_for_tuple[(tcn_val, lcn_val)] = color_palette[idx % len(color_palette)]

    unique_arms =  [f"{arm}{suffix}" for arm in range(1,23) for suffix in ['p','q']]
    # Rows = arms, Columns = samples (the input parameter)
    df_combo = pd.DataFrame(
        data=0,
        index=unique_arms,  # each arm is a row
        columns=samples      # each sample is a column
    )

    # 2c) For each (tcn_one, lcn_one) target, assign an integer ID = idx
    for idx, (tcn_one, lcn_one) in enumerate(tcn_lcn_targets, start=1):
        # Filter rows that match this (tcn, lcn) AND have frac_of_arm >= 0.5
        mask = (
            (df_arm_level['tcn'] == tcn_one) &
            (df_arm_level['lcn'] == lcn_one) &
            (df_arm_level['frac_of_arm'] >= 0.5)
        )
        subset = df_arm_level[mask]

        # For each row in that subset, fill df_combo with idx
        # Only fill if the cell is currently 0 (i.e., not already assigned)
        for row in subset.itertuples(index=False):
            # row => (sample, arm, tcn, lcn, cn_length, arm_length, frac_of_arm, cn_state, ...)
            if row.sample in df_combo.columns and row.arm in df_combo.index:
                if df_combo.at[row.arm, row.sample] == 0:
                    df_combo.at[row.arm, row.sample] = idx

    if df_combo is None:
        return

    # Remove acrocentric arms if present
    acrocentric_arms = ['13p', '14p', '15p', '21p', '22p']
    arms_to_drop = [arm for arm in acrocentric_arms if arm in df_combo.columns]
    df_combo.drop(columns=arms_to_drop, errors='ignore', inplace=True)

    # -------------------------------------------------------------------
    # 3) SORT COLUMNS (NEW oncoprint-style LOGIC FROM plot_oncoprint_arm_level_cna)
    # -------------------------------------------------------------------
    final_sorted_columns = []

    # 3a) If the specified arm is in the index, sort by that first (descending)
    if sort_arm in df_combo.index:
        sorted_columns = df_combo.loc[:, :].sort_values(by=sort_arm, axis=1, ascending=False).columns.tolist()
        final_sorted_columns += [col for col in sorted_columns if df_combo.loc[sort_arm, col] != 0]
    else:
        # If '17p' doesn't exist, just use the full list
        sorted_columns = df_combo.columns.tolist()

    # 3b) For each arm except the specified sort_arm, sort by that arm in descending order
    df_combo_index_without_sort_arm = [arm for arm in df_combo.index if arm != sort_arm]
    for arm in df_combo_index_without_sort_arm:
        remaining_columns = [col for col in sorted_columns if col not in final_sorted_columns]
        # sort descending by that arm
        sorted_columns_arm = df_combo.loc[arm, remaining_columns].sort_values(ascending=False).index.tolist()
        # add columns with arm != 0
        final_sorted_columns += [col for col in sorted_columns_arm if df_combo.loc[arm, col] != 0]

    # 3c) Finally, add columns that have 0 for all arms
    #     i.e., columns not yet included => no CNA in any arm
    no_change_columns = [col for col in sorted_columns 
                        if col not in final_sorted_columns 
                        and df_combo.loc[df_combo_index_without_sort_arm, col].sum() == 0]
    final_sorted_columns += no_change_columns

    # If user-supplied order is provided, respect that instead
    if use_order:
        final_sorted_columns = use_order

    # Reindex columns in the new order
    df_combo = df_combo[final_sorted_columns]

    # Move the specified arm to the top row if present
    if sort_arm in df_combo.index:
        new_index = [sort_arm] + [arm for arm in df_combo.index if arm != sort_arm]
        df_combo = df_combo.reindex(new_index)
    # -------------------------------------------------------------------

    # -----------------------------
    # 4) Compute % of non-zero in each row
    # -----------------------------
    alteration_percentages = (df_combo != 0).mean(axis=1) * 100

    # -----------------------------
    # 5) Plot
    # -----------------------------
    arms = df_combo.index.tolist()
    samples = df_combo.columns.tolist()
    n_rows, n_cols = len(arms), len(samples)

    patches_list = []
    facecolors = []
    for i, arm in enumerate(arms):
        for j, sample in enumerate(samples):
            val = df_combo.at[arm, sample]
            if val == 0:
                color = 'white'
            else:
                combo_idx = val - 1  # index into tcn_lcn_targets
                combo_key = tcn_lcn_targets[combo_idx]
                color = color_for_tuple[combo_key]
            rect = Rectangle((j, i), 1, 1)
            patches_list.append(rect)
            facecolors.append(color)

    pc = PatchCollection(
        patches_list, 
        facecolor=facecolors,
        linewidths=0.5
    )
    ax.add_collection(pc)

    if add_line is not None:
        ax.axvline(x=add_line, color='red', linestyle='-', linewidth=1)

    ax.set_xlim(0, n_cols)
    ax.set_ylim(0, n_rows)
    ax.invert_yaxis()
    ax.set_yticks(np.arange(n_rows) + 0.5)
    ax.set_yticklabels(arms)

    # Put x-tick every 250 samples
    tick_positions = np.arange(0, n_cols, 250)
    ax.set_xticks(tick_positions)
    ax.set_xticklabels(tick_positions)

    # Add alteration percentages
    for i, (arm, pct) in enumerate(zip(arms, alteration_percentages)):
        ax.text(-len(df_combo.columns) * 0.03, i+0.5, f"{int(pct)}%", ha="right", va="center")

    ax.set_aspect('auto')

    # -----------------------------
    # 6) Add legend at the bottom
    # -----------------------------
    legend_handles = []
    for (tcn_val, lcn_val) in tcn_lcn_targets:
        c = color_for_tuple[(tcn_val, lcn_val)]
        label_txt = f"({tcn_val}, {lcn_val})"
        legend_patch = mpatches.Patch(facecolor=c, edgecolor='white', label=label_txt)
        legend_handles.append(legend_patch)

    ax.legend(
        handles=legend_handles,
        loc='upper center',
        bbox_to_anchor=(0.5, -0.05),  
        fancybox=True,
        shadow=False,
        ncol=len(tcn_lcn_targets)  
    )

    if return_order:
        return df_combo.columns.tolist()
    

# --------------------------------------------
# Example color map, adapted from your variables
# --------------------------------------------
COLOR_MAP = {
    # Mutation colors
    'Truncating Mutation (putative driver)': '#000000',
    'Truncating Mutation (putative passenger)': '#708090',
    'Missense Mutation (putative driver)': '#008000',
    'Missense Mutation (putative passenger)': '#53D400',
    'Inframe Mutation (putative driver)': '#993404',
    'Inframe Mutation (putative passenger)': '#a68028',
    'Splice Mutation (putative driver)': '#e5802b',
    'Splice Mutation (putative passenger)': '#f0b87b',
    'promoter': '#00B7CE',
    'promoter (putative passenger)': '#8cedf9',
    'other': '#cf58bc',
    'other (putative passenger)': '#f96ae3',
    # Structural variant colors
    'Structural Variant (putative driver)': '#8B00C9',
    'Structural Variant (putative passenger)': '#ce92e8',
    # CNA colors
    'Amplification': '#ff0000',
    'Deep Deletion': '#0000ff',
    # fallback if not found
    'default': '#cccccc'
}

# --------------------------------------------
# Helper functions to classify events
# --------------------------------------------
def is_mutation(event_name: str) -> bool:
    return any(word in event_name.lower() for word in [
        'truncating', 'missense', 'inframe', 'splice', 'promoter', 'other'
    ])

def is_structural_variant(event_name: str) -> bool:
    return 'structural' in event_name.lower() 

def is_cna(event_name: str) -> bool:
    return event_name in ['Amplification', 'Deep Deletion', 'Gain', 'Het Loss']

def get_event_color(event_name: str) -> str:
    """Return color from the COLOR_MAP, or default if not found."""
    return COLOR_MAP.get(event_name, COLOR_MAP['default'])


# --------------------------------------------
# Optimized Oncoprint Using PatchCollection
# --------------------------------------------
def plot_oncoprint_fast(df: pd.DataFrame, 
                        genes: list,
                        sample_col: str = 'sample', 
                        track_col: str = 'track', 
                        event_col: str = 'event',
                        row_gap=0.1, 
                        sample_order=None,
                        ax=None,
                        figsize=(20,6)):

    # If no ax is provided, create a new figure
    if ax is None:
        fig, ax = plt.subplots(figsize=figsize)
    else:
        fig = ax.figure

    # ----------------------------
    # 1) Filter and pivot to include all samples and specified genes
    # ----------------------------
    df_filtered = df[df[track_col].isin(genes)]
    
    # Ensure all samples are represented, filling missing entries with empty lists
    if sample_order is None:
        print ('No Sample order provided, using the order in the data')
        sample_order = df_filtered[sample_col].unique()

    # 1) Pivot to get df_pivot
    df_pivot = df_filtered.pivot_table(
        index=track_col,
        columns=sample_col,
        values=event_col,
        aggfunc=lambda x: x
    )

    # 2) Reindex df_pivot columns to match sample_order
    df_pivot = df_pivot.reindex(columns=sample_order, fill_value='')

    # 3) Melt and pivot again
    df_list = (
        df_pivot
        .stack()
        .reset_index(name=event_col)
        .pivot(index=track_col, columns=sample_col, values=event_col)
    )

    # 4) Reindex rows & columns to enforce order
    df_list = (
        df_list
          .reindex(index=genes, fill_value='')          # Row order (genes)
          .reindex(columns=sample_order, fill_value='') # Column order
          .fillna('')                                   # Replace remaining NaNs
    )

    # 5) Replace NaN with empty strings
    df_list = df_list.fillna('')
        
    tracks = df_list.index.tolist()
    samples = df_list.columns.tolist()
    n_rows = len(tracks)
    n_cols = len(samples)

    # Prepare lists of patches (and colors) for each layer
    background_patches = []
    background_colors = []

    cna_patches = []
    cna_colors = []

    sv_patches = []
    sv_colors = []

    mut_patches = []
    mut_colors = []
            
    # 2) Collect patches for each (track, sample) cell
    for i, track in enumerate(tracks):
        # Instead of each row being at y = i,
        # we space them out by row_gap. For example, i*(1 + row_gap).
        row_y = i * (1 + row_gap)

        for j, sample in enumerate(samples):
            # (a) Coordinates for the cell
            x = j
            y = row_y
            
            # (b) Gather events
            cell_events = df_list.loc[track, sample]

            if isinstance(cell_events, str):
                cell_events = [cell_events]  # Wrap single string in a list

            elif isinstance(cell_events, np.ndarray):
                cell_events = cell_events.tolist()  # Convert array to list

            # Ensure cell_events is a list
            if not isinstance(cell_events, list):
                cell_events = [cell_events]

            # Filter events
            cna_events = [ev for ev in cell_events if ev and is_cna(ev)]
            sv_events  = [ev for ev in cell_events if ev and is_structural_variant(ev)]
            mut_events = [ev for ev in cell_events if ev and is_mutation(ev)]

            # Use only the first event
            cna_ev = cna_events[0] if cna_events else None
            sv_ev  = sv_events[0] if sv_events else None
            mut_ev = mut_events[0] if mut_events else None

            # 2.1 Background patch (light-grey)
            rect_bg = patches.Rectangle((x, y), 1, 1)
            background_patches.append(rect_bg)
            background_colors.append('lightgrey')

            # 2.2 CNA patch (full cell)
            if cna_ev is not None:
                cna_color = get_event_color(cna_ev)
                rect_cna = patches.Rectangle((x, y), 1, 1)
                cna_patches.append(rect_cna)
                cna_colors.append(cna_color)

            # 2.3 SV patch (2/3 rectangle in the vertical direction)
            if sv_ev is not None:
                sv_color = get_event_color(sv_ev)
                rect_x = x
                rect_y = y + 0.15
                rect_w = 1
                rect_h = 0.7
                rect_sv = patches.Rectangle((rect_x, rect_y), rect_w, rect_h)
                sv_patches.append(rect_sv)
                sv_colors.append(sv_color)

            # 2.4 Mutation patch (0.4 high band, full width)
            if mut_ev is not None:
                mut_color = get_event_color(mut_ev)
                band_height = 0.4
                rect_x = x
                rect_y = y + (1 - band_height) / 2
                rect_w = 1
                rect_h = band_height
                rect_mut = patches.Rectangle((rect_x, rect_y), rect_w, rect_h)
                mut_patches.append(rect_mut)
                mut_colors.append(mut_color)
    
    # 3) Create PatchCollections
    #    Provide edgecolor='white' so we can see the grid lines
    #    and line up to see black shapes distinctly.
    pc_bg = mc.PatchCollection(background_patches, facecolor=background_colors, linewidths=0.5, zorder=0)
    pc_cna = mc.PatchCollection(cna_patches, facecolor=cna_colors, linewidths=0.5, zorder=1)
    pc_sv = mc.PatchCollection(sv_patches, facecolor=sv_colors, linewidths=0.5, zorder=2)
    pc_mut = mc.PatchCollection(mut_patches, facecolor=mut_colors, linewidths=0.5, zorder=3)

    # 4) Add them to the axis in the correct order
    ax.add_collection(pc_bg)
    ax.add_collection(pc_cna)
    ax.add_collection(pc_sv)
    ax.add_collection(pc_mut)

    # 5) Formatting
    ax.set_xlim(0, n_cols)
    # The total height is n_rows * (1 + row_gap). We'll make y-limit that.
    total_height = n_rows * (1 + row_gap)
    ax.set_ylim(0, total_height)

    ax.invert_yaxis()

    # Because we've spaced rows, the "center" of row i is i*(1+row_gap)+0.5
    # We'll set yticks accordingly:
    row_centers = [i*(1+row_gap) + 0.5 for i in range(n_rows)]
    ax.set_yticks(row_centers)
    ax.set_yticklabels(tracks)

    tick_positions = np.arange(0, n_cols, 250) 
    ax.set_xticks(tick_positions)               
    ax.set_xticklabels(tick_positions)    

    # Hide spines
    for spine in ax.spines.values():
        spine.set_visible(False)
    
    ax.set_aspect('auto')
    plt.tight_layout()

    # Add alteration percentages
    alteration_percentages = []
    for track in tracks:
        # count how many samples are non-empty in df_list.loc[track, sample]
        # i.e., at least 1 event
        count_nonempty = sum(len(df_list.loc[track, s]) > 0 for s in samples)
        pct = 100.0 * count_nonempty / n_cols
        alteration_percentages.append(pct)

    # Place text to the left of each row, similar to how you do in the CNA code
    for i, (track, pct) in enumerate(zip(tracks, alteration_percentages)):
        # row_y center is row_centers[i]
        y_text = row_centers[i]
        # let's offset x by -0.5 or so
        x_text = -len(samples) * 0.04
        ax.text(x_text, y_text, f"{int(pct)}%", ha="right", va="center")

    # --------------------------------------------
    # Add legend at the bottom
    # --------------------------------------------
    # We'll build custom patches that visually mimic each shape in your oncoprint:
    import matplotlib.patches as mpatches

    # 1. Background (lightgrey), full cell
    legend_bg = mpatches.Rectangle((0,0), 1, 1, facecolor='lightgrey', edgecolor='white',
                                label='No Alteration')

    # 2. CNA (full cell, color example = red)
    #    We'll use a full square to represent a CNA. For the legend color, pick one typical CNA color
    legend_cna_amplification = mpatches.Rectangle((0,0), 1, 1, facecolor=COLOR_MAP['Amplification'], edgecolor='white',
                                label='Amplification') 
    legend_cna_deletion = mpatches.Rectangle((0,0), 1, 1, facecolor=COLOR_MAP['Deep Deletion'], edgecolor='white', 
                                label='Deep Deletion')

    # 3. Structural Variant (2/3 rectangle)
    #    We'll use a rectangle 0.7 high to mimic your main plot. Color example = purple.
    legend_sv = mpatches.Rectangle((0,0.15), 1, 0.7, facecolor=COLOR_MAP['Structural Variant (putative driver)'], edgecolor='white',
                                label='Structural Variant (putative driver)') 

    # 4. Mutation (0.4 high band)
    #    We'll use a rectangle 0.4 high. Color example = black for a truncating mutation.
    legend_mut_truncating = mpatches.Rectangle((0,0.3), 1, 0.4, facecolor=COLOR_MAP['Truncating Mutation (putative driver)'], edgecolor='white',
                                label='Truncating Mutation (putative driver)') 
    legend_mut_missense = mpatches.Rectangle((0,0.3), 1, 0.4, facecolor=COLOR_MAP['Missense Mutation (putative driver)'], edgecolor='white', 
                                label='Missense Mutation (putative driver)') 
    legend_mut_inframe = mpatches.Rectangle((0,0.3), 1, 0.4, facecolor=COLOR_MAP['Inframe Mutation (putative driver)'], edgecolor='white',
                                label='Inframe Mutation (putative driver)')
    legend_mut_splice = mpatches.Rectangle((0,0.3), 1, 0.4, facecolor=COLOR_MAP['Splice Mutation (putative driver)'], edgecolor='white',
                                label='Splice Mutation (putative driver)')    

    # Now add them all to a legend. 
    # We'll place it below the plot using bbox_to_anchor, for example.
    legend_handles = [legend_bg, legend_cna_amplification, legend_cna_deletion, legend_sv, 
                      legend_mut_truncating, legend_mut_missense, legend_mut_inframe, legend_mut_splice] 
    ax.legend(handles=legend_handles,
            loc='upper center',        # place legend at the bottom by adjusting bbox_to_anchor
            bbox_to_anchor=(0.5, -0.05),
            fancybox=True,
            shadow=False,
            ncol=5)  
    
    # plt.show()


# now, plot the clustering arm level cna first, then plot the oncoprint with the same order 
def plot_cna_and_oncoprint_combined(
    df_arm_level,
    df_oncoprint,
    samples, 
    genes, 
    tcn_lcn_targets,
    color_all_cn_states=False, 
    use_order=None, 
    add_line=None, 
    sort_arm='17p'
):
    """
    1) Plots the CNA heatmap (arm-level) on the top axis using `plot_clustering_arm_level_cna_on_axis_patches`.
    2) Retrieves the sample order from that plot.
    3) Plots the oncoprint on the bottom axis with the same sample order.
    """

    # -- Dynamically compute the figure height for the bottom plot --
    # Each gene => about 0.3 height units, then round
    bottom_height = round(0.3 * len(genes))  
    # We'll add a fixed height (e.g., 5) for the top plot
    top_height = 6
    fig_height = top_height + bottom_height

    # Create the figure and two subplots
    fig, (ax_top, ax_bottom) = plt.subplots(
        nrows=2,
        figsize=(20, fig_height),
        gridspec_kw={"height_ratios": [top_height, bottom_height]}
)

    # 1) Plot CNA on the top axis; get the final sorted sample order
    sorted_cols = plot_oncoprint_arm_level_cna_on_axis_patches(
        df_arm_level=df_arm_level,
        samples=samples,
        tcn_lcn_targets=tcn_lcn_targets,
        color_all_cn_states=color_all_cn_states,
        ax=ax_top,
        return_order=True,  # so we get the columns
        use_order=use_order,
        add_line=add_line, 
        sort_arm=sort_arm
    )

    # 2) Plot the oncoprint on the bottom axis with the final sorted sample order
    #    Note: We do NOT pass figsize here, since we're already using subplots()
    plot_oncoprint_fast(
        df=df_oncoprint,
        genes=genes, 
        sample_col='sample',
        track_col='track',
        event_col='event',
        row_gap=0.1,
        sample_order=sorted_cols,  # align columns
        ax=ax_bottom  # draw on bottom axis
    )

    plt.tight_layout()
    plt.show()