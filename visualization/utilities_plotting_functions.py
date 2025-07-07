def plot_ch_presence_absence_bars(muts_df, path_sample_information, color_dict, ax, annotate_what, lu_pts_n_total = None, caba_pts_n_total = None, add_p_values = False):
    """
    Plots the presence / absence of CH in the cohort as bar charts. Separates into two groups.
    """
    sample_info = pd.read_csv(path_sample_information, sep = "\t")    
    
    if lu_pts_n_total is None and caba_pts_n_total is None:
        lu_pts_n_total = sample_info[sample_info["Arm"] == "LuPSMA"]["Patient_id"].drop_duplicates().shape[0]
        caba_pts_n_total = sample_info[sample_info["Arm"] == "Cabazitaxel"]["Patient_id"].drop_duplicates().shape[0]
        
    if annotate_what.lower() == "chip":
        col_to_filter = "VAF_n"
    else:
        col_to_filter = "VAF_t"
    
    results_dict = {}
    for min_vaf in [0, 2, 10]:
        muts_df_filtered  = muts_df[muts_df[col_to_filter] >= min_vaf]
        
        lu_pts = muts_df_filtered[muts_df_filtered["Arm"] == "LuPSMA"]["Patient_id"].unique().shape[0]
        lu_perc = round((lu_pts/lu_pts_n_total)*100)
        
        caba_pts = muts_df_filtered[muts_df_filtered["Arm"] == "Cabazitaxel"]["Patient_id"].unique().shape[0]
        caba_perc = round((caba_pts/caba_pts_n_total)*100)
        
        # save results
        results_dict[min_vaf] = {"LuPSMA perc": lu_perc, "Cabazitaxel perc": caba_perc, "LuPSMA npts": lu_pts, "Cabazitaxel npts": caba_pts}
        # sum_perc = (lu_pts + caba_pts)/(group1_ntotal+group2_ntotal)
        # print(f"Min vaf={min_vaf}, {sum_perc}")
    
    # Now plotting the bar charts    
    df = pd.DataFrame.from_dict(results_dict, orient='index').reset_index().rename(columns = {"index": "min_vaf"})
    for i, row in df.iterrows():
        ax.bar(i-0.2, row["LuPSMA perc"], color = color_dict["LuPSMA"], edgecolor = None, width = 0.4)
        ax.bar(i+0.2, row["Cabazitaxel perc"], color = color_dict["Cabazitaxel"], edgecolor = None, width = 0.4)
        
        # Annotate the number of patients on top of the bars
        ax.text(i - 0.2, row["LuPSMA perc"]+1, str(row["LuPSMA npts"]), ha='center', va='bottom', fontsize=7, color='black')
        ax.text(i + 0.2, row["Cabazitaxel perc"]+1, str(row["Cabazitaxel npts"]), ha='center', va='bottom', fontsize=7, color='black')
    
    # Z test
    if add_p_values:
        p_values_ztest = []
        nobs = np.array([lu_pts_n_total, caba_pts_n_total])
        for i, (group1_count, group2_count) in enumerate(zip(df["LuPSMA npts"], df["Cabazitaxel npts"])):
            counts = np.array([group1_count, group2_count])
            _, p_val = proportions_ztest(counts, nobs, alternative="larger")  # alternative="larger" for one-sided test
            p_values_ztest.append(p_val)
        
        for i, p_val in enumerate(p_values_ztest):
            # Get the bladder count, the p values will be annotated right on top of it.
            count_value = df.loc[i, "LuPSMA perc"]
            if p_val < 0.001:
                ax.text(i + 0.2, count_value+3, "***", ha='center', va='bottom', color='black')
            elif p_val < 0.01:
                ax.text(i + 0.2, count_value+3, "**", ha='center', va='bottom', color='black')
            elif p_val < 0.05:
                ax.text(i + 0.2, count_value+3, "*", ha='center', va='bottom', color='black')
            
    ax.set_xlim((-0.7, 2.5))
    ax.set_xticks([0, 1, 2])
    ax.tick_params(axis='x', bottom=False)
    ax.tick_params(axis='x', pad=2)
    ax.set_xticklabels(["≥0.25%", "≥2%", "≥10%"], fontsize = 8)
    ax.set_ylabel("% of patients in arm")
    ax.set_xlabel("Minimum CH VAF%")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    
    ax.set_ylim((0, 100))
    ax.set_yticks([0, 25, 50, 75, 100])
    ax.set_yticklabels(["0", "25", "50", "75", "100"])
    
    # ADD LEGEND
    legend_colors = color_dict.values()
    legend_labels = color_dict.keys()
    legend_handles = [plt.Line2D([0], [0], marker='s', color=color, label=label, markersize=5, linestyle='') for color, label in zip(legend_colors, legend_labels)]
    ax.legend(handles=legend_handles, loc="upper right", frameon=False, fontsize = 8, handlelength=2, handletextpad = 0.1)
    return(ax)

def plot_ch_presence_absence_bars_not_grouped(muts_df, ax, ntotal_pts, annotate_what):
    """
    Plots the presence / absence of CH in the cohort as bar charts.
    DOESN'T SEPARATE INTO GROUPS. JUST ONE BLACK BAR.
    """
            
    if annotate_what.lower() == "chip":
        col_to_filter = "VAF_n"
    else:
        col_to_filter = "VAF_t"
    
    results_dict = {}
    for min_vaf in [0, 2, 10]:
        muts_df_filtered = muts_df[muts_df[col_to_filter] >= min_vaf].drop_duplicates("Patient_id")
        npts = muts_df_filtered.shape[0]
        perc = npts/ntotal_pts*100
        results_dict[min_vaf] = {"perc": perc, "npts": npts}
    
    # Now plotting the bar charts    
    df = pd.DataFrame.from_dict(results_dict, orient='index').reset_index().rename(columns = {"index": "min_vaf"})
    for i, row in df.iterrows():
        ax.bar(i, row["perc"], color = "goldenrod", edgecolor = None, width = 0.4)
        ax.text(i, row["perc"]+1, str(row["npts"]), ha='center', va='bottom', fontsize=7, color='black')
                
    ax.set_xlim((-0.7, 2.5))
    ax.set_xticks([0, 1, 2])
    ax.tick_params(axis='x', bottom=False)
    ax.tick_params(axis='x', pad=2)
    ax.set_xticklabels(["≥0.25%", "≥2%", "≥10%"], fontsize = 8)
    ax.set_ylabel("% of patients")
    ax.set_xlabel("Minimum CH VAF%")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    
    ax.set_ylim((0, 100))
    ax.set_yticks([0, 25, 50, 75, 100])
    ax.set_yticklabels(["0", "25", "50", "75", "100"])
    
    return(ax)


def plot_per_patient_counts_grouped_bar_chart(muts_df, grouping_col_name, color_dict, group1_ntotal, group2_ntotal, ax):
    """
    In a given ax plots the fraction of the bladder and kidney cohorts separately that have CH mutations.
    """
    # Extract names of groups
    group1_name = "LuPSMA"
    group2_name = "Cabazitaxel"
    
    # Determine colors for each group
    group1_color = color_dict[group1_name]
    group2_color = color_dict[group2_name]
    
    # get mut counts
    mutation_counts = muts_df.groupby([grouping_col_name, 'Patient_id']).size().reset_index(name='Mutation count')
    count_table = mutation_counts.groupby([grouping_col_name, 'Mutation count']).size().reset_index(name='Number of patients')
    pivot_table = count_table.pivot(index='Mutation count', columns=grouping_col_name, values='Number of patients').fillna(0).astype(int)
    pivot_table.columns = [f'Number of patients {col.lower()}' for col in pivot_table.columns]
    pivot_table.reset_index(inplace=True)
    pivot_table[f"{group1_name}_fraction"] = pivot_table[f"Number of patients {group1_name.lower()}"] / group1_ntotal
    pivot_table[f"{group2_name}_fraction"] = pivot_table[f"Number of patients {group2_name.lower()}"] / group2_ntotal
    
    # Get the number of patients with 0 mutations
    n_ch_positive_group1 = muts_df[muts_df[grouping_col_name] == group1_name]["Patient_id"].unique().shape[0]
    n_ch_positive_group2 = muts_df[muts_df[grouping_col_name] == group2_name]["Patient_id"].unique().shape[0]
    n_0_group1 = group1_ntotal - n_ch_positive_group1
    n_0_group2 = group2_ntotal - n_ch_positive_group2
    
    pts_with_0_muts_dict = {
        "Mutation count": 0, 
        f"Number of patients {group1_name.lower()}": n_0_group1, 
        f"Number of patients {group2_name.lower()}": n_0_group2,
        f"{group1_name}_fraction": n_0_group1/group1_ntotal,
        f"{group2_name}_fraction": n_0_group2/group2_ntotal}
    
    pts_with_0_muts_df = pd.DataFrame.from_dict(pts_with_0_muts_dict, orient='index').T
    pivot_table = pd.concat([pivot_table, pts_with_0_muts_df]).sort_values(by = "Mutation count").reset_index(drop = True)
    pivot_table["xpos"] = pivot_table.index
    # pivot_table.loc[pivot_table["xpos"] == 12, "xpos"] = 13
    
    # For the swarm plot get the vafs of all muts
    muts_df = muts_df.merge(mutation_counts, how = "left")
    group1_muts = muts_df[muts_df[grouping_col_name] == group1_name].reset_index(drop = True)
    group2_muts = muts_df[muts_df[grouping_col_name] == group2_name].reset_index(drop = True)
    
    # Plot them
    ax2 = ax.twinx()
    for i, row in pivot_table.iterrows():
        ax.bar(row["xpos"]-0.2, row[f"{group1_name}_fraction"], color=group1_color, width = 0.4, edgecolor = "None")
        ax.text(row["xpos"] - 0.2, 0.3, int(row[f"Number of patients {group1_name.lower()}"]), ha='center', va='center', fontsize=5, color='black')
        
        ax.bar(row["xpos"]+0.2, row[f"{group2_name}_fraction"], color=group2_color, width = 0.4, edgecolor = "None")
        ax.text(row["xpos"] + 0.2, 0.3, int(row[f"Number of patients {group2_name.lower()}"]), ha='center', va='center', fontsize=5, color='black')
    
    # plot the vafs in logscale
    muts_df["VAF_n_log"] = np.log10(muts_df["VAF_n"].replace(0, np.nan))
    for i, row in muts_df.iterrows():
        if row[grouping_col_name] == group1_name: 
            offset = 0.2
        else:
            offset = -0.2
        jitter = np.random.uniform(-0.08, 0.08, 1)
        ax2.scatter(row["Mutation count"]+offset+jitter[0], row["VAF_n_log"], color="black", s = 0.15, alpha = 0.7)
    
    # plot the vafs for bladder
    # for i, row in bladder_muts.iterrows():
    #     jitter = np.random.normal(-0.05, 0.05, 1)
    #     print(row["Mutation count"]-0.2+jitter, row["VAF_n"])
    #     ax2.scatter(row["Mutation count"]-0.2+jitter, row["VAF_n"], color="black", s = 2)
    
    # Aes
    for a in [ax, ax2]:
        a.spines["top"].set_visible(False)
    # x ticks
    ax.set_xticks(pivot_table["xpos"])
    ax.set_xticklabels(pivot_table["xpos"])
    ax.set_xlabel("Number of CH mutations")
    ax.set_ylabel("% of patients in arm")
    ax.tick_params(axis='x', bottom=False)
    ax2.tick_params(axis='x', pad=-5)
    # ax.tick_params(axis="both", direction="out", which="both", left=True, bottom=True , colors='k')    
    ax.set_ylim((0, 0.3))
    ax.set_yticks([0, 0.1, 0.2, 0.3])
    ax.set_yticklabels(["0", "10", "20", "30"])
    ax2.set_ylim((np.log10(0.25), np.log10(60)))
    ax2.set_yticks([np.log10(0.25), np.log10(1), np.log10(2), np.log10(10), np.log10(50)])
    ax2.set_yticklabels([".25", "1", "2", "10", "50"])
    ax2.set_ylabel("WBC VAF %")
    
    # Add legend
    legend_colors = color_dict.values()
    legend_labels = color_dict.keys()
    legend_handles = [plt.Line2D([0], [0], marker='s', color=color, label=label, markersize=5, linestyle='') for color, label in zip(legend_colors, legend_labels)]
    ax.legend(handles=legend_handles, loc="upper right", frameon=False, fontsize = 8, handlelength=2, handletextpad = 0.5)
    return([ax, ax2])

def plot_gene_counts_grouped_single_ax(muts_df, grouping_col_name, color_dict, ax, gene_list = None, horizontal = False, group1_ntotal = 95, group2_ntotal = 81, use_entire_denom = True, perform_z_test = False):
    """
    Plots the number of mutations each individual gene has in the cohort as a BAR CHART. If unique, non unique mutations in the same gene are counted only once.
    Plots kidney and bladder separatelyin grouped bar charts. 
    """
    # Extract names of groups
    group1_name = "LuPSMA"
    group2_name = "Cabazitaxel"
    
    # Determine colors for each group
    group1_color = color_dict[group1_name]
    group2_color = color_dict[group2_name]
    
    # Generate our denominators, which will be the number of CH+ people.
    if use_entire_denom:
        n_group1 = group1_ntotal
        n_group2 = group2_ntotal
    else:
        n_group1 = muts_df[muts_df[grouping_col_name] == group1_name]["Patient_id"].unique().shape[0]
        n_group2 = muts_df[muts_df[grouping_col_name] == group2_name]["Patient_id"].unique().shape[0]
    
    if gene_list is not None: 
        df = muts_df[muts_df["Gene"].isin(gene_list)].reset_index(drop = True)
    
    # subset
    group1_df = df[df[grouping_col_name] == group1_name].reset_index(drop = True)[["Patient_id", "Gene"]].drop_duplicates()["Gene"].value_counts().reset_index().rename(columns = {"index": "Gene", "Gene": f"Counts_{group1_name}"})
    group2_df = df[df[grouping_col_name] == group2_name].reset_index(drop = True)[["Patient_id", "Gene"]].drop_duplicates()["Gene"].value_counts().reset_index().rename(columns = {"index": "Gene", "Gene": f"Counts_{group2_name}"})
    
    # Go from counts to fractions
    group1_df[f"frac_{group1_name}"] = (group1_df[f"Counts_{group1_name}"]/n_group1)*100
    group2_df[f"frac_{group2_name}"] = (group2_df[f"Counts_{group2_name}"]/n_group2)*100
    
    combined = group1_df.merge(group2_df, on = "Gene", how = "outer").reset_index().fillna(0) # gene order will be determined based on bladder df
    gene_rank = {gene: rank for rank, gene in enumerate(gene_list)}
    combined['index'] = combined['Gene'].map(gene_rank)
    combined = combined.sort_values(by = "index")
    combined[f"Counts_{group1_name}"] = combined[f"Counts_{group1_name}"].astype(int)
    combined[f"Counts_{group2_name}"] = combined[f"Counts_{group2_name}"].astype(int)
        
    # plotting
    for i, row in combined.iterrows():
        # Plot bars
        if horizontal:
            ax.barh(row["index"]+0.2, row[f"frac_{group1_name}"], color=group1_color, height=0.4, edgecolor="None")
            ax.barh(row["index"]-0.2, row[f"frac_{group2_name}"], color=group2_color, height=0.4, edgecolor="None")
            
            ax.text(row[f"frac_{group1_name}"]+4, row["index"]+0.2, row[f"Counts_{group1_name}"], va='center', color='black', fontsize=5)
            ax.text(row[f"frac_{group2_name}"]+4, row["index"]-0.2, row[f"Counts_{group2_name}"], va='center', color='black', fontsize=5)
        else:
            ax.bar(row["index"]+0.2, row[f"frac_{group1_name}"], color=group2_color, width = 0.4, edgecolor = "None")
            ax.bar(row["index"]-0.2, row[f"frac_{group2_name}"], color=group2_color, width = 0.4, edgecolor = "None")
            
            ax.text(row["index"]+0.2, row[f"frac_{group1_name}"]+2, row[f"Counts_{group1_name}"], ha='center', color='black', fontsize = 5)
            ax.text(row["index"]-0.2, row[f"frac_{group2_name}"]+2, row[f"Counts_{group2_name}"], ha='center', color='black', fontsize = 5)    
    
    # aesthetics
    if horizontal:
        # ax.tick_params(axis='x', bottom=False, labelbottom=False)
        ax.tick_params(axis='y', left=False, labelleft=False)
        ax.spines[["top", "left"]].set_visible(False)
        tick_df = combined[["index", "Gene"]].sort_values(by = "index")
        ax.set_yticks(tick_df["index"].tolist())
        ax.set_yticklabels(tick_df["Gene"].tolist(), rotation = 0, fontstyle = "italic", va = "center", fontsize = 8, ha = "center")
        ax.tick_params(axis='y', which='both', pad=20)  # Adjust padding between tick labels and axis
        ax.tick_params(axis = "y", which='both', length=0)
        ax.set_xlabel("% mutated in CH+ pts")
        ax.set_xlim((0, 60))
        ax.set_xticks([0, 15, 30, 45, 60])
        ax.set_xticklabels(["0", "15", "30", "45", "60"])
        ax.set_ylim((-0.6, 14.6))
    else:
        ax.spines[["top", "right"]].set_visible(False)
        tick_df = combined[["index", "Gene"]].sort_values(by = "index")
        ax.set_xticks(tick_df["index"].tolist())
        ax.set_xticklabels(tick_df["Gene"].tolist(), rotation = 90, fontstyle = "italic", va = "center", fontsize = 8)
        ax.tick_params(axis='x', which='both', pad=20)  # Adjust padding between tick labels and axis
        ax.tick_params(axis = "x", which='both', length=0)
        ax.set_ylabel("% mutated in CH+ pts")
        ax.set_ylim((0, 60))
        ax.set_yticks([0, 15, 30, 45, 60])
        ax.set_yticklabels(["0", "15", "30", "45", "60"])
        ax.set_xlim((-0.6, 14.6))
    
    # Calculate proportions for each gene in bladder and kidney cancers
    combined[f'Proportion_{group1_name}'] = combined[f"Counts_{group1_name}"] / combined[f"Counts_{group1_name}"].sum()
    combined[f'Proportion_{group2_name}'] = combined[f"Counts_{group2_name}"] / combined[f"Counts_{group2_name}"].sum()
    
    # Perform the Z-test for each gene
    if perform_z_test:
        p_values_ztest = []
        nobs = np.array([group1_ntotal, group2_ntotal])
        for i, (group1_count, group2_count) in enumerate(zip(combined[f"Counts_{group1_name}"], combined[f"Counts_{group2_name}"])):
            counts = np.array([group1_count, group2_count])
            _, p_val = proportions_ztest(counts, nobs, alternative="larger")  # alternative="larger" for one-sided test
            p_values_ztest.append(p_val)
        
        adjusted_p_values_ztest = multipletests(p_values_ztest, method="fdr_bh")[1]
        adjusted_p_values_ztest = [float(p) for p in adjusted_p_values_ztest]
        # 3.064581261712751e-05
        # Plotting with significance marks
        # max_val_ax = ax.get_ylim()[1]
        for i, p_val in enumerate(adjusted_p_values_ztest):
            # Get the bladder count, the p values will be annotated right on top of it.
            i = i+1
            group2_count = combined[combined["index"] == i][f"frac_{group2_name}"].iloc[0]
            if horizontal:
                if p_val < 0.05:
                    formatted_p_value = round(float(p_val), 3)
                    formatted_text = f"p={formatted_p_value}"
                    ax.text(group2_count+4, i-0.2, formatted_text, ha='right', va='bottom', color='black', fontsize = 6, fontdict={'family': 'Arial'})
                elif p_val < 0.001: 
                    coefficient = f"{value:.2e}".split('e')[0]  # e.g., 2.59
                    exponent = int(f"{value:.2e}".split('e')[1])  # e.g., -5
                    formatted_p_value = rf"${coefficient} \times 10^{{{exponent}}}$"
                    formatted_text = f"p={formatted_p_value}"
                    ax.text(group2_count+4, i-0.2, formatted_text, ha='right', va='bottom', color='black', fontsize = 6, fontdict={'family': 'Arial'})
            else:
                if p_val < 0.001:
                    ax.text(i + 0.2, group2_count+5, "***", ha='center', va='bottom', color='black')
                elif p_val < 0.01:
                    ax.text(i + 0.2, group2_count+5, "**", ha='center', va='bottom', color='black')
                elif p_val < 0.05:
                    ax.text(i + 0.2, group2_count+5, "*", ha='center', va='bottom', color='black')
        print(p_values_ztest)
    
    # Generate the gene_order dict for the subsequent plot
    gene_order = combined[["Gene", "index"]]
    
    # Add legend
    legend_colors = color_dict.values()
    legend_labels = color_dict.keys()
    legend_handles = [plt.Line2D([0], [0], marker='s', color=color, label=label, markersize=5, linestyle='') for color, label in zip(legend_colors, legend_labels)]
    ax.legend(handles=legend_handles, loc="lower left", frameon=False, fontsize = 8, handlelength=2, handletextpad = 0.5)
    
    if horizontal:
        ax.invert_xaxis()
        ax.set_xlabel("% mutated in CH+ pts")
    return(ax, gene_order)

def plot_vafs_swarm(muts_df, grouping_col_name, ax_swarm, gene_list = None, revert_vafs = False, gene_order = None, add_violins = True, horizontal = False):
    """
    Plots the VAFs as swarm.
    """
    group1_name = "LuPSMA"
    group2_name = "Cabazitaxel"
    
    # Determine colors
    color_dict = {"LuPSMA": "crimson", "Cabazitaxel": "yellowgreen"}
    group1_color = color_dict[group1_name]
    group2_color = color_dict[group2_name]
    
    # Subset to genes of interest
    if gene_list is not None: 
        df = muts_df[muts_df["Gene"].isin(gene_list)].reset_index(drop = True)
    
    if "VAF_n" in df.columns:
        vaf_col = "VAF_n"
    else: 
        vaf_col = "VAF"
    df = df[["Gene", vaf_col, grouping_col_name]]
    
    if gene_order is not None:
        df = df.merge(gene_order, how = "left")
    
    # if revert_vafs:
    #     df["VAF_n"] = df["VAF_n"]*-1
    # Convert VAFs to a log scale, keeping them negative
    df["VAF_n_log"] = np.log10(df[vaf_col].replace(0, np.nan))  # Replace 0 with NaN to avoid log(0)
    
    group1_vafs = df[df[grouping_col_name] == group1_name]
    group2_vafs = df[df[grouping_col_name] == group2_name]
        
    for i, row in group1_vafs.iterrows():
        jitter = np.random.uniform(-0.05, 0.05, 1)
        ax_swarm.scatter(row["VAF_n_log"], row["index"]+0.2+jitter, color="black", alpha = 0.7, s = 0.15, zorder = 1) # group 1
    
    for i, row in group2_vafs.iterrows():
        jitter = np.random.normal(-0.05, 0.05, 1)
        ax_swarm.scatter(row["VAF_n_log"], row["index"]-0.2+jitter, color="black", alpha = 0.7, s = 0.15, zorder = 1) # group 2
    
    # ax_swarm.set_yscale('log')
    ax_swarm.spines[["right", "top"]].set_visible(False)
    ax_swarm.set_xlim((np.log10(0.25), np.log10(50)))
    ax_swarm.set_xticks([np.log10(0.25), np.log10(1), np.log10(2), np.log10(10), np.log10(50)])
    ax_swarm.set_xticklabels([".25", "1", "2", "10", "50"])
    ax_swarm.tick_params(axis='y', labelleft=False, left = False)
    ax_swarm.set_xlabel("WBC VAF%")
    
    # Shades behind the swarmplot
    if add_violins:
        positions_group1 = np.arange(len(gene_list)) + 0.2
        positions_group2 = np.arange(len(gene_list)) - 0.2
        width = 0.3 
        x_min, x_max = ax_swarm.get_xlim() # Get the current y limits (assuming your axis has been set up already
        bar_height = x_max - x_min #  Calculate the height for the bars
        
        # Add rectangles
        for pos in positions_group1:
            ax_swarm.barh(pos, bar_height, height = width, left = x_min, color=group1_color, alpha=0.1, edgecolor = None, zorder = 0)
        
        for pos in positions_group2:
            ax_swarm.barh(pos, bar_height, height = width, left = x_min, color=group2_color, alpha=0.1, edgecolor = None, zorder = 0)
    
    return(ax_swarm)

def annotate_mutation_status_lu(mutations_df, cohort, PATH_sample_information, annotate_what, timepoint, annotate_gene = False): 
    """
    Given a mutation mutations_df and a list of all pts in the cohort, annotate the status of patients in the mutation list.
    Intended for survival analysis.
    annotate_what: Either provide ctDNA or CHIP.
    annotate_gene: If a gene name is provided it will annotate the status of that gene. A list can also be given.
    """
    if isinstance(annotate_gene, str):
        mutations_df = mutations_df[mutations_df["Gene"] == annotate_gene].reset_index(drop = True)
    elif isinstance(annotate_gene, list):
        mutations_df = mutations_df[mutations_df["Gene"].isin(annotate_gene)].reset_index(drop=True)
    
    mutations_df = mutations_df[["Patient_id", "Timepoint_t", "Sample_name_t", "Gene", "VAF_n", "Protein_annotation"]]
    mutations_df[annotate_what + " status"] = "Positive"
    mutations_df = mutations_df[["Patient_id", annotate_what + " status"]].drop_duplicates().reset_index(drop = True)
    all_pts = pd.read_csv(PATH_sample_information, sep="\t")
    all_pts = all_pts[(all_pts["Timepoint"] == timepoint.capitalize()) & (all_pts["Cohort"] == cohort)].rename(columns = {"Timepoint": "Timepoint_t"})
    mutations_df = all_pts[["Patient_id", "Cohort", "Timepoint_t"]].drop_duplicates().merge(mutations_df, how = "left")
    if cohort != "Both":
        mutations_df = mutations_df[mutations_df["Cohort"] == cohort]
    mutations_df[annotate_what + " status"] = mutations_df[annotate_what + " status"].fillna("Negative")
    mutations_df = mutations_df.drop_duplicates().reset_index(drop = True)[["Patient_id", "Cohort", "Timepoint_t", annotate_what+" status"]]
    return(mutations_df)

def age_vs_CH_presence_grouped_ch_thresholds(muts_df, grouping_col_name, dict_box_colors, color_dict, age_df, PATH_sample_information, ax, fontsize = 10):
    """
    VERY SIMILAR TO age_vs_CH_presence. ONLY DIFFERENCE IS PLOTS BLADDER AND KIDNEY TOGETHER.
    Boxplot. Plotting the ages on the y-axis of patients, grouped into CH+ and CH-. Give a df consisting of only one timepoint, baseline or OT. 
    test_to_use: use either 'MWU' or 'T test'
    """
    group1_name = muts_df[grouping_col_name].unique()[0]
    group2_name = muts_df[grouping_col_name].unique()[1]
    
    # Determine colors
    group1_color_box = dict_box_colors[group1_name]
    group2_color_box = dict_box_colors[group2_name]
    
    group1_color = color_dict[group1_name]
    group2_color = color_dict[group2_name]
    
    dict_offset = {group1_name: -0.2, group2_name: 0.2}
    
    # dict_box_colors = {"Bladder": (143/255, 215/255, 239/255, 255/255), "Kidney": (239/255, 169/255, 143/255, 255/255)}
    # dict_scatter_colors = {"Bladder": "deepskyblue", "Kidney": "orangered"}
    # dict_offset = {"Bladder": -0.2, "Kidney": 0.2}
    
    # Data manipulation
    all_ch = annotate_mutation_status_lu(muts_df, "TheraP", PATH_sample_information, annotate_what = "CHIP", timepoint= "Baseline")
    all_ch = all_ch[all_ch["Timepoint_t"] == "Baseline"].reset_index(drop = True).merge(age_df, how = "left").dropna()
    
    # Add the arm information
    path_arm = "/groups/wyattgrp/users/amunzur/lu_chip/resources/clinical_tables/arm_info.tsv"
    arm_df = pd.read_csv(path_arm, sep = "\t")
    all_ch = all_ch.drop("Cohort", axis = 1).merge(arm_df)
    ch_neg = all_ch[(all_ch["CHIP status"] == "Negative")]
    
    max_ch_df = muts_df[["Patient_id", "VAF_n"]].groupby("Patient_id")["VAF_n"].max().reset_index().merge(arm_df)
    
    ch1 = max_ch_df[max_ch_df["VAF_n"] < 2].merge(age_df).assign(**{"CHIP status": "Positive"})
    ch2 = max_ch_df[(max_ch_df["VAF_n"] >= 2) & (max_ch_df["VAF_n"] < 10)].merge(age_df).assign(**{"CHIP status": "Positive"})
    ch3 = max_ch_df[max_ch_df["VAF_n"] >= 10].merge(age_df).assign(**{"CHIP status": "Positive"})
    
    flierprops_dict = dict(marker='o', markersize=5, markeredgecolor='black', linestyle='None')
    whiskerprops_dict =dict(color='black')
    medianprops_dict = dict(color='black')
    capprops_dict = dict(color='black')
    
    # Plot the boxes
    x_tick_labels_list = []
    x_ticks_list = []
    for i, df in enumerate([ch_neg, ch1, ch2, ch3]):
        df["CHIP status"] = df["CHIP status"].replace({"Negative": "CH-", "Positive": "CH+"})
        if i != 0:
            df = df[df["CHIP status"] == "CH+"]
        for group in [group1_name, group2_name]:
            df_diagnosis = df[df[grouping_col_name] == group]
            boxprops_dict = dict(facecolor=dict_box_colors[group], edgecolor='black', linewidth = 0.7)  
            boxplot = ax.boxplot(df_diagnosis["Age"], positions = [i+dict_offset[group]], flierprops = flierprops_dict, boxprops = boxprops_dict, medianprops = medianprops_dict, capprops = capprops_dict, widths = 0.3, showfliers = False, patch_artist = True)
            ax.scatter(np.random.uniform(i+dict_offset[group]-0.08, i+dict_offset[group]+0.08, len(df_diagnosis["Age"])), df_diagnosis["Age"], s = 3, color = color_dict[group], alpha = 1, zorder = 100)
            n_patients = len(df_diagnosis["Age"])
            x_tick_labels_list.append(f"n={n_patients}")
            x_ticks_list.append(i+dict_offset[group])
    # Aes
    ax.spines[["top", "right"]].set_visible(False)
    ax.set_ylabel("Age at baseline draw")
    ax.set_ylim((40, 90))
    ax.set_yticks([40, 50, 60, 70, 80, 90])
    ax.set_yticklabels(["40", "50", "60", "70", "80", "90"])
    ax.set_xticks(x_ticks_list)
    ax.set_xticklabels(x_tick_labels_list, fontsize = 7)
    
    # add text
    text_pos_dict = {0: "CH-", 1: "0.25%-2%", 2: "2%-10%", 3: ">10%"}
    for index, (key, value) in enumerate(text_pos_dict.items()):
        ax.text(key, 90, value, ha='center', va='top', fontsize=9) 
    
    # do MWU between various groups
    for group in [group1_name, group2_name]:
        age_ch_neg = ch_neg[ch_neg[grouping_col_name] == group]["Age"]
        age_ch1 = ch1[ch1[grouping_col_name] == group]["Age"]
        age_ch2 = ch2[ch2[grouping_col_name] == group]["Age"]
        age_ch3 = ch3[ch3[grouping_col_name] == group]["Age"]
        
        mwu1 = "{:.2e}".format(mannwhitneyu(age_ch_neg, age_ch1)[1])
        mwu2 = "{:.2e}".format(mannwhitneyu(age_ch_neg, age_ch2)[1])
        mwu3 = "{:.2e}".format(mannwhitneyu(age_ch_neg, age_ch3)[1])     
        
        # np.median(age_ch_neg)
        # np.median(age_ch1)
        # np.median(age_ch2)
        # np.median(age_ch3)
        
        print(group)
        print(f"CH neg vs 0.25%-2% p={mwu1}")
        print(f"CH neg vs 2%-10% p={mwu2}")
        print(f"CH neg vs >10% p={mwu3}")
        # "{:.2e}".format(mannwhitneyu(age_ch1, age_ch3)[1])
        # "{:.2e}".format(mannwhitneyu(age_ch_neg, age_ch3)[1])
    
    # Add legend
    legend_colors = color_dict.values()
    legend_labels = color_dict.keys()
    legend_handles = [plt.Line2D([0], [0], marker='s', color=color, label=label, markersize=5, linestyle='') for color, label in zip(legend_colors, legend_labels)]
    ax.legend(handles=legend_handles, loc="lower right", frameon=False, fontsize = 8, handlelength=2, handletextpad = 0.5)
    
    return(ax)

def CH_presence_line_thresholds(muts_df, age_df, PATH_sample_information, ax, axforest, fontsize = 10):
    """
    Plos the fraction of the cohort that is CH positive in each age.
    """
    import scipy.stats as stats
    sample_info = pd.read_csv(PATH_sample_information, sep = "\t")
    sample_info = sample_info[(sample_info["Timepoint"] == "Baseline")].reset_index(drop = True)
    
    # Add arm info
    path_arm = "/groups/wyattgrp/users/amunzur/lu_chip/resources/clinical_tables/arm_info.tsv"
    arm_df = pd.read_csv(path_arm, sep = "\t")
    
    # Generate dfs of varying thresholds and annotate their CHIP mutation status
    ch1 = annotate_mutation_status_lu(muts_df[(muts_df["VAF_n"] <= 2)], "TheraP", PATH_sample_information, annotate_what = "CHIP", timepoint= "Baseline").merge(age_df, how = "left").dropna().merge(arm_df)
    ch2 = annotate_mutation_status_lu(muts_df[(muts_df["VAF_n"] > 2) & (muts_df["VAF_n"] <= 10)], "TheraP", PATH_sample_information, annotate_what = "CHIP", timepoint= "Baseline").merge(age_df, how = "left").dropna().merge(arm_df)
    ch3 = annotate_mutation_status_lu(muts_df[muts_df["VAF_n"] > 10], "TheraP", PATH_sample_information, annotate_what = "CHIP", timepoint= "Baseline").merge(age_df, how = "left").dropna().merge(arm_df)
    
    ch1 = ch1[ch1["Timepoint_t"] == "Baseline"].reset_index(drop = True)
    ch2 = ch2[ch2["Timepoint_t"] == "Baseline"].reset_index(drop = True)
    ch3 = ch3[ch3["Timepoint_t"] == "Baseline"].reset_index(drop = True)
    
    colors = {"0.25%-2%": "limegreen", "2%-10%": "#709bd0", ">10%": "#224193"}
    
    bins = [40, 60, 70, 80, 90]
    labels = ['40-59', '60-69', '70-79', '80-89']
    
    for j, (df, df_annot) in enumerate(zip([ch1, ch2, ch3], ["0.25%-2%", "2%-10%", ">10%"])):
        # Bin the ages
        df['Age_bin'] = pd.cut(df['Age'], bins=bins, labels=labels, right=False)
        df['Age_bin'] = pd.Categorical(df['Age_bin'], categories=labels, ordered=True)
        
        # Calculate fraction that is CH+
        fractions = []
        for i, group in df.groupby("Age_bin"):
            n_pos = group[group["CHIP status"] == "Positive"].shape[0]
            n_neg = group[group["CHIP status"] == "Negative"].shape[0]
            perc_pos = n_pos/(n_pos+n_neg)
            fractions.append(perc_pos)
        
        # LOGISTIC REGRESSOIN
        df['CHIP_binary'] = df['CHIP status'].apply(lambda x: 1 if x == 'Positive' else 0)
        df['Age_bin'] = pd.Categorical(df['Age_bin'], ordered=True) # If 'Age_Bin' is already categorical, convert it to numeric codes (ordered)
        df['Age_bin_numeric'] = df['Age_bin'].cat.codes  # Convert Age Bin to numeric
        
        # Step 1: Prepare X (predictor) and y (outcome)
        X = df['Age_bin_numeric']  # Age bins as predictor
        y = df['CHIP_binary']      # Binary CHIP status as outcome
        
        # Step 2: Add constant to X for the intercept
        X = sm.add_constant(X)
        
        # Step 3: Fit logistic regression model
        logit_model = sm.Logit(y, X)
        result = logit_model.fit()
        
        odds_ratios = pd.DataFrame({"Odds Ratio": np.exp(result.params),"p-value": result.pvalues})
        print(f"CH category: {df_annot}")
        print(odds_ratios)
        
        # Step 4: Calculate Odds Ratios and Confidence Intervals
        odds_ratios = pd.DataFrame({
            "Odds Ratio": np.exp(result.params),
            "p-value": result.pvalues,
            "Lower 95% CI": np.exp(result.params - 1.96 * result.bse),  # Lower bound
            "Upper 95% CI": np.exp(result.params + 1.96 * result.bse)   # Upper bound
        })
        OR = odds_ratios[odds_ratios.index == "Age_bin_numeric"]["Odds Ratio"][0]
        p_value = odds_ratios[odds_ratios.index == "Age_bin_numeric"]["p-value"][0]
        lower_ci = odds_ratios[odds_ratios.index == "Age_bin_numeric"]["Lower 95% CI"][0]
        upper_ci = odds_ratios[odds_ratios.index == "Age_bin_numeric"]["Upper 95% CI"][0]
        
        # Step 4. Plotting LR results as an inset plot
        age_bins_numeric_range = np.linspace(df['Age_bin_numeric'].min(), df['Age_bin_numeric'].max(), 100)
        X_plot = sm.add_constant(age_bins_numeric_range)
        y_pred_probs = result.predict(X_plot)
        ax.scatter(sorted(df['Age_bin_numeric'].unique()), fractions, color=colors[df_annot], label='Actual Fraction Positive', s = 10) #Plot the actual fraction of CH-positive individuals by age bin
        ax.plot(age_bins_numeric_range, y_pred_probs, color=colors[df_annot], label='Logistic Regression Fit', linewidth = 0.8) # Step 4: Plot the logistic regression predicted probabilities
        
        # Calculating error bars for each dot
        age_bin_stats = df.groupby('Age_bin').agg(count_positive=('CHIP_binary', 'sum'),count_total=('CHIP_binary', 'size')).reset_index()
        age_bin_stats['proportion_positive'] = age_bin_stats['count_positive'] / age_bin_stats['count_total']
        conf_intervals = []
        for index, row in age_bin_stats.iterrows():
            n = row['count_total']
            p = row['proportion_positive']
            
            if n > 0:
                ci_low, ci_high = stats.binom.interval(0.95, n=n, p=p)
                # Calculate the error margins
                conf_intervals.append((p - ci_low/n, ci_high/n - p))
            else:
                conf_intervals.append((0, 0))  # No error if no total
        conf_intervals = np.array(conf_intervals).T
        
        # Plotting the error bar
        df_plotting = pd.DataFrame({'Age_bin': age_bin_stats['Age_bin'],'Positive': age_bin_stats['proportion_positive']})
        ax.errorbar(df_plotting.index, df_plotting['Positive'], yerr=conf_intervals, fmt='none', color=colors[df_annot], capsize=2, elinewidth=0.5, capthick = 0.5)
        
        # Plotting the forest
        if axforest is not None:
            axforest.scatter(OR, j, color=colors[df_annot], s = 15)
            axforest.text(3, j, round(p_value, 4), ha='center', va='center', fontsize=7, color = "black")
            axforest.errorbar(OR, j, xerr = [[OR - lower_ci], [upper_ci - OR]], fmt='none', color=colors[df_annot], capsize=2, elinewidth=1, capthick = 1)
    
    # Annotate the number of people in each age group below the x-tick labels
    patient_ages = sample_info.merge(age_df).drop_duplicates(["Patient_id", "Age"])
    patient_ages["Age_bin"] = pd.cut(patient_ages['Age'], bins=bins, labels=labels, right=False)
    patient_ages['Age_bin'] = pd.Categorical(patient_ages['Age_bin'], categories=labels, ordered=True)
    patient_ages = patient_ages["Age_bin"].value_counts().reset_index()
    patient_ages['index'] = pd.Categorical(patient_ages['index'], categories=labels, ordered=True)
    patient_ages = patient_ages.sort_values('index')
    x_ticklabellist = patient_ages["index"].astype(str) + "\nn=" + patient_ages["Age_bin"].astype(str).tolist()    
    ax.set_xticks([0, 1, 2, 3])
    ax.set_xticklabels(x_ticklabellist)
    # AES
    ax.set_ylim((-0.1, 1.1))
    ax.set_yticks([0, 0.2, 0.4, 0.6, 0.8, 1])
    ax.set_yticklabels(["0", "20", "40", "60", "80", "100"], fontsize=10)
    ax.set_xlabel('Age at baseline draw')
    ax.set_ylabel('% of patients with CH')
    ax.spines[["top", "right"]].set_visible(False)
    
    # add legend
    legend_colors = colors.values()
    legend_labels = colors.keys()
    legend_handles = [plt.Line2D([0], [0], marker='o', color=color, label=label, markersize=5, linestyle='') for color, label in zip(legend_colors, legend_labels)]
    ax.legend(handles=legend_handles, loc="upper right", frameon=False, fontsize = 8, handletextpad=0.1, title = "CH VAF interval")
    
    if axforest is not None:
        axforest.axvline(1, color='black', linestyle='--', linewidth = 0.5)
        axforest.set_ylim((-0.1, 2.1))
        axforest.set_yticks([0, 1, 2])
        axforest.tick_params(axis='y', left=False, labelleft=False)
        axforest.spines[["top", "right", "left"]].set_visible(False)
        axforest.set_yticklabels(["0.25%-2%", "2%-10%", ">10%"])
        axforest.set_xlabel("OR")
        return(ax, axforest)
    else:
        return(ax)

def plot_gene_counts_grouped_single_ax_progression(muts_df, path_sample_info, grouping_col_name, color_dict, ax, gene_list = None, horizontal = False, perform_z_test = False):
    """
    Plots the number of mutations each individual gene has in the cohort as a BAR CHART. If unique, non unique mutations in the same gene are counted only once.
    Plots kidney and bladder separatelyin grouped bar charts. 
    """
    sample_info = pd.read_csv(path_sample_info, sep = "\t")
    progression_samples = sample_info[(sample_info["Cohort"] == "TheraP") & (sample_info["Timepoint"] == "FirstProgression")]
    
    ntotal_caba_progression = progression_samples[progression_samples["Arm"] == "Cabazitaxel"].shape[0]
    ntotal_lu_progression = progression_samples[progression_samples["Arm"] == "LuPSMA"].shape[0]
    
    group1_name = "LuPSMA"
    group2_name = "Cabazitaxel"
        
    if gene_list is not None: 
        df = muts_df[muts_df["Gene"].isin(gene_list)].reset_index(drop = True)
    
    # subset
    group1_df = df[df["Arm"] == group1_name].reset_index(drop = True)[["Patient_id", "Gene"]].drop_duplicates()["Gene"].value_counts().reset_index().rename(columns = {"index": "Gene", "Gene": f"Counts_{group1_name}"})
    group2_df = df[df["Arm"] == group2_name].reset_index(drop = True)[["Patient_id", "Gene"]].drop_duplicates()["Gene"].value_counts().reset_index().rename(columns = {"index": "Gene", "Gene": f"Counts_{group2_name}"})
    
    # Go from counts to fractions
    group1_df[f"frac_{group1_name}"] = (group1_df[f"Counts_{group1_name}"]/ntotal_lu_progression)*100
    group2_df[f"frac_{group2_name}"] = (group2_df[f"Counts_{group2_name}"]/ntotal_caba_progression)*100
    
    combined = group1_df.merge(group2_df, on = "Gene", how = "outer").fillna(0).sort_values(by = "frac_LuPSMA").reset_index(drop=True).reset_index()
    gene_rank = {gene: rank for rank, gene in enumerate(gene_list)}
    combined['index'] = combined['Gene'].map(gene_rank)
    combined = combined.sort_values(by = "index")
    combined["Counts_Cabazitaxel"] = combined["Counts_Cabazitaxel"].astype(int)
    combined["Counts_LuPSMA"] = combined["Counts_LuPSMA"].astype(int)
    
    # plotting
    for i, row in combined.iterrows():
            ax.barh(row["index"]+0.2, row[f"frac_{group1_name}"], color=color_dict[group1_name], height=0.4, edgecolor="None")
            ax.barh(row["index"]-0.2, row[f"frac_{group2_name}"], color=color_dict[group2_name], height=0.4, edgecolor="None")
            
            ax.text(row[f"frac_{group1_name}"]+2, row["index"]+0.2, row[f"Counts_{group1_name}"], va='center', color='black', fontsize=5)
            ax.text(row[f"frac_{group2_name}"]+2, row["index"]-0.2, row[f"Counts_{group2_name}"], va='center', color='black', fontsize=5)
    
    # aesthetics
    # ax.tick_params(axis='x', bottom=False, labelbottom=False)
    ax.set_xlabel("% of patients with progression samples")
    ax.tick_params(axis='y', left=False, labelleft=False)
    ax.spines[["top", "left"]].set_visible(False)
    tick_df = combined[["index", "Gene"]].sort_values(by = "index")
    ax.set_yticks(tick_df["index"].tolist())
    ax.set_yticklabels(tick_df["Gene"].tolist(), rotation = 0, fontstyle = "italic", va = "center", fontsize = 8, ha = "center")
    ax.tick_params(axis='y', which='both', pad=20)  # Adjust padding between tick labels and axis
    ax.tick_params(axis = "y", which='both', length=0)
    ax.set_xlim((0, 45))
    ax.set_xticks([0, 15, 30, 45])
    ax.set_xticklabels(["0", "15", "30", "45"])
    ax.set_ylim((-0.6, 14.6))
    
    # Calculate proportions for each gene in bladder and kidney cancers
    combined[f'Proportion_{group1_name}'] = combined[f"Counts_{group1_name}"] / combined[f"Counts_{group1_name}"].sum()
    combined[f'Proportion_{group2_name}'] = combined[f"Counts_{group2_name}"] / combined[f"Counts_{group2_name}"].sum()
    
    # Perform the Z-test for each gene
    if perform_z_test:
        p_values_ztest = []
        nobs = np.array([ntotal_caba_progression, ntotal_lu_progression])
        for i, (group1_count, group2_count) in enumerate(zip(combined[f"Counts_{group1_name}"], combined[f"Counts_{group2_name}"])):
            counts = np.array([group1_count, group2_count])
            _, p_val = proportions_ztest(counts, nobs, alternative="larger")  # alternative="larger" for one-sided test
            p_values_ztest.append(p_val)
        
        adjusted_p_values_ztest = multipletests(p_values_ztest, method="fdr_bh")[1]
        adjusted_p_values_ztest = [float(p) for p in adjusted_p_values_ztest]
        
        for i, p_val in enumerate(adjusted_p_values_ztest):
            if p_val <= 0.05:
                group2_count = combined[combined["index"] == i][f"frac_{group1_name}"].iloc[0]
                if p_val == 5.269907520177619e-05:
                    formatted_text = "p=0.000053"
                else:
                    formatted_p_value = round(float(p_val), 3)
                    formatted_text = f"p={formatted_p_value}"
                ax.text(group2_count+2, i-0.2, formatted_text, ha='right', va='bottom', color='black', fontsize = 6, fontdict={'family': 'Arial'})
        print(adjusted_p_values_ztest)
    
    # Generate the gene_order dict for the subsequent plot
    gene_order = combined[["Gene", "index"]]
    
    # Add legend
    legend_colors = color_dict.values()
    legend_labels = color_dict.keys()
    legend_handles = [plt.Line2D([0], [0], marker='s', color=color, label=label, markersize=5, linestyle='') for color, label in zip(legend_colors, legend_labels)]
    ax.legend(handles=legend_handles, loc="lower left", frameon=False, fontsize = 8, handlelength=2, handletextpad = 0.5)
    
    ax.invert_xaxis()
    return(ax, gene_order)

def fraction_ch_pos_bar_chart(muts, path_sample_information, color_dict, ax, ax_title, which_timepoint = "Baseline"):
    """
    Shows fraction of patients that are CH+ and CH-. Uses the entire cohort size as denominator.
    """
    sample_info = pd.read_csv(path_sample_info, sep = "\t")
    sample_info = sample_info[sample_info["Timepoint"] == which_timepoint]
    
    arms = ['LuPSMA', 'Cabazitaxel']
    
    for i, arm in enumerate(arms):
        all_pts_in_arm = sample_info[(sample_info["Cohort"] == "TheraP") & (sample_info["Arm"] == arm)]["Patient_id"].unique()
        muts_subset = muts[muts["Patient_id"].isin(all_pts_in_arm)]
        CH_pos_pts = muts_subset["Patient_id"].unique()
        
        n_total = all_pts_in_arm.shape[0]
        n_ch_pos = CH_pos_pts.shape[0]
        
        fr_ch_pos = n_ch_pos/n_total*100
        
        # Plotting
        ax.bar(i, fr_ch_pos, color = color_dict[arm], edgecolor = None, width = 0.7)
        
        # Annotate on the bars the number of samples
        y_pos = 95
        ax.text(i, y_pos, f'{n_ch_pos}/{n_total}', ha='center', va='center', color="black", fontsize=10)
    
    # AES
    ax.set_ylim((0, 100))
    ax.set_yticks([0, 25, 50, 75, 100])
    ax.set_yticklabels(["0", "25", "50", "75", "100"], fontsize=10)
    ax.set_xlabel('Arm')
    ax.set_ylabel('% CH+')
    ax.set_xticks([0, 1])
    ax.set_xticklabels(["LuPSMA", "Cabazitaxel"])
    ax.spines[["top", "right"]].set_visible(False)
    ax.set_title(ax_title)
    
    return(ax)

def plot_mutcounts_per_gene(muts_df, genes_list, ax, ntotal_lu, ntotal_caba):
    """
    In a given muts list plots the number of patients having n mutations
    """
    muts_df_subset = muts_df[muts_df["Gene"].isin(genes_list)]
    muts_count = muts_df_subset[["Patient_id", "Gene", "Arm"]].value_counts().reset_index()
    muts_count.columns = ["Patient_id", "Gene", "Arm", "Counts"]
    muts_count.loc[muts_count["Counts"] > 1, "Counts"] = 2 # More than 1 gets 2
    
    muts_count = muts_count[["Gene", "Arm", "Counts"]].value_counts().reset_index()
    muts_count.columns = ["Gene", "Arm", "Count_per_pt", "pt_count"]
    
    ntotals_list = {'LuPSMA': ntotal_lu, 'Cabazitaxel': ntotal_caba}
    
    # Plotting
    arm_offset = {"LuPSMA": -0.2, "Cabazitaxel": 0.2}
    color_dict_counts = {
        "LuPSMA":{
            2: (220/255, 20/255, 60/255, 1),
            1: (220/255, 20/255, 60/255, 0.5)
            }, 
        "Cabazitaxel": {
            2: (154/255, 205/255, 50/255, 1),
            1: (154/255, 205/255, 50/255, 0.5)
            }
        }
    
    for i, gene in enumerate(genes_list):
        gene_df = muts_count[muts_count["Gene"] == gene]
        for arm in ["LuPSMA", "Cabazitaxel"]:
            arm_ntotal = ntotals_list[arm]
            arm_df = gene_df[gene_df["Arm"] == arm].sort_values(by = "Count_per_pt")
            arm_df["Fraction"] = (arm_df["pt_count"]/arm_ntotal)*100
            arm_colors = color_dict_counts[arm]
            
            # Initialize y val to make stacked bar chart
            y_val_sum = 0
            for j, row in arm_df.iterrows():
                count_color = arm_colors[row["Count_per_pt"]]
                ax.bar(i+arm_offset[arm], row["Fraction"], bottom=y_val_sum, color = count_color, edgecolor = None, width = 0.4)
                
                # Annotate patient numbers
                ax.text(i+arm_offset[arm], row["Fraction"]+y_val_sum-1, str(row["pt_count"]), ha='center', va='center', fontsize=7, color='black')
                y_val_sum += row["Fraction"]
    
    # AES
    # ax.set_ylim((0, 100))
    # ax.set_yticks([0, 20, 40, 60, 80, 100])
    # ax.set_yticklabels(["0", "20", "40", "60", "80", "100"], fontsize=10)
    # ax.set_xlabel('Arm')
    ax.set_ylabel('% of patients in arm')
    ax.set_xticks(range(0, len(genes_list)))
    ax.set_xticklabels(genes_list, rotation = 90, fontstyle="italic")
    ax.tick_params(axis='x', bottom=False, pad = 1)
    ax.spines[["top", "right"]].set_visible(False)
    
    # Add legend
    legend_colors = list(color_dict_counts["LuPSMA"].values()) + list(color_dict_counts["Cabazitaxel"].values())
    legend_labels = ["LuPSMA >1", "LuPSMA 1 mut", "Caba >1 mut", "Caba 1 mut"]
    legend_handles = [plt.Line2D([0], [0], marker='s', color=color, label=label, markersize=5, linestyle='') for color, label in zip(legend_colors, legend_labels)]
    ax.legend(handles=legend_handles, loc="upper right", ncol = 2, frameon=False, fontsize = 8, handlelength=2, handletextpad = 0.1)
    
    return(ax)

def plot_new_mutation_fraction_at_progression(df, path_sample_info, ax, color_dict): 
    """
    Plots fraction patients at each arm that develops a new CH mutation at progression
    """
    sample_info = pd.read_csv(path_sample_info, sep = "\t")
    progression_samples = sample_info[(sample_info["Cohort"] == "TheraP") & (sample_info["Timepoint"] == "FirstProgression")]
    progression_samples = progression_samples[["Patient_id", "Arm"]].drop_duplicates()
    
    # Return total counts
    ntotal_caba_progression = progression_samples[progression_samples["Arm"] == "Cabazitaxel"].shape[0]
    ntotal_lu_progression = progression_samples[progression_samples["Arm"] == "LuPSMA"].shape[0]
    
    # Return counts from mutations
    n_caba_progression = df[df["Arm"] == "Cabazitaxel"][["Patient_id", "Arm"]].drop_duplicates().shape[0]
    n_lu_progression = df[df["Arm"] == "LuPSMA"][["Patient_id", "Arm"]].drop_duplicates().shape[0]
    
    group1_name = "Cabazitaxel"
    group2_name = "LuPSMA"
        
    perc_new_caba = n_caba_progression/ntotal_caba_progression*100
    perc_new_lu = n_lu_progression/ntotal_lu_progression*100
    
    # Plotting
    ax.bar(1, perc_new_caba, color = color_dict["Cabazitaxel"], edgecolor = None, width = 0.6)
    ax.bar(2, perc_new_lu, color = color_dict["LuPSMA"], edgecolor = None, width = 0.6)
    
    ax.set_yticks([0, 25, 50, 75, 100])
    ax.set_yticklabels(["0", "25", "50", "75", "100"])
    ax.set_ylabel("% patients with new CH mutations\nat progression")
    
    ax.set_xticks([1, 2])
    ax.set_xticklabels(["Caba", "LuPSMA"])
    ax.set_xlabel("Arm")
    
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    
    # Annotate patient numbers
    ax.text(1, 95, f'n={arm_counts["Cabazitaxel"]}', fontsize=10, color='black', ha = "center")
    ax.text(2, 95, f'n={arm_counts["LuPSMA"]}', fontsize=10, color='black', ha = "center")
    
    return(ax)

def plot_PPM1D_vaf_at_two_timepoints(baseline_called, color_dict, gene, ax): 
    """
    Plots the median VAF of PPM1D mutations at baseline and progression as a bar chart in two arms. 
    """
    # Mutations in PPM1D present in both timepoints
    df = baseline_called[(baseline_called["Gene"] == gene) & (~baseline_called["Progression_VAF"].isna())]
    lu_vaf_baseline = df[df["Arm"] == "LuPSMA"]["VAF_t"].median()
    lu_vaf_prog = df[df["Arm"] == "LuPSMA"]["Progression_VAF"].median()
    caba_vaf_baseline = df[df["Arm"] == "Cabazitaxel"]["VAF_t"].median()
    caba_vaf_prog = df[df["Arm"] == "Cabazitaxel"]["Progression_VAF"].median()
    
    # Plot them in scatter
    ax.scatter(0, lu_vaf_baseline, s=25, edgecolor=None, color=color_dict["LuPSMA"], zorder = 100)
    ax.scatter(0, caba_vaf_baseline, s=25, edgecolor=None, color=color_dict["Cabazitaxel"], zorder = 100)
    ax.scatter(1, lu_vaf_prog, s=25, edgecolor=None, color=color_dict["LuPSMA"], zorder = 100)
    ax.scatter(1, caba_vaf_prog, s=25, edgecolor=None, color=color_dict["Cabazitaxel"], zorder = 100)
    
    # Plot lines in between
    ax.plot([0, 1], [lu_vaf_baseline, lu_vaf_prog], color=color_dict["LuPSMA"], linewidth=2)
    ax.plot([0, 1], [caba_vaf_baseline, caba_vaf_prog], color=color_dict["Cabazitaxel"], linewidth=2)
        
    ax.set_yticks([0, 2, 4, 6, 8, 10])
    ax.set_yticklabels(["0", "2", "4", "6", "8", "10"])
    ax.set_ylabel("Median VAF%")
    
    ax.set_xticks([0, 1])
    ax.set_xticklabels(["Base", "Prog"])
    ax.set_xlabel("Timepoint")
    
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    
    # ax.set_title("Median VAF change\nPPM1D")
    
    # Add legend
    legend_colors = color_dict.values()
    legend_labels = color_dict.keys()
    legend_handles = [plt.Line2D([0], [0], marker='o', color=color, label=label, markersize=5, linestyle='') for color, label in zip(legend_colors, legend_labels)]
    ax.legend(handles=legend_handles, loc="upper left", frameon=False, fontsize = 8, handlelength=2, handletextpad = 0.1)
    
    return(ax)

def plot_gene_vafs_as_boxp(df, vaf_col_name, color_dict_boxp, color_dict_scatter, ax, gene = None):
    """
    Plots the VAFs from  given col in a given gene as boxplot and swarm plot.
    """
    # Convert VAFs to log 
    df["VAF_log"] = np.log10(df[vaf_col_name].replace(0, np.nan))  # Replace 0 with NaN to avoid log(0)
    
    if gene is not None:
        df_gene = df[df["Gene"] == gene]
        
    # Plot boxplot
    flierprops_dict = dict(marker='o', markersize=5, markeredgecolor='black', linestyle='None')
    whiskerprops_dict =dict(color='black', linewidth = 0.7)
    medianprops_dict = dict(color='black')
    capprops_dict = dict(color='black')
    
    medians_list = []
    for i, arm in enumerate(["LuPSMA", "Cabazitaxel"]):
        boxprops_dict = dict(facecolor=color_dict_boxp[arm], edgecolor='black', linewidth = 0.7)  
        df_vafs = df_gene[df_gene["Arm"] == arm]["VAF_log"]
        ax.boxplot(df_vafs, positions = [i], flierprops = flierprops_dict, boxprops = boxprops_dict, medianprops = medianprops_dict, capprops = capprops_dict, widths = 0.3, showfliers = False, patch_artist = True)
        ax.scatter(np.random.uniform(i+0.2, i-0.2, len(df_vafs)), df_vafs, s = 3, color = color_dict_scatter[arm], alpha = 1, zorder = 100)
        
        # Annotate the ns
        n_muts = df_vafs.shape[0]
        ax.text(i, 6, f'n={n_muts}', fontsize=10, color='black', ha = "center")
        
    # ax.set_yticks([0, 1, 2, 3, 4, 5, 6])
    # ax.set_yticklabels(["0", "1", "2", "3", "4", "5", "6"])
    ax.set_ylabel("VAF%", fontsize = 10)
    
    ax.set_xticks([0, 1])
    ax.set_xticklabels(["LuPSMA", "Cabazitaxel"])
    
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
       
    return(ax)

def progression_age_boxp(path_sample_information, path_therap_clinical_data, ax, timepoint, scatter_color_dict, boxp_color_dict):
    """
    Plots ages at progression or baseline as boxp.
    """
    clinical_data = pd.read_csv(path_therap_clinical_data, sep = "\t", skiprows = 1)
    clinical_data = clinical_data[["Patient", "Age", "Arm received"]].rename(columns = {"Patient": "Patient_id", "Age": "Age at baseline"})
    
    if timepoint == "Baseline":
        lupsma_ages = clinical_data[clinical_data["Arm received"] == "LuPSMA"]["Age at baseline"]
        caba_ages = clinical_data[clinical_data["Arm received"] == "Cabazitaxel"]["Age at baseline"]
    else: 
        sample_info = pd.read_csv(path_sample_information, sep = "\t")
        sample_info = sample_info[sample_info["Cohort"] == "TheraP"]
        
        # Calculate age at progression
        sample_info_prog = sample_info[sample_info["Timepoint"] == "FirstProgression"][["Patient_id", "Date_collected"]].rename(columns = {"Date_collected": "Date_collected_progression"})
        sample_info_base = sample_info[sample_info["Timepoint"] == "Baseline"][["Patient_id", "Date_collected"]].rename(columns = {"Date_collected": "Date_collected_baseline"})
        sample_info_prog = sample_info_prog.merge(sample_info_base, how = "inner")
        
        sample_info_prog["Date_collected_progression"] = pd.to_datetime(sample_info_prog["Date_collected_progression"], format='%Y%b%d')
        sample_info_prog["Date_collected_baseline"] = pd.to_datetime(sample_info_prog["Date_collected_baseline"], format='%Y%b%d')
        sample_info_prog["Timediff"] = (sample_info_prog["Date_collected_progression"] - sample_info_prog["Date_collected_baseline"]).dt.days / 365.25
        sample_info_prog = sample_info_prog.merge(clinical_data, how = "left")
        sample_info_prog["Age at prog"] = sample_info_prog["Age at baseline"] + sample_info_prog["Timediff"]
        lupsma_ages = sample_info_prog[sample_info_prog["Arm received"] == "Cabazitaxel"]["Age at prog"]
        caba_ages = sample_info_prog[sample_info_prog["Arm received"] == "LuPSMA"]["Age at prog"]
    
    # PLOTTING
    flierprops_dict = dict(marker='o', markersize=5, markeredgecolor='black', linestyle='None')
    whiskerprops_dict =dict(color='black')
    medianprops_dict = dict(color='black')
    capprops_dict = dict(color='black')
    
    ax.boxplot(lupsma_ages, positions = [0], flierprops = flierprops_dict, boxprops = dict(facecolor=boxp_color_dict["LuPSMA"], edgecolor='black', linewidth = 0.7)  , medianprops = medianprops_dict, capprops = capprops_dict, widths = 0.3, showfliers = False, patch_artist = True)
    ax.boxplot(caba_ages, positions = [1], flierprops = flierprops_dict, boxprops = dict(facecolor=boxp_color_dict["Cabazitaxel"], edgecolor='black', linewidth = 0.7)  , medianprops = medianprops_dict, capprops = capprops_dict, widths = 0.3, showfliers = False, patch_artist = True)
    ax.scatter(np.random.uniform(-0.3, 0.3, len(lupsma_ages)), lupsma_ages, s = 6, color = scatter_color_dict["LuPSMA"], zorder = 100)
    ax.scatter(np.random.uniform(0.7, 1.3, len(caba_ages)), caba_ages, s = 6, color = scatter_color_dict["Cabazitaxel"], zorder = 100)
    
    # Annotate medians
    lu_median = round(np.median(lupsma_ages), 1)
    caba_median = round(np.median(caba_ages), 1)
    ax.text(0, lu_median, str(lu_median), ha='center', va='bottom', fontsize=7, color='black', zorder = 999)
    ax.text(1, caba_median, str(caba_median), ha='center', va='bottom', fontsize=7, color='black', zorder = 999)
    
    # AES
    ax.set_ylabel("Age at first progression", fontsize = 10)
    ax.set_xticks([0, 1])
    ax.set_xticklabels(["LuPSMA", "Cabazitaxel"], fontsize = 10)
    
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    return(ax)

def make_survival_curve_therap(surv_df_merged, stratify_by, event_col, duration_col, output_path, plot_title, print_HRs = True, pos_curve_label_positions = None, neg_curve_label_positions = None, xmax = None, add_legend = True, print_survival = True, xlabel = None, ax = None):
    """
    stratify_by = Name of the column in df to stratify the curves by.
    """
    import math
    label_keyword = re.search("(ctDNA|CHIP|CH)", stratify_by, flags=re.IGNORECASE).group()
    kmf_positive = KaplanMeierFitter()
    kmf_negative = KaplanMeierFitter()
    groups = surv_df_merged[stratify_by]
    ix_positive = (groups == 'Positive') | (groups == 1)
    
    neg_curve_color = "black"
    if label_keyword == "ctDNA": 
        pos_curve_color = "forestgreen"
    else:
        pos_curve_color = "red"
    if ax is None:
        fig, ax = plt.subplots(figsize=(5, 5))
    
    # Negative group
    kmf_negative.fit(durations=surv_df_merged[duration_col][~ix_positive], event_observed=surv_df_merged[event_col][~ix_positive], label=f"{label_keyword}(-)")
    kmf_negative.plot_survival_function(ax=ax, ci_show=False, show_censors=True, linewidth=0.8, censor_styles={'marker': '|', 'ms': 5, 'markerfacecolor': neg_curve_color, 'mew': 0.8}, color = neg_curve_color)
    
    # Positive group
    kmf_positive.fit(durations=surv_df_merged[duration_col][ix_positive], event_observed=surv_df_merged[event_col][ix_positive], label=f"{label_keyword}(+)")
    kmf_positive.plot_survival_function(ax=ax, ci_show=False, show_censors=True, linewidth=0.8, censor_styles={'marker': '|', 'ms': 5, 'markerfacecolor': pos_curve_color, 'mew': 0.8}, color = pos_curve_color)
        
    if print_survival is not None:
        if isinstance(print_survival, bool) and print_survival:
            ax_to_plot = ax
        elif isinstance(print_survival, plt.Axes):
            ax_to_plot = print_survival
            if kmf_positive.median_survival_time_ == math.inf: 
                median_survival_positive = "Not reached"
            else:
                median_survival_positive = round(kmf_positive.median_survival_time_)
            if kmf_negative.median_survival_time_ == math.inf:
                median_survival_negative = "Not reached"
            else:
                median_survival_negative = round(kmf_negative.median_survival_time_)
                # check if the median survivals are significantly different using log rank test
                results = logrank_test(surv_df_merged[duration_col][~ix_positive], surv_df_merged[duration_col][ix_positive], event_observed_A=surv_df_merged[event_col][~ix_positive], event_observed_B=surv_df_merged[event_col][ix_positive])
                p_value_logrank = f"{results.p_value:.2e}"
                ax_to_plot.text(0.3, 0.37, f"Log rank p={p_value_logrank}", transform=ax_to_plot.transAxes, fontsize=10, ha='left', va='center')
            
            ax_to_plot.text(1.1, 0.8, "Median", transform=ax_to_plot.transAxes, fontsize=10, ha='center', va='center')
            ax_to_plot.text(1.1, 0.67, "mo", transform=ax_to_plot.transAxes, fontsize=10, ha='center', va='center', fontstyle = "italic")
            
            ax_to_plot.text(0.3, 0.57, f'{label_keyword}+', transform=ax_to_plot.transAxes, fontsize=10, ha='left', va='center')
            ax_to_plot.text(0.3, 0.47, f'{label_keyword}-', transform=ax_to_plot.transAxes, fontsize=10, ha='left', va='center')
            
            ax_to_plot.text(1.1, 0.57, f'{median_survival_positive}', transform=ax_to_plot.transAxes, fontsize=10, ha='center', va='center')
            ax_to_plot.text(1.1, 0.47, f'{median_survival_negative}', transform=ax_to_plot.transAxes, fontsize=10, ha='center', va='center')
            
            ax_to_plot.text(0.3, 0.17, "Hazard ratio,", transform=ax_to_plot.transAxes, fontsize=10, ha='left', va='center')
    
    # ax.annotate(f'p value: {p_value}', xy=(1, 0.5), xytext=(100, 80), textcoords='offset points', fontsize=12)
    if xlabel is None:
        ax.set_xlabel("Months since cfDNA collection")
    else: 
        ax.set_xlabel(xlabel)
    ax.set_ylabel("Survival fraction")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.legend(loc='best', frameon=False)
    ax.set_ylim((0, 1))
    ax.set_yticks((0, 0.2, 0.4, 0.6, 0.8, 1))
    ax.set_yticklabels(("0", "0.2", "0.4", "0.6", "0.8", "1"))
    if xmax is None:
        xmax = ax.get_xlim()[1]
    ax.set_xlim((0, xmax))
    ax.set_xticks(np.arange(0, xmax + 10, 10))
    if xmax >= 100:
        ax.set_xticks(np.arange(0, xmax + 10, 20))
    if xmax < 40:
        ax.set_xticks(np.arange(0, xmax + 10, 5))
    
    # add risk counts
    add_at_risk_counts(kmf_positive, kmf_negative, ax=ax, rows_to_show = ["At risk"])
    
    if not add_legend:
        ax.legend_.remove()
    ax.set_title(plot_title)
    
    # label the curves
    if neg_curve_label_positions is not None and pos_curve_label_positions is not None:
        ax.text(neg_curve_label_positions[0], 
                   neg_curve_label_positions[1], 
                   f'{label_keyword}-', transform=ax.transAxes, fontsize=10, ha='center', va='center', color = neg_curve_color)
        ax.text(pos_curve_label_positions[0], 
                   pos_curve_label_positions[1], 
                   f'{label_keyword}+', transform=ax.transAxes, fontsize=10, ha='center', va='center', color = pos_curve_color)
    
    # Run CPH here and print the hazard ratios to the plot.
    if print_HRs:
        try:
            ax = run_cox_proportional_hazards(surv_df_merged, duration_col, event_col, stratify_by, ax_to_plot)
        except Exception as e:
            print(f"An error occurred while running Cox Proportional Hazards.")
    
    if output_path is not None:
        fig.tight_layout()
        fig.savefig(output_path)
    return ax

def plot_vaf_scatter(df_chip, ax1, annotate_genes = True, fontsize = 10):
    """
    Plots a scatter plot comparing the cfDNA VAF to WBC VAF.
    """  
    # helps remove dependent calls with vaf 0 for both wbc and cfdna
    if "Dependent" in df_chip.columns:
        df_chip = df_chip[df_chip["Dependent"] == False]
    
    df_chip = df_chip.sort_values(by="VAF_n", ascending=False)
    df_chip["DTA genes"] = df_chip["Gene"].isin(["DNMT3A", "TET2", "ASXL1"])
    df_chip["zorder"] = np.random.randint(1, 3, df_chip.shape[0])
    df_chip["color"] = df_chip["DTA genes"].map({True: "blue", False: "black"})
    
    # Plotting
    sns.regplot(x="VAF_n", y="VAF_t", data=df_chip, ax=ax1, scatter = False, line_kws={'color': 'red', 'linewidth': 0.5, "linestyle": "--"})
    ax1.scatter(df_chip["VAF_n"], df_chip["VAF_t"], s = 1, c = df_chip["color"].tolist(), alpha = 0.6)
    
    # Run spearman's corr
    import scipy.stats as stats
    spearman_corr, p_value = stats.spearmanr(df_chip["VAF_n"], df_chip["VAF_t"])
    spearman_corr = round(spearman_corr, 2)
    p_value_str = f"{p_value:.2e}"
    ########################################################################################################
    # Add an inset ax to show the zoomed in values between 0-2.
    from mpl_toolkits.axes_grid1.inset_locator import inset_axes
    # inset_ax = inset_axes(ax1, width="35%", height="35%", loc=2, bbox_to_anchor=(0.15, 0.9, 0.35, 0.35), bbox_transform=ax1.transAxes)
    inset_ax = inset_axes(ax1, width="55%", height="55%", loc='upper left', bbox_to_anchor=(0.1, 0.5, 0.5, 0.5), bbox_transform=ax1.transAxes)
    
    # Move the inset plot to the right by adjusting the bbox_to_anchor
    inset_ax.set_position([inset_ax.get_position().x0 + 0.25,  # Shift right by increasing x0
                           inset_ax.get_position().y0,         # Keep the same y0
                           inset_ax.get_position().width, 
                           inset_ax.get_position().height])
    
    inset_ax.tick_params(axis='both', which='major', labelsize=5, width = 0.25, length = 1)  # Set tick size
    inset_ax.spines[['top', 'right']].set_visible(False)
    inset_ax.spines[['bottom', 'left']].set_linewidth(0.25)
    
    df_chip_2_perc = df_chip[(df_chip["VAF_t"] < 2) & (df_chip["VAF_n"] < 2)]
    sns.regplot(x="VAF_n", y="VAF_t", data=df_chip_2_perc, ax=inset_ax, scatter = False, line_kws={'color': 'red', 'linewidth': 0.5, "linestyle": "--"})
    inset_ax.scatter(df_chip_2_perc["VAF_n"], df_chip_2_perc["VAF_t"], s = 0.5, c = df_chip_2_perc["color"].tolist(), alpha = 0.6)
    inset_ax.set_ylim((0, 2))
    inset_ax.set_xlim((0, 2))
    inset_ax.set_yticks([0, 1, 2])
    inset_ax.set_xticks([0, 1, 2])
    inset_ax.set_yticklabels(["0", "1", "2"], fontsize = 6)
    inset_ax.set_xticklabels(["0", "1", "2"], fontsize = 6)
    inset_ax.set_ylabel("")
    inset_ax.set_xlabel("")
    inset_ax.xaxis.set_tick_params(pad=2)  # Reduces padding for x-ticks
    inset_ax.yaxis.set_tick_params(pad=2)  # Reduces padding for y-ticks
    ########################################################################################################
    
    # ax1.text(0.45, 0.85, f'Spearman\'s\ncorr: {spearman_corr:.2f}\np={p_value_str}', transform=ax1.transAxes, fontsize=8, color='red')
    ax1.spines["top"].set_visible(False)
    ax1.spines["right"].set_visible(False)
    ax1.set_xlabel("WBC VAF%", fontsize=fontsize)
    ax1.set_ylabel("Plasma cfDNA VAF%", fontsize=fontsize)
    
    if annotate_genes:
        top_points = df_chip.nlargest(3, 'VAF_t').reset_index(drop = True)
        texts = []
        line_length = 10
        for i, row in top_points.iterrows():
            y_pos_text = row["VAF_t"]-3
            if i < 2:
                # Top 2 points will have a line extending to the left
                x_right = row["VAF_n"]
                x_left = x_right-line_length
                ax1.text(x_left, y_pos_text, f'$\it{{{row["Gene"]}}}$', fontsize=8, color='black', ha = "right")
            else:
                x_left = row["VAF_n"]
                x_right = x_right+line_length
                ax1.text(x_right, y_pos_text, f'$\it{{{row["Gene"]}}}$', fontsize=8, color='black', ha = "left")
            ax1.plot([x_left, x_right], [row["VAF_t"], row["VAF_t"]], color='black', linewidth=0.5)
    
        # Couple of misc points to label, do it manually
        ax1.text(44.8+10, 34.3-3, f'$\it{"DNMT3A"}$', fontsize=8, color='blue', ha = "left")
        ax1.plot([44.8, 44.8+10], [34.3, 34.3], color='blue', linewidth=0.5)
        
        ax1.text(38.55-10, 41.48-3, f'$\it{"ASXL1"}$', fontsize=8, color='blue', ha = "right")
        ax1.plot([38.55, 38.55-10], [41.48, 41.48], color='blue', linewidth=0.5)
    
    ax1.set_ylim([0, 60])
    ax1.set_xlim([0, 60])
    ax1.set_yticks([0, 20, 40, 60])
    ax1.set_yticklabels(["0", "20", "40", "60"], fontsize = fontsize)
    ax1.set_xticks([0, 20, 40, 60])
    ax1.set_xticklabels(["0", "20", "40", "60"], fontsize = fontsize)
    
    # add legend for the colors
    legend_colors = ["blue", "black"]
    legend_labels = ["DTA", "non-DTA"]
    legend_handles = [plt.Line2D([0], [0], marker='o', color=color, label=label, markersize=2, linestyle='') for color, label in zip(legend_colors, legend_labels)]
    ax1.legend(handles=legend_handles, loc="lower right", frameon=False, fontsize = 8, handlelength=2, handletextpad=0.1)
    return(ax1)


def plot_per_patient_counts_grouped_bar_chart_UNGROUPED(muts_df, grouping_col_name, color, group1_ntotal, group2_ntotal, ax):
    """
    In a given ax plots the fraction of the bladder and kidney cohorts separately that have CH mutations.
    UNGROUPED - MEANING NO ARMS, NO DIAGNOSIS ETC USED FOR PRINCE.
    """
    # Extract names of groups
    group1_name = "LuPSMA"
    group2_name = "Cabazitaxel"
    
    # Determine colors for each group
    group1_color = color_dict[group1_name]
    group2_color = color_dict[group2_name]
    
    # get mut counts
    mutation_counts = muts_df.groupby([grouping_col_name, 'Patient_id']).size().reset_index(name='Mutation count')
    count_table = mutation_counts.groupby([grouping_col_name, 'Mutation count']).size().reset_index(name='Number of patients')
    pivot_table = count_table.pivot(index='Mutation count', columns=grouping_col_name, values='Number of patients').fillna(0).astype(int)
    pivot_table.columns = [f'Number of patients {col.lower()}' for col in pivot_table.columns]
    pivot_table.reset_index(inplace=True)
    pivot_table[f"{group1_name}_fraction"] = pivot_table[f"Number of patients {group1_name.lower()}"] / group1_ntotal
    pivot_table[f"{group2_name}_fraction"] = pivot_table[f"Number of patients {group2_name.lower()}"] / group2_ntotal
    
    # Get the number of patients with 0 mutations
    n_ch_positive_group1 = muts_df[muts_df[grouping_col_name] == group1_name]["Patient_id"].unique().shape[0]
    n_ch_positive_group2 = muts_df[muts_df[grouping_col_name] == group2_name]["Patient_id"].unique().shape[0]
    n_0_group1 = group1_ntotal - n_ch_positive_group1
    n_0_group2 = group2_ntotal - n_ch_positive_group2
    
    pts_with_0_muts_dict = {
        "Mutation count": 0, 
        f"Number of patients {group1_name.lower()}": n_0_group1, 
        f"Number of patients {group2_name.lower()}": n_0_group2,
        f"{group1_name}_fraction": n_0_group1/group1_ntotal,
        f"{group2_name}_fraction": n_0_group2/group2_ntotal}
    
    pts_with_0_muts_df = pd.DataFrame.from_dict(pts_with_0_muts_dict, orient='index').T
    pivot_table = pd.concat([pivot_table, pts_with_0_muts_df]).sort_values(by = "Mutation count").reset_index(drop = True)
    pivot_table["xpos"] = pivot_table.index
    # pivot_table.loc[pivot_table["xpos"] == 12, "xpos"] = 13
    
    # For the swarm plot get the vafs of all muts
    muts_df = muts_df.merge(mutation_counts, how = "left")
    group1_muts = muts_df[muts_df[grouping_col_name] == group1_name].reset_index(drop = True)
    group2_muts = muts_df[muts_df[grouping_col_name] == group2_name].reset_index(drop = True)
    
    # Plot them
    ax2 = ax.twinx()
    for i, row in pivot_table.iterrows():
        ax.bar(row["xpos"]-0.2, row[f"{group1_name}_fraction"], color=group1_color, width = 0.4, edgecolor = "None")
        ax.text(row["xpos"] - 0.2, 0.3, int(row[f"Number of patients {group1_name.lower()}"]), ha='center', va='center', fontsize=5, color='black')
        
        ax.bar(row["xpos"]+0.2, row[f"{group2_name}_fraction"], color=group2_color, width = 0.4, edgecolor = "None")
        ax.text(row["xpos"] + 0.2, 0.3, int(row[f"Number of patients {group2_name.lower()}"]), ha='center', va='center', fontsize=5, color='black')
    
    # plot the vafs in logscale
    muts_df["VAF_n_log"] = np.log10(muts_df["VAF_n"].replace(0, np.nan))
    for i, row in muts_df.iterrows():
        if row[grouping_col_name] == group1_name: 
            offset = 0.2
        else:
            offset = -0.2
        jitter = np.random.uniform(-0.08, 0.08, 1)
        ax2.scatter(row["Mutation count"]+offset+jitter[0], row["VAF_n_log"], color="black", s = 0.15, alpha = 0.7)
    
    # plot the vafs for bladder
    # for i, row in bladder_muts.iterrows():
    #     jitter = np.random.normal(-0.05, 0.05, 1)
    #     print(row["Mutation count"]-0.2+jitter, row["VAF_n"])
    #     ax2.scatter(row["Mutation count"]-0.2+jitter, row["VAF_n"], color="black", s = 2)
    
    # Aes
    for a in [ax, ax2]:
        a.spines["top"].set_visible(False)
    # x ticks
    ax.set_xticks(pivot_table["xpos"])
    ax.set_xticklabels(pivot_table["xpos"])
    ax.set_xlabel("Number of CH mutations")
    ax.set_ylabel("% of patients in arm")
    ax.tick_params(axis='x', bottom=False)
    ax2.tick_params(axis='x', pad=-5)
    # ax.tick_params(axis="both", direction="out", which="both", left=True, bottom=True , colors='k')    
    ax.set_ylim((0, 0.3))
    ax.set_yticks([0, 0.1, 0.2, 0.3])
    ax.set_yticklabels(["0", "10", "20", "30"])
    ax2.set_ylim((np.log10(0.25), np.log10(60)))
    ax2.set_yticks([np.log10(0.25), np.log10(1), np.log10(2), np.log10(10), np.log10(50)])
    ax2.set_yticklabels([".25", "1", "2", "10", "50"])
    ax2.set_ylabel("WBC VAF %")
    
    # Add legend
    legend_colors = color_dict.values()
    legend_labels = color_dict.keys()
    legend_handles = [plt.Line2D([0], [0], marker='s', color=color, label=label, markersize=5, linestyle='') for color, label in zip(legend_colors, legend_labels)]
    ax.legend(handles=legend_handles, loc="upper right", frameon=False, fontsize = 8, handlelength=2, handletextpad = 0.5)
    return([ax, ax2])


def calculate_OR_and_p(muts, gene, ntotal_pt_lu, n_total_pt_caba):
    """
    For a given gene calculates the OR of having a mutation in that gene in LuPSMA vs Cabazitaxel.
    """
    if gene is not None:
        muts_gene = muts[muts["Gene"] == gene]
    else: 
        muts_gene=muts.copy()
    
    # Number of patients mutated vs not in each arm
    npts_mutated_lu = muts_gene[muts_gene["Arm"] == "LuPSMA"]["Patient_id"].drop_duplicates().shape[0]
    npts_mutated_caba = muts_gene[muts_gene["Arm"] == "Cabazitaxel"]["Patient_id"].drop_duplicates().shape[0]      
    
    npts_nonmutated_lu=ntotal_pt_lu-npts_mutated_lu
    npts_nonmutated_caba=n_total_pt_caba-npts_mutated_caba    
    
    # Create the contingency table
    contingency_table = [[npts_mutated_lu, npts_nonmutated_lu], [npts_mutated_caba, npts_nonmutated_caba]]     
    
    # Calculate Odds Ratio (OR)
    OR = (npts_mutated_lu / npts_nonmutated_lu) / (npts_mutated_caba / npts_nonmutated_caba)
    _, p_value = stats.fisher_exact(contingency_table)
    result_dict = {"OR": OR, "p Fisher": p_value}
    
    return(result_dict)


