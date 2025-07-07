def plot_vaf_scatter(df_chip, ax1, annotate_genes = True, fontsize = 10, vafcolname="VAF_n"):
    """
    Plots a scatter plot comparing the cfDNA VAF to WBC VAF.
    """  
    # helps remove dependent calls with vaf 0 for both wbc and cfdna
    if "Dependent" in df_chip.columns:
        df_chip = df_chip[df_chip["Dependent"] == False]
    
    df_chip = df_chip.sort_values(by=vafcolname, ascending=False)
    df_chip["DTA genes"] = df_chip["Gene"].isin(["DNMT3A", "TET2", "ASXL1"])
    df_chip["zorder"] = np.random.randint(1, 3, df_chip.shape[0])
    df_chip["color"] = df_chip["DTA genes"].map({True: "blue", False: "black"})
    
    # Plotting
    sns.regplot(x=vafcolname, y="VAF_t", data=df_chip, ax=ax1, scatter = False, line_kws={'color': 'red', 'linewidth': 0.5, "linestyle": "--"})
    ax1.scatter(df_chip[vafcolname], df_chip["VAF_t"], s = 1, c = df_chip["color"].tolist(), alpha = 0.6)
    
    # Run spearman's corr
    import scipy.stats as stats
    spearman_corr, p_value = stats.spearmanr(df_chip[vafcolname], df_chip["VAF_t"])
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
    
    df_chip_2_perc = df_chip[(df_chip["VAF_t"] < 2) & (df_chip[vafcolname] < 2)]
    sns.regplot(x=vafcolname, y="VAF_t", data=df_chip_2_perc, ax=inset_ax, scatter = False, line_kws={'color': 'red', 'linewidth': 0.5, "linestyle": "--"})
    inset_ax.scatter(df_chip_2_perc[vafcolname], df_chip_2_perc["VAF_t"], s = 0.5, c = df_chip_2_perc["color"].tolist(), alpha = 0.6)
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
                x_right = row[vafcolname]
                x_left = x_right-line_length
                ax1.text(x_left, y_pos_text, f'$\it{{{row["Gene"]}}}$', fontsize=8, color='black', ha = "right")
            else:
                x_left = row[vafcolname]
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


def plot_per_patient_counts_grouped_bar_chart_UNGROUPED(muts, ntotal_pts, ax, barcolor, vafcolname):
    """
    In a given ax plots the fraction of the bladder and kidney cohorts separately that have CH mutations.
    UNGROUPED - MEANING NO ARMS, NO DIAGNOSIS ETC USED FOR PRINCE.
    """
        
    # get mut counts
    mutation_counts = muts.groupby('Patient_id').size().reset_index(name='Mutation count')
    count_table = mutation_counts.groupby('Mutation count').size().reset_index(name='Number of patients')
    count_table["Fraction"]=(count_table["Number of patients"]/ntotal_pts)*100
    
    # Get the number of patients with 0 mutations
    n_ch_positive = muts["Patient_id"].unique().shape[0]
    n_ch_neg = ntotal_pts - n_ch_positive
    
    pts_with_0_muts_dict = {
        "Mutation count": 0, 
        f"Number of patients": n_ch_neg, 
        f"Fraction": n_ch_neg/ntotal_pts*100}
    
    pts_with_0_muts_df = pd.DataFrame.from_dict(pts_with_0_muts_dict, orient='index').T
    count_table = pd.concat([count_table, pts_with_0_muts_df]).sort_values(by = "Mutation count").reset_index(drop = True)
    # count_table["xpos"] = count_table.index
    
    # For the swarm plot get the vafs of all muts
    muts = muts.merge(mutation_counts, how = "left")
    
    # Plot them
    ax2 = ax.twinx()
    for i, row in count_table.iterrows():
        row["Fraction"]
        ax.bar(row["Mutation count"], row["Fraction"], color=barcolor, width = 0.4, edgecolor = "None")
        ax.text(row["Mutation count"], 39, int(row["Number of patients"]), ha='center', va='center', fontsize=5, color='black')
        
    # plot the vafs in logscale
    muts["VAF_n_log"] = np.log10(muts[vafcolname].replace(0, np.nan))
    for i, row in muts.iterrows():
        jitter = np.random.uniform(-0.08, 0.08, 1)
        ax2.scatter(row["Mutation count"]+jitter[0], row["VAF_n_log"], color="black", s = 0.15, alpha = 0.7)
    
    # Aes
    for a in [ax, ax2]:
        a.spines["top"].set_visible(False)
    # x ticks
    ax.set_xticks(range(0, int(count_table["Mutation count"].max()+1)))
    ax.set_xticklabels(range(0, int(count_table["Mutation count"].max()+1)))
    ax.set_xlabel("Number of CH mutations")
    ax.set_ylabel("% of samples")
    ax.tick_params(axis='x', bottom=False)
    ax2.tick_params(axis='x', pad=-5)
    # ax.tick_params(axis="both", direction="out", which="both", left=True, bottom=True , colors='k')    
    ax.set_ylim((0, 40))
    ax.set_yticks([0, 10, 20, 30, 40])
    ax.set_yticklabels(["0", "10", "20", "30", "40"])
    ax2.set_ylim((np.log10(0.25), np.log10(60)))
    ax2.set_yticks([np.log10(0.25), np.log10(1), np.log10(2), np.log10(10), np.log10(50)])
    ax2.set_yticklabels([".25", "1", "2", "10", "50"])
    ax2.set_ylabel("WBC VAF %")
    
    return([ax, ax2])

def plot_ch_presence_absence_bars_UNGROUPED(muts, ax, ntotal_pts, vafcolname, barcolor):
    """
    Plots the presence / absence of CH in the cohort as bar charts. Separates into two groups.
    """
        
    results_dict = {}
    for min_vaf in [0, 2, 10]:
        muts_df_filtered  = muts[muts[vafcolname] >= min_vaf]
        
        pts = muts_df_filtered["Patient_id"].unique().shape[0]
        perc = round((pts/ntotal_pts)*100)
        
        # save results
        results_dict[min_vaf] = {"perc": perc, "pts": pts}
    
    # Now plotting the bar charts    
    df = pd.DataFrame.from_dict(results_dict, orient='index').reset_index().rename(columns = {"index": "min_vaf"})
    for i, row in df.iterrows():
        ax.bar(i, row["perc"], color = barcolor, edgecolor = None, width = 0.7)
        ax.text(i, row["perc"]+1, str(row["pts"]), ha='center', va='bottom', fontsize=7, color='black')
                
    ax.set_xlim((-0.7, 2.5))
    ax.set_xticks([0, 1, 2])
    ax.tick_params(axis='x', bottom=False)
    ax.tick_params(axis='x', pad=2)
    ax.set_xticklabels(["≥0.25%", "≥2%", "≥10%"], fontsize = 8)
    ax.set_ylabel("% of samples")
    ax.set_xlabel("Minimum CH VAF%")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    
    ax.set_ylim((0, 100))
    ax.set_yticks([0, 25, 50, 75, 100])
    ax.set_yticklabels(["0", "25", "50", "75", "100"])
    return(ax)

def plot_gene_counts_UNGROUPED_single_ax(muts, ax, ntotal_pts, gene_list = None, use_entire_denom = True, barcolor="gray"):
    """
    Plots the number of mutations each individual gene has in the cohort as a BAR CHART. If unique, non unique mutations in the same gene are counted only once.
    """
        
    if gene_list is not None: 
        df = muts[muts["Gene"].isin(gene_list)].reset_index(drop = True)
    
    # subset
    df=df[["Patient_id", "Gene"]].drop_duplicates()["Gene"].value_counts().reset_index().rename(columns = {"index": "Gene", "Gene": "counts"}).reset_index()
    if use_entire_denom:
        df["frac"]=df["counts"]/ntotal_pts*100
    else:
        n_ch_pos_pts=muts["Patient_id"].drop_duplicates().shape[0]
        df["frac"]=df["counts"]/ntotal_pts*100
    
    # plotting
    for i, row in df.iterrows():
        # Plot bars
            ax.barh(row["index"], row[f"frac"], color=barcolor, height=0.4, edgecolor="None")
            ax.text(row[f"frac"]+1, row["index"], row["counts"], va='center', ha="right", color="black", fontsize=7)
    
    # aesthetics
    ax.tick_params(axis='y', left=False, labelleft=False)
    ax.spines[["top", "left"]].set_visible(False)
    tick_df = df[["index", "Gene"]].sort_values(by = "index")
    ax.set_yticks(tick_df["index"].tolist())
    ax.set_yticklabels(tick_df["Gene"].tolist(), rotation = 0, fontstyle = "italic", va = "center", fontsize = 8, ha = "center")
    ax.tick_params(axis='y', which='both', pad=20)  # Adjust padding between tick labels and axis
    ax.tick_params(axis = "y", which='both', length=0)
    ax.set_xlabel("% mutated in CH+ pts")
    ax.set_xlim((0, 80))
    ax.set_xticks([0, 20, 40, 60, 80])
    ax.set_xticklabels(["0", "20", "40", "60", "80"])
    ax.set_ylim((-0.6, 14.6))
        
    # Generate the gene_order dict for the subsequent plot
    gene_order = df[["Gene", "index"]]
    ax.invert_xaxis()
    ax.set_xlabel("% mutated in CH+ pts")
    
    return(ax, gene_order)

def plot_vafs_swarm_UNGROUPED(muts_df, ax_swarm, vafcolname, gene_list = None, gene_order = None):
    """
    Plots the VAFs as swarm.
    """    
    # Subset to genes of interest
    if gene_list is not None: 
        df = muts_df[muts_df["Gene"].isin(gene_list)].reset_index(drop = True)
    
    df = df[["Gene", vafcolname]]
    
    if gene_order is not None:
        df = df.merge(gene_order, how = "left")
    
    df["VAF_n_log"] = np.log10(df[vafcolname].replace(0, np.nan))  # Replace 0 with NaN to avoid log(0)
            
    for i, row in df.iterrows():
        jitter = np.random.uniform(-0.05, 0.05, 1)
        ax_swarm.scatter(row["VAF_n_log"], row["index"]+jitter, color="black", alpha = 0.7, s = 0.15, zorder = 1) # group 1
        
    # ax_swarm.set_yscale('log')
    ax_swarm.spines[["right", "top"]].set_visible(False)
    ax_swarm.set_xlim((np.log10(0.25), np.log10(50)))
    ax_swarm.set_xticks([np.log10(0.25), np.log10(1), np.log10(2), np.log10(10), np.log10(50)])
    ax_swarm.set_xticklabels([".25", "1", "2", "10", "50"])
    ax_swarm.tick_params(axis='y', labelleft=False, left = False)
    ax_swarm.set_xlabel("WBC VAF%")
    
    return(ax_swarm)

def plot_mutcounts_per_gene(muts, ntotal_pts, gene_list, ax, color_dict):
    """
    In a given muts list plots the number of patients having n mutations
    """
    muts_subset = muts[muts["Gene"].isin(gene_list)]
    muts_count = muts_subset[["Patient_id", "Gene"]].value_counts().reset_index()
    muts_count.columns = ["Patient_id", "Gene", "Counts"]
    muts_count.loc[muts_count["Counts"] > 1, "Counts"] = 2 # More than 1 gets 2
    
    muts_count = muts_count[["Gene", "Counts"]].value_counts().reset_index()
    muts_count.columns = ["Gene", "Count_per_pt", "pt_count"]
    muts_count["barcolor"]=muts_count["Count_per_pt"].map(color_dict)
    
    # Plotting    
    for i, gene in enumerate(gene_list):
        gene_df = muts_count[muts_count["Gene"] == gene].sort_values(by = "Count_per_pt")
        gene_df["Fraction"] = (gene_df["pt_count"]/ntotal_pts)*100
        
        # Initialize y val to make stacked bar chart
        y_val_sum = 0
        for j, row in gene_df.iterrows():
            ax.bar(i, row["Fraction"], bottom=y_val_sum, color = row["barcolor"], edgecolor = None, width = 0.55)
            
            # Annotate patient numbers
            ax.text(i, row["Fraction"]+y_val_sum-2, str(row["pt_count"]), ha='center', va='center', fontsize=7, color='black')
            y_val_sum += row["Fraction"]
    
    # AES
    # ax.set_ylim((0, 100))
    # ax.set_yticks([0, 20, 40, 60, 80, 100])
    # ax.set_yticklabels(["0", "20", "40", "60", "80", "100"], fontsize=10)
    # ax.set_xlabel('Arm')
    ax.set_ylabel('% of samples')
    ax.set_xticks(range(0, len(gene_list)))
    ax.set_xticklabels(gene_list, rotation = 90, fontstyle="italic")
    ax.tick_params(axis='x', bottom=False, pad = 1)
    ax.spines[["top", "right"]].set_visible(False)
    
    # Add legend
    legend_colors = color_dict.values()
    legend_labels = ["1 mutation", ">1 mutation"]
    legend_handles = [plt.Line2D([0], [0], marker='s', color=color, label=label, markersize=5, linestyle='') for color, label in zip(legend_colors, legend_labels)]
    ax.legend(handles=legend_handles, loc="upper right", ncol = 1, frameon=False, fontsize = 8, handlelength=2, handletextpad = 0.1)
    
    return(ax)

