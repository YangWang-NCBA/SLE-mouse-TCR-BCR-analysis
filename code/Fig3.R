make_fig3 <- function(out_root) {
  folder <- file.path(out_root, "fig", "Fig3_umi_support")
  source_dir <- file.path(folder, "source_data")
  qc <- read_sd(file.path(DATA, "analysis/mixcr_library_qc.tsv"))
  qc <- ensure_group_metadata(qc, read_sd(file.path(DATA,"metadata/libraries.tsv")))
  stopifnot(!anyDuplicated(qc$analysis_library_id),
            all(is.finite(qc$umi_output_diversity)), all(qc$umi_output_diversity > 0))
  tissues <- c("blood","spleen","lymph_node","thymus","brain","kidney","liver","lung","heart","intestine","skin")
  tissue_labels <- c("Blood","Spleen","Lymph\nnodes","Thymus","Brain","Kidney","Liver","Lung","Heart","Intestine","Skin")
  qc$group_public <- factor_groups(qc$group_public)
  qc$receptor <- factor(qc$receptor,levels=c("BCR","TCR"))
  qc$tissue <- factor(qc$tissue,levels=tissues)
  qc <- qc[order(qc$receptor,qc$tissue,qc$group_public,qc$analysis_library_id),]
  write_source(qc,file.path(source_dir,"Fig3_umi_output_diversity.tsv"))
  # Calculate boxes from untransformed observed counts, then display them
  # on a log10 axis. Whiskers end at observed values within 1.5 raw-count IQRs.
  summary <- do.call(rbind,lapply(split(qc,interaction(qc$receptor,qc$tissue,qc$group_public,drop=TRUE)),function(z) {
    y <- z$umi_output_diversity
    q <- unname(quantile(y,c(.25,.5,.75),type=7));iq <- q[3]-q[1]
    data.frame(receptor=z$receptor[1],tissue=z$tissue[1],group_public=z$group_public[1],
      n_libraries=nrow(z),n_mice=length(unique(z$mouse_id)),
      lower=q[1],middle=q[2],upper=q[3],
      ymin=min(y[y>=q[1]-1.5*iq]),ymax=max(y[y<=q[3]+1.5*iq]),
      total_umi_output=sum(y))
  }))
  write_source(summary,file.path(source_dir,"Fig3_umi_output_diversity_summary.tsv"))
  p <- ggplot() +
    geom_boxplot(data=summary,aes(tissue,group=interaction(tissue,group_public),
                 ymin=ymin,lower=lower,middle=middle,upper=upper,ymax=ymax,
                 fill=group_public,colour=group_public),stat="identity",width=0.56,
                 position=position_dodge(width=0.66),alpha=0.24,linewidth=yw_lw(0.75)) +
    geom_point(data=qc,aes(tissue,umi_output_diversity,colour=group_public,group=group_public),
                 shape=16,size=1.05,alpha=0.8,
                 position=position_jitterdodge(jitter.width=0.11,jitter.height=0,dodge.width=0.66,seed=20261012),
                 show.legend=FALSE) +
    facet_wrap(~receptor,ncol=1) +
    scale_x_discrete(limits=tissues,labels=tissue_labels,drop=FALSE) +
    scale_y_log10(breaks=10^(3:7),labels=count_axis_labels,limits=c(900,1e7),expand=expansion(mult=c(0,.03))) +
    scale_colour_manual(values=GROUP_COLORS,breaks=GROUP_ORDER,labels=GROUP_LABELS,name=NULL) +
    scale_fill_manual(values=GROUP_COLORS,breaks=GROUP_ORDER,labels=GROUP_LABELS,name=NULL) +
    guides(colour=guide_legend(nrow=1),fill=guide_legend(nrow=1)) +
    labs(x=NULL,y="UMI output diversity") + theme_sd(base_size=7) +
    theme(legend.position="top",legend.justification="right",legend.direction="horizontal",
          legend.key.width=grid::unit(3,"mm"),legend.key.height=grid::unit(3,"mm"),
          legend.margin=margin(0,0,1,0),legend.box.spacing=grid::unit(1,"mm"),
          panel.spacing=grid::unit(5,"mm"),plot.margin=margin(3,3,3,3))
  note <- paste0("Fig. 3 | Output UMI diversity by specimen source. Boxplots show the UMI output diversity recorded in MiXCR refineTagsAndSort reports for BCR (upper row) and TCR (lower row) libraries, grouped by ten endpoint tissues and peripheral blood. Control and IMQ are shown separately within each source. Each filled circle represents one specimen-library. Each endpoint tissue includes five specimens per group; blood includes ten control and fifteen IMQ collections from five mice per group. Boxes show the median and 25th–75th percentiles; whiskers end at the most extreme observations within 1.5 times the interquartile range. Summaries are calculated from original counts and displayed on a log10 axis.")
  export_panel(p,folder,"Fig3_umi_output_diversity",note)
  writeLines(note,file.path(folder,"Fig3_panel_legends.txt"),useBytes=TRUE)
}
