make_panel <- function(
    scale_key,
    method,
    caption_prefix   = "Concensus Model:",  # set to NULL to drop prefix
    legend_position  = NULL,                # e.g. "bottom" or NULL
    axes_collect_y   = FALSE,
    drop_guides      = FALSE                # TRUE for panel A only
) {
  # Extract plots from mc_efa_results
  cooccurrence_p <- mc_efa_results[[scale_key]][["cooccurrence"]][[method]][["plot_heatmap"]]
  stability_p    <- mc_efa_results[[scale_key]][["cooccurrence"]][[method]][["plot_stability"]]
  
  # Build caption from the existing caption in stability_p
  cap_raw       <- gsub("\n\n", "|", c(stability_p@labels$caption))
  cap_extracted <- stringr::word(cap_raw, 2, sep = "\\|")
  caption_text  <- if (is.null(caption_prefix)) {
    cap_extracted
  } else {
    paste(caption_prefix, cap_extracted)
  }
  
  # Start with base stability plot adjustments
  stab_layered <- stability_p &
    ggsci::scale_fill_jama() &
    ggplot2::labs(caption = caption_text)
  
  # Optionally drop guides (for first panel)
  if (drop_guides) {
    stab_layered <- stab_layered &
      guides(color = "none", fill = "none")
  }
  
  # Layout: optionally collect y-axes
  layout_args <- if (axes_collect_y) list(axes = "collect_y") else list()
  
  # Combine stability + cooccurrence with shared theming
  (stab_layered +
      (cooccurrence_p & ggplot2::labs(subtitle = NULL))) +
    do.call(patchwork::plot_layout, layout_args) &
    ggplot2::theme_minimal(base_size = 16) &
    ggplot2::theme(
      legend.position = legend_position,
      axis.text.x = ggplot2::element_text(
        angle = 45,
        vjust = .85,
        hjust = 1
      )
    )
}