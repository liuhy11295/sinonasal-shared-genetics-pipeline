PROJECT_ROOT <- if (dir.exists(file.path(getwd(), "figure_final"))) {
  normalizePath(getwd(), mustWork = TRUE)
} else {
  normalizePath(file.path(getwd(), ".."), mustWork = TRUE)
}

source(file.path(PROJECT_ROOT, "figure_final", "scripts",
                 "build_main_figure_components.R"))
source(file.path(PROJECT_ROOT, "figure_final", "scripts",
                 "build_figure4_directional_panel.R"))

p4a_arc <- p4_candidate_arc_network +
  labs(
    title = "Genetically prioritized drug-target arc network",
    subtitle = paste(
      "Approved drugs connect to prioritized Gene-A/Gene-B target genes;",
      "line width reflects disease-pair support"
    )
  )

p4c_wrapped <- wrap_elements(full = directional_preview)

figure4_top <- (p4a_arc | p4b) +
  plot_layout(widths = c(1.28, 1.00), guides = "collect") &
  theme(legend.position = "bottom")

figure4_assembled <- figure4_top / p4c_wrapped +
  plot_layout(heights = c(1.00, 0.58), guides = "collect") +
  plot_annotation(tag_levels = "A") &
  theme(
    plot.tag = element_text(face = "bold", size = 10,
                            colour = COL["graphite"]),
    plot.margin = margin(4, 4, 4, 4)
  )

save_artifact(p4a_arc, "Figure4A",
              file.path("main", "panels"),
              183, 132, "Figure4", "A", "panel", dpi = 600)
save_artifact(directional_preview, "Figure4C",
              file.path("main", "panels"),
              183, 69, "Figure4", "C", "panel", dpi = 600)
save_artifact(figure4_assembled, "Figure4_revised",
              file.path("main", "composites"),
              183, 198, "Figure4", "ABC", "composite", dpi = 600)
copy_composite_to_root("Figure4_revised")

cat("Assembled Figure 4 panels: top row arc-network A + matrix B, bottom row C\n")
