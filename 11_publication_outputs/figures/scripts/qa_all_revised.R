#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE)

PROJECT_ROOT <- if (dir.exists(file.path(getwd(), "figure_final"))) {
  normalizePath(getwd(), mustWork = TRUE)
} else {
  normalizePath(file.path(getwd(), ".."), mustWork = TRUE)
}
FIG <- file.path(PROJECT_ROOT, "figure_final")

png_files <- sort(list.files(
  c(file.path(FIG, "main"), file.path(FIG, "supplementary")),
  pattern = "\\.png$", recursive = TRUE, full.names = TRUE
))
pdf_files <- sort(list.files(
  c(file.path(FIG, "main"), file.path(FIG, "supplementary")),
  pattern = "\\.pdf$", recursive = TRUE, full.names = TRUE
))

identify_resolution <- function(path) {
  out <- system2(
    "identify",
    c("-format", shQuote("%x|%y|%U"), shQuote(path)),
    stdout = TRUE, stderr = TRUE
  )
  parts <- strsplit(out[1], "|", fixed = TRUE)[[1]]
  x <- suppressWarnings(as.numeric(parts[1]))
  y <- suppressWarnings(as.numeric(parts[2]))
  unit <- parts[3]
  dpi_x <- if (identical(unit, "PixelsPerCentimeter")) x * 2.54 else x
  dpi_y <- if (identical(unit, "PixelsPerCentimeter")) y * 2.54 else y
  c(dpi_x = dpi_x, dpi_y = dpi_y)
}

pdf_text_chars <- function(path) {
  out <- system2("pdftotext", c(shQuote(path), "-"),
                 stdout = TRUE, stderr = FALSE)
  nchar(gsub("\\s+", "", paste(out, collapse = "")))
}

pdf_info_value <- function(path, field) {
  out <- system2("pdfinfo", shQuote(path), stdout = TRUE, stderr = FALSE)
  hit <- out[grepl(paste0("^", field, ":"), out)]
  if (!length(hit)) return(NA_character_)
  trimws(sub("^[^:]+:", "", hit[1]))
}

png_qa <- do.call(rbind, lapply(png_files, function(path) {
  resolution <- identify_resolution(path)
  data.frame(
    artifact = sub(paste0("^", FIG, "/"), "", path),
    format = "png",
    file_size = file.info(path)$size,
    dpi_x = round(resolution["dpi_x"], 2),
    dpi_y = round(resolution["dpi_y"], 2),
    text_chars = NA_integer_,
    page_count = NA_integer_,
    encrypted = NA_character_,
    page_structure_ok = NA,
    status = ifelse(
      file.info(path)$size > 15000 &&
        resolution["dpi_x"] >= 599 &&
        resolution["dpi_y"] >= 599,
      "pass", "review"
    )
  )
}))

pdf_qa <- do.call(rbind, lapply(pdf_files, function(path) {
  text_chars <- pdf_text_chars(path)
  page_count <- suppressWarnings(as.integer(pdf_info_value(path, "Pages")))
  encrypted <- tolower(pdf_info_value(path, "Encrypted"))
  page_structure_ok <- page_count == 1 ||
    (basename(path) == "Supplementary_all_pair_circular_loci.pdf" &&
       page_count > 1)
  data.frame(
    artifact = sub(paste0("^", FIG, "/"), "", path),
    format = "pdf",
    file_size = file.info(path)$size,
    dpi_x = NA_real_,
    dpi_y = NA_real_,
    text_chars = text_chars,
    page_count = page_count,
    encrypted = encrypted,
    page_structure_ok = page_structure_ok,
    status = ifelse(file.info(path)$size > 10000 && text_chars > 0 &&
                      page_structure_ok && encrypted == "no",
                    "pass", "review")
  )
}))

qa <- rbind(png_qa, pdf_qa)
write.table(
  qa, file.path(FIG, "qa", "final_figure_QA.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE, na = ""
)

source_files <- list.files(file.path(FIG, "source_data"),
                           pattern = "\\.tsv$", full.names = FALSE)
supplementary_content <- do.call(rbind, lapply(seq_len(18), function(i) {
  stem <- paste0("SupplementaryFigure", i, "_revised")
  source_pattern <- if (i == 4) {
    "^Supplementary_Figure4_.*\\.tsv$"
  } else {
    paste0("^SupplementaryFigure", i, "_.*\\.tsv$")
  }
  data.frame(
    supplementary_figure = i,
    pdf_present = file.exists(file.path(
      FIG, "supplementary", "composites", paste0(stem, ".pdf")
    )),
    png_present = file.exists(file.path(
      FIG, "supplementary", "composites", paste0(stem, ".png")
    )),
    source_tables = sum(grepl(source_pattern, source_files)),
    stringsAsFactors = FALSE
  )
}))
supplementary_content$status <- ifelse(
  supplementary_content$pdf_present &
    supplementary_content$png_present &
    supplementary_content$source_tables > 0,
  "pass", "review"
)
write.table(
  supplementary_content,
  file.path(FIG, "qa", "supplementary_content_QA.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

summary <- data.frame(
  check = c("png_files", "pdf_files", "png_600dpi", "pdf_editable_text",
            "pdf_page_structure_unencrypted",
            "supplementary_figures_complete",
            "all_files_pass"),
  observed = c(
    nrow(png_qa),
    nrow(pdf_qa),
    sum(png_qa$status == "pass"),
    sum(pdf_qa$text_chars > 0),
    sum(pdf_qa$page_structure_ok & pdf_qa$encrypted == "no"),
    sum(supplementary_content$status == "pass"),
    sum(qa$status == "pass")
  ),
  expected = c(
    nrow(png_qa),
    nrow(pdf_qa),
    nrow(png_qa),
    nrow(pdf_qa),
    nrow(pdf_qa),
    18,
    nrow(qa)
  )
)
write.table(
  summary, file.path(FIG, "qa", "final_figure_QA_summary.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

if (any(qa$status != "pass") ||
    any(supplementary_content$status != "pass")) {
  print(qa[qa$status != "pass", ])
  print(supplementary_content[supplementary_content$status != "pass", ])
  stop("Figure QA failed", call. = FALSE)
}

print(summary)
