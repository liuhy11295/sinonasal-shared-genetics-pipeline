options(stringsAsFactors = FALSE)

find_figure_root <- function() {
  candidates <- c(getwd(), file.path(getwd(), "figure_final"),
                  file.path(getwd(), ".."))
  for (candidate in candidates) {
    candidate <- normalizePath(candidate, mustWork = FALSE)
    if (dir.exists(file.path(candidate, "main", "composites")) &&
        dir.exists(file.path(candidate, "source_data"))) {
      return(candidate)
    }
  }
  stop("Cannot locate figure_final project root")
}

FIG <- find_figure_root()
PDF <- file.path(FIG, "main", "composites", "Figure1.pdf")
PNG <- file.path(FIG, "main", "composites", "Figure1.png")
SRC <- file.path(FIG, "source_data", "Figure1_key_numbers.tsv")
OUT <- file.path(FIG, "qa", "Figure1_graphical_abstract_QA.tsv")
dir.create(dirname(OUT), recursive = TRUE, showWarnings = FALSE)

run <- function(command, args) {
  paste(system2(command, args, stdout = TRUE, stderr = TRUE), collapse = "\n")
}

pdf_info <- run("pdfinfo", shQuote(PDF))
pdf_fonts <- run("pdffonts", shQuote(PDF))
pdf_text <- run("pdftotext", c(shQuote(PDF), "-"))
png_info <- run("identify", c("-verbose", shQuote(PNG)))
src <- read.delim(SRC, sep = "\t", quote = "", check.names = FALSE,
                  comment.char = "", na.strings = "")

check <- function(name, pass, observed, expected) {
  data.frame(
    check = name,
    status = if (isTRUE(pass)) "pass" else "fail",
    observed = observed,
    expected = expected,
    stringsAsFactors = FALSE
  )
}

required_text <- c(
  "GWAS summary", "Genetic comorbidity", "Shared SNP and",
  "Gene and pathway", "Downstream", "31", "13", "18",
  "cross-trait SNP and locus evidence",
  "exploratory drug-target annotation",
  "secondary directional evidence",
  "Fig. 2", "Fig. 3", "Fig. 4"
)
forbidden_text <- c("FUMA", "CellChat", "PWAS", "single-cell")

fixed <- src[src$status == "final", ]
unfixed <- src[src$status == "intentionally_unnumbered", ]
fixed_expected <- c(
  nasal_phenotypes = 4,
  non_nasal_comparators = 57,
  positive_disease_pairs = 31,
  allergic_rhinitis_pairs = 13,
  nasal_polyps_pairs = 18
)
fixed_observed <- setNames(as.integer(fixed$value), fixed$metric)

qa <- rbind(
  check("pdf_exists", file.exists(PDF) && file.info(PDF)$size > 10000,
        if (file.exists(PDF)) file.info(PDF)$size else 0,
        "file size > 10000 bytes"),
  check("png_exists", file.exists(PNG) && file.info(PNG)$size > 10000,
        if (file.exists(PNG)) file.info(PNG)$size else 0,
        "file size > 10000 bytes"),
  check("pdf_single_page", grepl("Pages:[[:space:]]+1", pdf_info),
        sub(".*(Pages:[[:space:]]+[0-9]+).*", "\\1", pdf_info),
        "Pages: 1"),
  check("pdf_unencrypted", grepl("Encrypted:[[:space:]]+no", pdf_info),
        sub(".*(Encrypted:[[:space:]]+[^\\n]+).*", "\\1", pdf_info),
        "Encrypted: no"),
  check("pdf_vector_text", nchar(pdf_text) > 500,
        paste0(nchar(pdf_text), " extracted characters"),
        "> 500 extracted characters"),
  check("pdf_embedded_arial", grepl("Arial", pdf_fonts),
        if (grepl("Arial", pdf_fonts)) "Arial embedded" else "Arial absent",
        "Arial embedded"),
  check("png_dimensions", grepl("Geometry: 9600x3960", png_info),
        sub(".*(Geometry: [^\\n]+).*", "\\1", png_info),
        "Geometry: 9600x3960"),
  check("png_600_dpi", grepl("Resolution: 236\\.22x236\\.22", png_info),
        sub(".*(Resolution: [^\\n]+).*", "\\1", png_info),
        "236.22 px/cm (600 dpi)"),
  check("required_text_present",
        all(vapply(required_text, grepl, logical(1), x = pdf_text,
                   fixed = TRUE)),
        paste(required_text[vapply(required_text, grepl, logical(1),
                                   x = pdf_text, fixed = TRUE)],
              collapse = "; "),
        paste(required_text, collapse = "; ")),
  check("forbidden_text_absent",
        !any(vapply(forbidden_text, grepl, logical(1), x = pdf_text,
                    fixed = TRUE)),
        paste(forbidden_text[vapply(forbidden_text, grepl, logical(1),
                                    x = pdf_text, fixed = TRUE)],
              collapse = "; "),
        "none"),
  check("fixed_numbers_match",
        all(fixed_observed[names(fixed_expected)] == fixed_expected),
        paste(names(fixed_observed), fixed_observed, sep = "=",
              collapse = "; "),
        paste(names(fixed_expected), fixed_expected, sep = "=",
              collapse = "; ")),
  check("unfinalized_counts_omitted",
        all(is.na(unfixed$value)),
        paste(unfixed$metric, ifelse(is.na(unfixed$value), "blank",
                                    unfixed$value),
              sep = "=", collapse = "; "),
        "all intentionally unnumbered values blank"),
  check("pair_subtotals_consistent",
        fixed_observed["allergic_rhinitis_pairs"] +
          fixed_observed["nasal_polyps_pairs"] ==
          fixed_observed["positive_disease_pairs"],
        paste0(fixed_observed["allergic_rhinitis_pairs"], " + ",
               fixed_observed["nasal_polyps_pairs"], " = ",
               fixed_observed["positive_disease_pairs"]),
        "13 + 18 = 31")
)

write.table(qa, OUT, sep = "\t", quote = FALSE, row.names = FALSE, na = "")
if (any(qa$status != "pass")) {
  print(qa[qa$status != "pass", ], row.names = FALSE)
  stop("Figure 1 QA failed")
}
message("All ", nrow(qa), " Figure 1 QA checks passed")
