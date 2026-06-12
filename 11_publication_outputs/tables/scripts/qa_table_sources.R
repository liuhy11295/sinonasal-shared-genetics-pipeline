options(stringsAsFactors = FALSE)

find_root <- function() {
  configured <- Sys.getenv("PAPER_ROOT")
  if (nzchar(configured)) return(normalizePath(configured, mustWork = TRUE))
  candidates <- c(getwd(), dirname(getwd()))
  for (candidate in candidates) {
    if (dir.exists(file.path(candidate, "Table")) &&
        dir.exists(file.path(candidate, "figure_final"))) {
      return(normalizePath(candidate))
    }
  }
  stop("Cannot locate project root")
}

root <- find_root()
manifest_path <- file.path(root, "Table", "source_data",
                           "table_source_manifest.tsv")
manifest <- read.delim(manifest_path, sep = "\t", quote = "",
                       check.names = FALSE)

resolve_source <- function(path) {
  if (grepl("^/", path)) path else file.path(root, path)
}

inspect_source <- function(path) {
  resolved <- resolve_source(path)
  exists <- file.exists(resolved) || dir.exists(resolved)
  kind <- if (dir.exists(resolved)) "directory" else "file"
  rows <- NA_integer_
  columns <- NA_integer_

  if (file.exists(resolved) && grepl("\\.(tsv|csv)$", resolved,
                                     ignore.case = TRUE)) {
    separator <- if (grepl("\\.csv$", resolved, ignore.case = TRUE)) "," else "\t"
    con <- file(resolved, open = "r")
    header <- readLines(con, n = 1, warn = FALSE)
    close(con)
    columns <- if (length(header)) {
      length(strsplit(header, separator, fixed = TRUE)[[1]])
    } else {
      0L
    }
    rows <- max(0L, length(readLines(resolved, warn = FALSE)) - 1L)
  }

  data.frame(
    source = path,
    resolved_source = resolved,
    source_type = kind,
    exists = exists,
    data_rows = rows,
    data_columns = columns
  )
}

inventory <- do.call(rbind, lapply(unique(manifest$authoritative_source),
                                  inspect_source))
inventory <- merge(manifest, inventory,
                   by.x = "authoritative_source", by.y = "source",
                   all.x = TRUE, sort = FALSE)

out <- file.path(root, "Table", "qa", "table_source_inventory.tsv")
write.table(inventory, out, sep = "\t", quote = FALSE, row.names = FALSE,
            na = "")

if (any(!inventory$exists)) {
  missing <- unique(inventory$authoritative_source[!inventory$exists])
  stop("Missing table sources: ", paste(missing, collapse = ", "))
}

cat("Table source QA passed:", nrow(inventory), "manifest entries checked\n")
