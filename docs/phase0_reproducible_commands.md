# Phase 0 Reproducible Commands

Working directory:

```powershell
Set-Location "D:\生信"
```

## 1. Re-check inventory and queue

```powershell
$inventory = Import-Csv .\gwas_inventory.tsv -Delimiter "`t"
$queue = Import-Csv .\gwas_download_queue.tsv -Delimiter "`t"
$inventory.Count
$queue | Group-Object first_pass,action | Select-Object Name,Count
```

Expected:

```text
11
yes, download_for_phase0: 7
yes, already_downloaded_provisional: 1
no, defer_download: 3
```

## 2. Create download directories

```powershell
New-Item -ItemType Directory -Force -Path data\gwas_raw,data\gwas_metadata,data\gwas_munged,results\phase0_shared_genetics | Out-Null
```

## 3. Download first-pass GWAS files

This intentionally uses `gwas_download_queue.tsv` so the queue controls what is downloaded.

```powershell
$queue = Import-Csv .\gwas_download_queue.tsv -Delimiter "`t"
$targets = $queue | Where-Object { $_.first_pass -eq "yes" -and $_.action -eq "download_for_phase0" }
foreach ($t in $targets) {
  $out = Join-Path "data\gwas_raw" ([IO.Path]::GetFileName($t.download_url))
  if (-not (Test-Path $out)) {
    Invoke-WebRequest -Uri $t.download_url -OutFile $out -TimeoutSec 0
  }
}
```

## 4. Verify downloaded sizes

```powershell
Get-ChildItem data\gwas_raw | Select-Object Name,Length,LastWriteTime | Sort-Object Name
```

## 5. Inspect headers before munging

```powershell
Get-ChildItem data\gwas_raw\*.gz | ForEach-Object {
  Write-Output "==== $($_.Name)"
  gzip -dc $_.FullName | Select-Object -First 1
}
Get-Content data\sources\gwas_catalog_metadata\GCST006358\IGE.V1_GWAS.txt -TotalCount 1
```

If `gzip` is unavailable on Windows, use R or Python gzip readers instead.

## 6. Build decision

Current build situation:

- FinnGen R12 disease and expansion rows are GRCh38.
- GWAS Catalog IgE / eosinophil / neutrophil / CRP rows are GRCh37.

Do not run LDSC across these files until the build strategy is fixed. The likely first-pass choice is to harmonize all inputs to GRCh37 for compatibility with common LDSC European reference assets.

