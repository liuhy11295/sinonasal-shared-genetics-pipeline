# GitHub publication checklist

Before pushing this repository:

- Confirm no `.slurm` scheduler templates are present.
- Confirm no raw GWAS, reference panels, result tables, logs, archives, or private
  server paths are present.
- Run `python tools/audit_github_code_export.py .`.
- Run `python -m pytest` if the Python test environment is available.
- Add a GitHub remote only after the safety audit passes.

Suggested first push:

```bash
git remote add origin https://github.com/<user>/<repo>.git
git branch -M main
git push -u origin main
```

