# Result table contract

This project expects result tables to keep stable identifiers across methods.

## Required pair fields

- `pair_id`
- `trait1`
- `trait2`

## Required locus fields

- `pair_id`
- `trait1`
- `trait2`
- `locus_id`
- `CHR`
- `START`
- `END`

## Required SNP fields

- `pair_id`
- `trait1`
- `trait2`
- `SNP`
- `CHR`
- `BP`
- `EA`
- `OA`

## Genome-wide PLACO output

- `pair_id`
- `trait1`
- `trait2`
- `SNP`
- `CHR`
- `BP`
- `EA`
- `OA`
- `Z_trait1`
- `Z_trait2`
- `P_trait1`
- `P_trait2`
- `PLACO_stat`
- `PLACO_p`
- `PLACO_bh_q`
- `significance_class`

## Genome-wide CPASSOC output

- `pair_id`
- `trait1`
- `trait2`
- `SNP`
- `CHR`
- `BP`
- `EA`
- `OA`
- `P_trait1`
- `P_trait2`
- `SHet`
- `SHet_p`
- `SHet_bh_q`
- `SHom`
- `SHom_p`
- `SHom_bh_q`
- `significance_class`

Genome-wide PLACO/CPASSOC retained rows are restricted to `P < 5e-8` or
`BH-FDR < 0.05`.

