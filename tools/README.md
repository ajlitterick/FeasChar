# Tools

*Written by an AI system (Claude, Anthropic), October 2026, under the author's direction; see
"Who wrote what" in the main README.*

## Checking

| File | Purpose |
|---|---|
| `verify_row.g` | Checks any row of the updated tables from its witness, without the search code. See the main README. |
| `fg_dump.g`, `check_classdata.py` | Compare this code's semisimple classes (traces on the adjoint and minimal modules, all orders 2–17) with the legacy Magma data `legacy/n.M`. |
| `check_dims.py` | Checks that every row of the Memoir's tables, as parsed, has the right module dimensions. This is how the typographical errors listed in `tables/DIFFERENCES.md` were found. |

## Regenerating the updated tables

The Memoir's tables were turned into data in four steps:
1. `parse_tables_gen.py` reads the F₄, E₆ and E₇ tables, and `parse_tables_e8.py` the E₈ ones.
2. Both write JSON (`tables_gen.json`, `e8tables.json`).
3. `memoir_numbers.py` numbers every table as LaTeX does (6.1 – 6.352).
4. `build_index.py` combines these into `memoir_index.g`, the Memoir's rows with their numbers.

To regenerate, copy `feasgen.g`, `lubeck_data.g`, `memoir_index.g` and `fg_tables.g` into one
folder, make a subfolder `out/tables`, and from that folder run, for each table index i from 1 to
352:

```
printf 'IDX := i;\nRead("fg_tables.g");\n' > run.g && gap -q -o 8g run.g
```

`fg_tables.sh` does this in parallel. Each run writes `out/tables/6.i.tex`, `.txt` and `.json`.

Then, with `memoir_numbers.json` in the same folder, run

```
python3 assemble_tables.py memoir_numbers.json out/tables <dir>
```

This writes `tables.tex`, `F4.tex` … `E8.tex`, `F4.json` … `E8.json`, `F4.g` … `E8.g` and the
per-table part of `DIFFERENCES.md` into `<dir>`.

Most tables take seconds. A few need much more time and memory, because their groups have
elements of large order: 2·Ru and SL₂(29) in E₇ (order 58), 3·J₃ in E₆ (order 57), and L₂(61)
in E₈ (order 61). These need the semisimple classes of G of that order, about 1.5 million for
E₇ and E₈. Allow about 20 GB.

## Data

| File | Contents |
|---|---|
| `lubeck2gap.py` | Converts F. Lübeck's web tables of weight multiplicities into `../lubeck_data.g`. |
| `memoir_index.g`, `memoir_numbers.json` | The Memoir's tables as data, numbered as in the Memoir. |
