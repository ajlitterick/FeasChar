#!/bin/bash
# Usage: fg_tables.sh <timeout-s> <parallel> idx1 idx2 ...   (memoir_index.g indices, 1..352)
# Writes out/tables/<6.N>.tex and .txt; GAPMEM (default 8g) sets GAP's memory limit.
# run from a folder containing feasgen.g, lubeck_data.g, memoir_index.g and fg_tables.g
TO=$1; PAR=$2; shift 2
mkdir -p out/tables out/tables/run
GAP=${GAP:-gap}
export TO GAP
one() {
  i=$1
  printf 'IDX := %d;\nRead("fg_tables.g");\n' "$i" > "out/tables/run/run_$i.g"
  timeout "$TO" $GAP -q -o ${GAPMEM:-8g} "out/tables/run/run_$i.g" < /dev/null > "out/tables/run/raw_$i.txt" 2>&1
  if [ "$?" = "124" ]; then echo "TIMEOUT idx $i" >> out/tables/timeouts.txt; fi
}
export -f one
printf '%s\n' "$@" | xargs -P "$PAR" -I{} bash -c 'one {}'
echo TABLES DONE
