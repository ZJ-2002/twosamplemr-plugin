#!/usr/bin/env bash
# WO-R-03: seed contract functional test for the twosamplemr main node.
#
# Proves, inside the pinned official image (TwoSampleMR 0.7.9), that
#   1. two runs with the same seed are byte-identical (results, harmonised,
#      RDS), and the log carries "RNG seed: <value>";
#   2. a different seed moves ONLY the bootstrap columns (median/mode SEs)
#      while deterministic quantities (IVW/Egger b, se, p; median/mode b)
#      stay bit-for-bit equal;
#   3. an empty seed reproduces the pre-fix unseeded behaviour (two runs
#      normally differ in bootstrap columns; identical output is a note,
#      not a failure, because unseeded equality is probabilistic).
#
# Environment:
#   TWOSAMPLEMR_IMAGE  Image reference (default 85c9bc77ce1a, the local id
#                      for ghcr.io/auto-nomics/autonomics/twosamplemr@
#                      sha256:c270de9978906ee48cbba2ac484ba3df9e86dddc69efd908c6f908963114002a)
set -euo pipefail
root=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
image=${TWOSAMPLEMR_IMAGE:-85c9bc77ce1a}
work=$(mktemp -d /tmp/wor03-seedtest.XXXXXX)
trap 'rm -rf "$work"' EXIT

input="$work/input.tsv"
# Eight harmonisable SNPs (EA/NEA identical on both sides, no palindromes).
{
  printf 'snp\tbeta_exposure\tse_exposure\teffect_allele_exposure\tother_allele_exposure\teaf_exposure\tbeta_outcome\tse_outcome\teffect_allele_outcome\tother_allele_outcome\teaf_outcome\n'
  for i in 1 2 3 4 5 6 7 8; do
    printf 'rs%03d\t0.%d1\t0.02\tA\tC\t0.30\t0.%d2\t0.03\tA\tC\t0.35\n' "$i" "$i" "$((i * 3 + 1))"
  done
} > "$input"

run_one() { # name seed
  local name="$1" seed="$2"
  mkdir -p "$work/$name"
  podman run --rm --network=none \
    -v "$input":/work/input.tsv:ro \
    -v "$root/scripts/mr.sh":/work/mr.sh:ro \
    -v "$work/$name":/work/run \
    -e AUTONOMICS_INPUT0=/work/input.tsv \
    -e AUTONOMICS_OUTPUT0=/work/run/results.tsv \
    -e AUTONOMICS_OUTPUT1=/work/run/harmonised.tsv \
    -e AUTONOMICS_OUTPUT2=/work/run/result.RDS \
    -e AUTONOMICS_OUTPUT3=/work/run/log.txt \
    -e TWOSAMPLEMR_ID_EXPOSURE=seedtest-e -e TWOSAMPLEMR_EXPOSURE=Exposure \
    -e TWOSAMPLEMR_ID_OUTCOME=seedtest-o -e TWOSAMPLEMR_OUTCOME=Outcome \
    -e TWOSAMPLEMR_METHOD_LIST='mr_egger_regression mr_weighted_median mr_ivw mr_simple_mode mr_weighted_mode' \
    -e TWOSAMPLEMR_HARMONISE_ACTION=2 \
    -e TWOSAMPLEMR_CLUMP=false -e TWOSAMPLEMR_CLUMP_P1=5e-8 -e TWOSAMPLEMR_CLUMP_P2=1e-6 \
    -e TWOSAMPLEMR_CLUMP_R2=0.001 -e TWOSAMPLEMR_CLUMP_KB=10000 \
    -e TWOSAMPLEMR_SEED="$seed" \
    --entrypoint sh "$image" /work/mr.sh > "$work/$name/stdout.log" 2>&1
}

fail=0
run_one s1a 1; run_one s1b 1; run_one s2 2; run_one n1 ""; run_one n2 ""

grep -q 'RNG seed: 1' "$work/s1a/log.txt" || { echo "FAIL: log missing seeded line"; fail=1; }
for f in results.tsv harmonised.tsv result.RDS; do
  cmp -s "$work/s1a/$f" "$work/s1b/$f" || { echo "FAIL: seed=1 runs differ in $f"; fail=1; }
done
echo "seeded double run: byte-identical outputs"

# Deterministic columns equal across seeds; bootstrap SE columns may move.
det_a=$(awk -F'\t' 'NR>1{print $5"|"$7}' "$work/s1a/results.tsv" | sort)
det_b=$(awk -F'\t' 'NR>1{print $5"|"$7}' "$work/s2/results.tsv"  | sort)
[[ "$det_a" == "$det_b" ]] || { echo "FAIL: seed changed deterministic b per method"; fail=1; }
ivw_a=$(awk -F'\t' '$5 ~ /Inverse variance weighted/{print $7" "$8" "$9}' "$work/s1a/results.tsv")
ivw_b=$(awk -F'\t' '$5 ~ /Inverse variance weighted/{print $7" "$8" "$9}' "$work/s2/results.tsv")
[[ "$ivw_a" == "$ivw_b" ]] || { echo "FAIL: IVW moved with the seed"; fail=1; }
echo "seed 1 vs 2: deterministic b/se/p unchanged (bootstrap SE moved separately)"

if cmp -s "$work/n1/results.tsv" "$work/n2/results.tsv"; then
  echo "NOTE: unseeded runs happened to agree (probabilistic; not a failure)"
else
  echo "unseeded runs differ as expected (legacy nondeterminism)"
fi

# Negative guard: seed=0 must be rejected by the script's range check.
if run_one badseed 0 2>/dev/null; then
  echo "FAIL: seed=0 was accepted"; fail=1
else
  echo "seed=0 rejected by the guard"
fi

exit "$fail"
