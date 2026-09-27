# twosamplemr plugin

Migrated from the legacy `twosamplemr_container` +
`twosamplemr_harmonise_container` wrapper pair in nodes-io. One directory =
one plugin family = one git-able unit; the family holds TWO node kinds
because they share one image (official R TwoSampleMR 0.7.9 plus the
official PLINK2 binary baked into the same image) and one panel binding.

## Layout

- `manifest.toml` — node kinds `twosamplemr` (main MR analysis) and
  `twosamplemr_harmonise` (two-file harmonisation): params, ports, the
  1000G EUR PLINK binary panel binding, image provenance
- `scripts/mr.sh` — the main node's hybrid shell + R execution script
  (PLINK2 clumping stage feeding the TwoSampleMR R stage), referenced
  relatively and inlined by the loader at startup
- `scripts/harmonise.sh` — the harmonise node's R script (pure R source
  executed by the `Rscript` interpreter, the coloc `abf.sh` pattern;
  the `.sh` suffix matches the wave's script-file convention)
- `Dockerfile` — image build provenance (moved verbatim from
  `containers/twosamplemr/`; build + push still via GHCR)
- `test_twosamplemr.sh` — image baseline: builds the image and
  reproduces the upstream `mr()` golden estimates for IVW (`0.4459`) and
  MR Egger (`0.5025`) from the official `test_commondata.RData`;
  `root=` repointed to this directory
- `container_README.md` — the original containers/twosamplemr/README.md
  (verbatim, modulo the regenerated `$AUTONOMICS_IMAGE_PREFIX` and
  wrapper-path spellings)

**Fixtures did not move.** The chr22 instrument fixture remains in the
autonomics repository at `containers/twosamplemr/fixtures/` because the
still-live nodes-io Rust integration test
(`crates/node-bundles/nodes-io/tests/container_file_flow.rs`, the
`#[ignore]`d `real_catalog_backed_official_twosamplemr_runs_with_container_backend`)
reads it from the repository-relative path and honors the
`AUTONOMICS_TWOSAMPLEMR_IT_SUMSTATS` override.

## Provenance

- Image:
  `ghcr.io/auto-nomics/autonomics/twosamplemr@sha256:c270de9978906ee48cbba2ac484ba3df9e86dddc69efd908c6f908963114002a`,
  tag `0.7.9`, from the moved `Dockerfile` (base `rocker/r-ver:4.5.1`).
- Upstream: official R
  [TwoSampleMR](https://github.com/MRCIEU/TwoSampleMR) 0.7.9 (source
  archive SHA-256 `6848c344c5eead601ff52e9a88b2c23c2b61b0ce4f848e331e7494b657be64e5`,
  revision `3d119f2`) on a 2026-09-13 Posit Package Manager CRAN
  snapshot, with MRMix `56afdb2`, RadialMR `a30ff11`, and MR-PRESSO
  `3e3c92d`; plus the official PLINK2 v2.0.0-a.6.26 binary
  (`plink2_linux_avx2.zip`, asset SHA-256
  `f578a450af382d7dd6665aecf0ca1d280971c2b3d5bb5556efbf9266c4c8da0f`).
- License: MIT per the image label (the vendored PLINK2 binary is
  GPL-3.0; see `container_README.md`).

## Two-node layout

- `twosamplemr` — one SNP-merged exposure/outcome sumstats File in;
  optional official PLINK2 clumping against `/panels/plink_ref`, then
  `TwoSampleMR::harmonise_data()` + `TwoSampleMR::mr()`. Four outputs:
  results TSV, harmonised TSV, full RDS result, execution log.
- `twosamplemr_harmonise` — separate exposure + outcome sumstats Files
  in; `read_exposure_data` + `read_outcome_data` + `harmonise_data`,
  with the native dotted column names rewritten to snake_case so the
  output feeds the `twosamplemr` node directly. Two outputs: merged TSV
  and run log.

## Migration parity

The golden test (`crates/container-plugin/tests/twosamplemr_migration.rs`)
compares each compiled `ContainerCommandSpec` against the legacy Rust
wrapper: image, outputs, panel bundle ids/mounts, resources, and timeout
are byte-exact per kind. Deliberate deltas, recorded for honesty:

- **Kind renames**: `twosamplemr_container` → `twosamplemr` and
  `twosamplemr_harmonise_container` → `twosamplemr_harmonise` (the
  `_container` suffix only ever distinguished a wrapper from a native
  Rust port, which is gone). The artifact prefixes follow the kinds
  (`/artifacts/twosamplemr_container` → `/artifacts/twosamplemr`), the
  wave-wide `/artifacts/{kind}` convention. DAG specs referencing the
  old kinds must be regenerated.
- **Command vectors**: the main node's legacy `["sh", "-c"]` becomes
  interpreter `sh` (the executor inserts the materialized script at
  index 1; the trailing `-c` was a no-op positional under
  `sh <script>` — same delta as the plink2 migration). The harmonise
  node's `["Rscript"]` is reproduced byte-for-byte.
- **`timeout_secs` / `artifact_prefix` are node-level constants**
  (1800 s / 600 s, `/artifacts/twosamplemr*`) instead of per-instance
  spec params; the legacy specs accepted per-node overrides, the plugin
  DSL does not. Submitting either name fails the
  `additionalProperties: false` gate.
- **Family-scoped panels**: the legacy harmonise wrapper attached no
  panels; the manifest's `[[panels]]` is family-level, so the harmonise
  node now carries (and its DAG must bind) the `plink_ref` bundle even
  though its script never reads `/panels` — the v0 DSL has no per-node
  panel opt-out (the `ldsc_munge` propagation precedent, inverted: here
  the extra binding is new, not parity-preserving). Mount is read-only
  and unused.
- **Validation moved into the scripts** (failure point moves from
  registry build to container start, the coloc/deseq2 precedent): the
  label-emptiness checks (`id_exposure cannot be empty`, ...), the
  `method_list` SUPPORTED_METHODS membership (no enum param type; the
  manifest carries `minItems = 1` for the empty case), and the
  cross-field `clump_p2 >= clump_p1` rule (checked with `awk` before
  the PLINK2 stage). The DSL carries what it can: `harmonise_action`
  1..=3, `clump_*` thresholds in (0, 1], `clump_kb` > 0, and
  `chr` 1..=22 via numeric bounds.
- **R string interpolation becomes env**: the legacy wrapper baked
  labels into the R source with Rust `{:?}` quoting; the plugin passes
  them through `TWOSAMPLEMR_*` environment variables read with
  `Sys.getenv()` — identical runtime values, no escaping layer.
- **Threshold spelling**: the legacy clump stage rendered
  `clump_p1`/`clump_p2` with `{:.0e}` (`5e-08`); the env renderer
  spells the same f64 values `5e-8` / `1e-6` (serde_json/ryu). Equal
  after float parsing; PLINK2 accepts both.
- **method_list travels space-joined** (`TWOSAMPLEMR_METHOD_LIST`);
  the R stage splits on a single space. Official method names never
  contain spaces, so the join is injective over legacy-legal input.
- **Ports**: output ports gain file-stem labels because the manifest
  pipeline always names output ports; the legacy ports were unlabeled —
  the same accepted delta as the mvmr/deseq2 migrations.

The scripts preserve the legacy semantics exactly: the clumping stage
keeps the official tokens (`--bfile /panels/plink_ref/1000G.EUR.QC.<chr>`,
`--clump /work/clump_input.tsv`, `--clump-id-field SNP
--clump-p-field P`, `--clump-p1/p2/r2/kb`, `--threads 1`, the
`seq 1 22` / `seq <chr> <chr>` loop with the optional `--chr`, per-
chromosome `.clumps` concatenation into `/work/all.clumps`, and the
plink2 stdout/stderr append into the log output), the R stage keeps the
`TwoSampleMR::harmonise_data()` / `TwoSampleMR::mr()` calls, the
pval-or-normal-approximation p-value derivation, the pre-clumped
`clump = false` path, and the four-artifact epilogue. There is no gzip
stanza in either script: neither legacy wrapper ever decompressed
inputs. The harmonise script keeps the 0-row harmonisation guard and
the `SNP` → `snp` / dotted-to-snake_case column rewrite verbatim.

## Install

```sh
export NODE_PLUGINS_ROOT=/mnt/projects/node-plugins
cargo test -p container-plugin --test twosamplemr_migration   # golden parity
```

Declare the family in `~/.autonomics/plugins.toml` (or the deployment
equivalent) once it is pushed to git:

```toml
[[plugin]]
name = "twosamplemr"
git = "git@github.com:auto-nomics/twosamplemr-plugin.git"
rev = "<pinned commit SHA>"
```
