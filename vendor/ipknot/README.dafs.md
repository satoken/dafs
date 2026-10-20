# Bundled reusable IPknot decoder

Upstream: https://github.com/satoken/ipknot, dev commit
`fb0071d6629b06027d936130c8d724a740826e40`. License: GPL-3.0-or-later.

The upstream files are copied without edits. `UPSTREAM.json` records each
file's SHA-256. This snapshot contains `decoder/`, the shared DP/solver source
files referenced by its CMake target, CMake solver finders, and `COPYING`.
It excludes the IPknot CLI and probability engines.

DAFS defaults to this source tree. `IPKNOT_SOURCE_DIR` can still override the
snapshot, or be empty to select an installed `IPknotDecoder` package.

To update, copy these same paths from the chosen upstream revision and refresh
`UPSTREAM.json`; run DAFS's decoder, certificate, isolated-pair, and pipeline
tests before accepting the new revision.
