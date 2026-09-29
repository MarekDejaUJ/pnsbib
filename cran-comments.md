# pnsbib 0.2.0: submission preparation

## Test environment

Apple Silicon macOS Tahoe 26.6.2, R 4.5.2, rzig 0.2.3 and Zig 0.16.0. The source archive was built and checked in a separate temporary directory with standard compiler flags. The final full `R CMD check --as-cran`, including the PDF manual and vignette rebuild, completed on 29 September 2026.

Result: 0 errors, 0 warnings, 4 notes. All 33,435 R assertions and 32 separate native kernel tests passed. All 32 exports have executable examples and method-reference documentation. The longest Rd example took 0.135 seconds on this machine. The installed workflow vignette does not include installation commands.

## Outstanding notes and checks

1. Incoming URL checks could not resolve external hosts in the checking environment. Online incoming checks must be repeated with network access.
2. The current-time check could not reach its external service.
3. The available HTML Tidy was too old for HTML validation. Repeat with a supported Tidy version.
4. Apple's developer tools produced `xcrun_db` in the temporary directory. Invoking the system `nm` utility independently also creates this file. Confirm the result on a clean macOS runner.

Current R release checks on macOS, Windows and Linux, and R-devel on Linux, are configured in the repository workflow but have not run for this version. The local R version is reported above without claiming it is the current release. This document is a preparation record, not a claim that the outstanding checks have passed.

## Installation and filesystem behavior

Source installation uses R and rzig to build the bundled Zig code. No Python script or manually installed Zig distribution is required. If a compatible compiler is absent, the R build helper obtains the pinned official Zig archive over HTTPS, with a 300-second timeout, size validation and SHA-256 verification. Default compiler storage is session-temporary and is removed after compilation. Persistent caching requires an explicitly supplied `PNSBIB_ZIG_CACHE` path. Build concurrency is limited to two jobs, and the native C wrapper respects R's compiler flags. No compiler acquisition occurs when the installed package is loaded or used.

The automatic setup was verified with the checksum-validated official archive supplied offline, with `ZIG` unset, and with downloader/error/cleanup regression tests. A first-time live HTTPS download remains to be tested. Acceptance of the external compiler requirement remains subject to CRAN review.

Export functions require explicit file or directory arguments. Examples, tests and the vignette write only to temporary paths and clean their generated outputs. No API keys or private data are required for the examples or aggregate reproduction scripts. Private individual-level records are not distributed.

## Scope of the release

Version 0.2.0 adds documentation, method citations, an installed tutorial, executable examples and installation safeguards. Existing R numerical expressions and Zig/C computational sources are unchanged from 0.0.22. Aggregate replays reproduce 570 bound rows, 600 simulation repetitions and 1,764 independently checked LP certificates, with seven result tables byte-identical to the previous release. License: GPL-3, with retained notices for the bundled rzig framework.

No submission to CRAN has been made.
