# Versioning

ICE follows [Semantic Versioning 2.0.0](https://semver.org/), as ORION, FLINT and ATLAS do:
a release is a tag `vX.Y.Z` on `main`, and its changes are listed in
[CHANGELOG.md](https://github.com/open-hydra/ICE/blob/main/CHANGELOG.md). Before 1.0.0 a breaking
change raises `Y`, anything else `Z`.

## How ICE requires its libraries

`CMakeLists.txt` names the releases ICE is tested with:

| Library | Requirement | How |
|---|---|---|
| ORION | `ICE_ORION_VERSION 2.1`: 2.1.0 and later, up to 3.0.0 excluded | `semver_dependency` of `cmake/SemVerDeps.cmake` |
| FiNeR | `ICE_FINER_TAG v2.3.0` | the submodule `lib/FiNeR`, checked against its `VERSION` file |

`cmake/SemVerDeps.cmake` is FLINT's file, copied unmodified: update it from FLINT, never here.
A requirement accepts the releases from its own up to the next change of its leftmost nonzero
number (`2.1` accepts 2.1.0 to 3.0.0 excluded, `0.3` accepts 0.3.x only). Only `vX.Y.Z` tags are
releases; a tag such as `v2.2.0-rc1` is never selected.

ICE requires the oldest release that has everything it uses. ORION `2.1`: ORION 2.0.0 removed
`tec_read_points_multivars`, and the property table is now read with
`tec_read_structured_multiblock`, whose Tecplot reader 2.1.0 fixed.

FiNeR's own dependencies are fetched as FiNeR's author does (`fobis fetch`: their default
branches, `ICE_FINER_DEPS=FREE`). If an upstream change breaks ICE, `ICE_FINER_DEPS=BACKUP`
checks out the releases of `ICE_FINER_BACKUP_PINS`, the ones ATLAS and MOSE share.

## Which release a build gets

The release of ORION is, in order:

1. the clone given with `-DFETCHCONTENT_SOURCE_DIR_ORION=<folder>` (`install.sh --include-orion`),
   as it is; the configure stops if it is not a release that ICE accepts;
2. the release in `deps.lock`, if it satisfies the requirement;
3. the newest release on GitHub that satisfies it.

The configure checks the release out in `lib/ORION` and writes it to `deps.lock`; both are git
ignored. `-DSEMVER_UPDATE=ALL` (`install.sh build --update-deps`) selects the newest one again.
`deps.lock` is not committed: in hydra, the top of the build tree, hydra chooses the releases.

## The requirements in hydra

Hydra builds ICE with its own ORION and FiNeR. When hydra selects ORION with SemVerDeps, ICE's
requirement is checked against that release, so:

- a new **major** ORION in ICE's requirement is a major (before 1.0.0, minor) release of ICE;
- a later **minor** ORION (from `2.1` to `2.2`) is a minor (before 1.0.0, patch) release;
- another FiNeR release is a major (before 1.0.0, minor) release of ICE.

## Making a release

1. Choose the number from the `[Unreleased]` entries of `CHANGELOG.md`.
2. Rename `[Unreleased]` to `[X.Y.Z] - YYYY-MM-DD` and open a new, empty `[Unreleased]` above it.
3. Merge into `main` through a pull request.
4. Tag the merge commit and push the tag:
   ```bash
   git fetch origin
   git tag vX.Y.Z origin/main
   git push origin vX.Y.Z
   ```

!!! warning "A published tag stays"
    A pushed tag is never moved or deleted: other projects' lock files name it.
