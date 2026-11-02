# Amyloid Modifier for ChimeraX

Amyloid Modifier provides a convenient graphical workflow for preparing multilayer amyloid fibril models in ChimeraX: identifying layers and protofilaments, trimming or expanding fibrils, renaming chains, reassessing secondary structure, and generating ISOLDE torsion-restraint commands.

**Author:** Songxian Ou, Ruben Hervas Millan's Lab, The University of Hong Kong.

[ToolShed page](https://cxtoolshed.rbvi.ucsf.edu/apps/chimeraxamyloidmodifier) · [Latest GitHub release](https://github.com/songxou0909/AmyloidModifier/releases/latest) · [Video tutorial](https://www.youtube.com/watch?v=X26wr-cCPcM)

## Install

Install Amyloid Modifier through the [ChimeraX ToolShed](https://cxtoolshed.rbvi.ucsf.edu/apps/chimeraxamyloidmodifier), or download a wheel below and use ChimeraX's `toolshed install` command with its local file path. Launch it from **Tools → Structure Editing → Amyloid Modifier**.

The wheels include the Python source and package metadata used by ChimeraX to identify the bundle and its dependencies. The original build-time `bundle_info.xml` is not needed to install a wheel.

## Citation

If you used Amyloid Modifier and would like to cite it, please use the reference for the version you used. **Citations are appreciated, but not expected.**

Use GitHub's **Cite this repository** button or the citation on the corresponding Zenodo record. Please also follow the citation guidance for ChimeraX and any underlying tools you use.

## Upload a new version (for the maintainer)

1. Build the new wheel with a **new, higher version number** in the bundle's build configuration. For example, after `2.2.2`, use `2.2.3` or `2.3.0`. Renaming the wheel alone does not change its internal version.
2. Open the [`wheels` folder on GitHub](https://github.com/songxou0909/AmyloidModifier/tree/main/wheels).
3. Choose **Add file → Upload files**, drag in the new `.whl`, and choose **Commit changes** directly to `main`. A short message such as “Add 2.2.3” is enough; detailed release notes are optional.
4. Open the [Actions tab](https://github.com/songxou0909/AmyloidModifier/actions) and wait for **Publish wheel release** to finish. It updates the version list, checksums, and citation metadata, then publishes a GitHub Release. Zenodo automatically archives that release and assigns its DOI.

Upload one new version at a time. Keep the older wheels: do not replace or delete a published file. If you find a mistake, build another version. A repeated upload of identical bytes is harmless; changing an already catalogued wheel makes validation stop.

No separate GitHub Release form or Zenodo upload is needed. Uploading here does **not** update the ChimeraX ToolShed; submit there separately if you also want ChimeraX users to receive the new version through ToolShed.

If an upload fails, check its Actions run. After correcting the issue, you can choose **Actions → Publish wheel release → Run workflow** to retry. Existing published releases are left unchanged.

## Version archive

The initial collection contains all 26 versions available from ToolShed on 2026-11-02. These original wheels are preserved byte for byte, with their download URLs and release notes in [`releases.json`](releases.json). [`SHA256SUMS.txt`](SHA256SUMS.txt) records their checksums.

The initial GitHub release is a snapshot containing the whole collection; the older wheels do not each have a separate GitHub Release or Zenodo DOI. Future releases archive the growing collection and identify the latest version.

<!-- WHEELS:START -->
Latest: **2.2.2**. 26 archived wheels.

| Version | Released | Download |
| --- | --- | --- |
| 2.2.2 | 2026-10-27 | [Wheel](wheels/chimerax_amyloidmodifier-2.2.2-py3-none-any.whl) |
| 2.2.1 | 2026-10-16 | [Wheel](wheels/ChimeraX_AmyloidModifier-2.2.1-py3-none-any.whl) |
| 2.2.0 | 2026-10-14 | [Wheel](wheels/ChimeraX_AmyloidModifier-2.2.0-py3-none-any.whl) |
| 2.1.3 | 2026-10-13 | [Wheel](wheels/chimerax_amyloidmodifier-2.1.3-py3-none-any.whl) |
| 2.1.2 | 2026-10-13 | [Wheel](wheels/chimerax_amyloidmodifier-2.1.2-py3-none-any.whl) |
| 2.1.1 | 2026-09-25 | [Wheel](wheels/chimerax_amyloidmodifier-2.1.1-py3-none-any.whl) |
| 2.1 | 2026-09-24 | [Wheel](wheels/chimerax_amyloidmodifier-2.1-py3-none-any.whl) |
| 2.0.1 | 2026-09-24 | [Wheel](wheels/ChimeraX_AmyloidModifier-2.0.1-py3-none-any.whl) |
| 2.0 | 2026-09-24 | [Wheel](wheels/ChimeraX_AmyloidModifier-2.0-py3-none-any.whl) |
| 1.16 | 2026-09-14 | [Wheel](wheels/chimerax_amyloidmodifier-1.16-py3-none-any.whl) |
| 1.15 | 2026-07-22 | [Wheel](wheels/chimerax_amyloidmodifier-1.15-py3-none-any.whl) |
| 1.14 | 2026-05-22 | [Wheel](wheels/chimerax_amyloidmodifier-1.14-py3-none-any.whl) |
| 1.13 | 2026-04-15 | [Wheel](wheels/ChimeraX_AmyloidModifier-1.13-py3-none-any.whl) |
| 1.12 | 2026-03-31 | [Wheel](wheels/ChimeraX_AmyloidModifier-1.12-py3-none-any.whl) |
| 1.11 | 2026-03-17 | [Wheel](wheels/ChimeraX_AmyloidModifier-1.11-py3-none-any.whl) |
| 1.10 | 2026-03-16 | [Wheel](wheels/ChimeraX_AmyloidModifier-1.10-py3-none-any.whl) |
| 1.9 | 2026-03-15 | [Wheel](wheels/chimerax_amyloidmodifier-1.9-py3-none-any.whl) |
| 1.8 | 2026-03-15 | [Wheel](wheels/chimerax_amyloidmodifier-1.8-py3-none-any.whl) |
| 1.7 | 2026-03-15 | [Wheel](wheels/ChimeraX_AmyloidModifier-1.7-py3-none-any.whl) |
| 1.6 | 2026-03-15 | [Wheel](wheels/ChimeraX_AmyloidModifier-1.6-py3-none-any.whl) |
| 1.5 | 2026-03-13 | [Wheel](wheels/ChimeraX_AmyloidModifier-1.5-py3-none-any.whl) |
| 1.4 | 2026-03-12 | [Wheel](wheels/chimerax_amyloidmodifier-1.4-py3-none-any.whl) |
| 1.3 | 2026-03-11 | [Wheel](wheels/ChimeraX_AmyloidModifier-1.3-py3-none-any.whl) |
| 1.2 | 2026-03-11 | [Wheel](wheels/chimerax_amyloidmodifier-1.2-py3-none-any.whl) |
| 1.1 | 2026-03-10 | [Wheel](wheels/chimerax_amyloidmodifier-1.1-py3-none-any.whl) |
| 1.0 | 2026-03-10 | [Wheel](wheels/chimerax_amyloidmodifier-1.0-py3-none-any.whl) |

<!-- WHEELS:END -->

## Using the tool

1. Open a model in ChimeraX and select **Load and Copy from ChimeraX**.
2. Adjust the Z-shift and XY-shift limits when needed; enable **Anti** for antiparallel pairs.
3. Select the desired central layers and trim or expand the model.
4. Rename chains, re-evaluate beta-sheet structure, and optionally export ISOLDE restraint commands.

Trimming modifies chain-related records; unrelated records such as AUTHOR and REMARKS are not rewritten. Expansion uses the model's calculated twist and rise, so inspect the result against your density map.

## Source and license

[MIT License](LICENSE). Copyright (c) 2026 Songxian.

The Python source is inside each wheel (a ZIP-format archive). Extract a wheel with a ZIP utility to inspect the `.py` files. Historical wheels retain their original embedded metadata and notices; the accompanying repository license uses the author's requested name, Songxian. Dependencies retain their own licenses.

This repository stores release distributions. Keep your development source and build configuration separately when building future wheels. The automated workflow and catalogue helper are in `.github/`.
