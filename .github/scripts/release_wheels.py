"""Catalogue distribution wheels without importing or executing their code."""

import argparse
import datetime
import email
import hashlib
import html
import json
from pathlib import Path
import re
import subprocess
import tempfile
import zipfile

ROOT = Path(__file__).resolve().parents[2]
TITLE = "Amyloid Modifier for ChimeraX"
REPO = "https://github.com/songxou0909/AmyloidModifier"
TOOLSHED = "https://cxtoolshed.rbvi.ucsf.edu/apps/chimeraxamyloidmodifier"
MESSAGE = ("If you used Amyloid Modifier and would like to cite it, please use the "
           "reference below for the version you used. Citations are appreciated, "
           "but not expected.")
ABSTRACT = ("Amyloid Modifier is a ChimeraX add-on that streamlines routine preparation "
            "and editing of multilayer amyloid fibril models. It provides a graphical "
            "interface for identifying layers and protofilaments, trimming and expanding "
            "fibril models, renaming chains, invoking secondary-structure reassignment, "
            "and generating ISOLDE torsion-restraint command files.")
KEYWORDS = ["ChimeraX", "amyloid", "structural biology", "cryo-EM", "molecular modeling"]
RESPONSIBILITY = ("Amyloid Modifier is solely the work of Songxian Ou. "
                  "It does not represent or imply endorsement by Ruben Hervas Millan's lab. "
                  "The lab is not responsible for the software, its use, or any resulting outcomes. "
                  "The software is provided under the terms of the MIT License.")


def version_key(version):
    if not re.fullmatch(r"\d+(?:\.\d+){1,3}", version):
        raise ValueError(f"Expected a numeric release version, got {version!r}")
    parts = tuple(map(int, version.split(".")))
    return parts + (0,) * (4 - len(parts))


def catalogue(root=ROOT):
    manifest = root / "releases.json"
    previous = json.loads(manifest.read_text(encoding="utf-8")) if manifest.exists() else []
    by_file = {r["filename"]: r for r in previous}
    wheel_paths = sorted((root / "wheels").glob("*.whl"))
    missing = set(by_file) - {p.name for p in wheel_paths}
    if missing:
        raise ValueError(f"Keep existing wheels; missing: {sorted(missing)}")
    result, versions = [], set()
    for path in wheel_paths:
        data = path.read_bytes()
        digest = hashlib.sha256(data).hexdigest()
        old = by_file.get(path.name, {})
        if old and old["sha256"] != digest:
            raise ValueError(f"Do not replace a released wheel: {path.name}. Build a new version.")
        with zipfile.ZipFile(path) as wheel:
            names = wheel.namelist()
            metadata_files = [n for n in names if n.endswith(".dist-info/METADATA")]
            if len(metadata_files) != 1:
                raise ValueError(f"Expected exactly one METADATA file in {path.name}")
            if wheel.getinfo(metadata_files[0]).file_size > 1024 * 1024:
                raise ValueError(f"Unexpectedly large METADATA in {path.name}")
            metadata = email.message_from_bytes(wheel.read(metadata_files[0]))
            if metadata.get("Name") != "ChimeraX-AmyloidModifier":
                raise ValueError(f"Wrong package in {path.name}")
            version = metadata["Version"]
            key = version_key(version)
            if key in versions:
                raise ValueError(f"Duplicate version: {version}")
            versions.add(key)
            if path.name.split("-")[0].lower() != "chimerax_amyloidmodifier" or version_key(path.name.split("-")[1]) != key:
                raise ValueError(f"Filename and wheel metadata disagree: {path.name}")
            if old and version_key(old["version"]) != key:
                raise ValueError(f"Version changed in {path.name}")
            if not old and previous and key <= max(version_key(r["version"]) for r in previous):
                raise ValueError(f"New uploads must have a newer version: {version}")
            python_files = [n for n in names if n.endswith(".py")]
            if not python_files:
                raise ValueError(f"No readable Python source in {path.name}")
            record = dict(old)
            record.update(version=version, filename=path.name, bytes=len(data), sha256=digest,
                          md5=hashlib.md5(data).hexdigest(), python_files=python_files,
                          compiled_files=[n for n in names if n.endswith((".pyd", ".so", ".dll", ".pyc"))],
                          license_files=[n for n in names if "license" in n.lower()],
                          requires_dist=metadata.get_all("Requires-Dist", []))
            record.setdefault("released", datetime.datetime.now(datetime.timezone.utc).date().isoformat())
            record.setdefault("notes", "")
            record.setdefault("download_url", f"{REPO}/releases/download/v{version}/{path.name}")
            result.append(record)
    if not result:
        raise ValueError("No wheels found")
    return sorted(result, key=lambda r: version_key(r["version"]))


def write_json(path, value):
    path.write_text(json.dumps(value, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")


def prepare(root=ROOT):
    records = catalogue(root)
    latest = records[-1]
    write_json(root / "releases.json", records)
    (root / "SHA256SUMS.txt").write_text("".join(f'{r["sha256"]}  wheels/{r["filename"]}\n' for r in records), encoding="utf-8")
    cff_path = root / "CITATION.cff"
    # Keep an assigned DOI only for the exact version to which it belongs.
    old_cff = cff_path.read_text(encoding="utf-8") if cff_path.exists() else ""
    old_version = re.search(r'^version: "([^"]+)"$', old_cff, re.M)
    old_doi = re.search(r'^doi: "([^"]+)"$', old_cff, re.M)
    doi_line = f'doi: "{old_doi[1]}"\n' if old_doi and old_version and old_version[1] == latest["version"] else ""
    cff_path.write_text(
        f'cff-version: 1.2.0\nmessage: {json.dumps(MESSAGE)}\ntype: software\n'
        f'title: "{TITLE}"\nauthors:\n  - family-names: Ou\n    given-names: Songxian\n'
        f'version: "{latest["version"]}"\ndate-released: {latest["released"][:10]}\n{doi_line}'
        f'url: "{TOOLSHED}"\nrepository-code: "{REPO}"\nlicense: MIT\n'
        f'abstract: {json.dumps(ABSTRACT)}\nkeywords:\n' + "".join(f"  - {k}\n" for k in KEYWORDS), encoding="utf-8")
    write_json(root / ".zenodo.json", {
        "title": TITLE, "upload_type": "software", "version": latest["version"],
        "publication_date": latest["released"][:10], "license": "MIT",
        "creators": [{"name": "Ou, Songxian"}],
        "description": f"<p>{html.escape(ABSTRACT)}</p><p>{html.escape(MESSAGE)}</p>"
                       f"<p>{html.escape(RESPONSIBILITY)}</p>"
                       "<p>Please also follow the citation guidance for ChimeraX and any underlying tools you use.</p>"
                       f"<p>This repository snapshot contains {len(records)} versioned wheels, including the latest version "
                       f"{latest['version']}. Each wheel contains Python source and installation metadata. "
                       "Historical ToolShed wheels are preserved byte for byte. Copyright (c) 2026 Songxian.</p>",
        "keywords": KEYWORDS,
        "related_identifiers": [{"identifier": TOOLSHED, "relation": "isSupplementTo", "scheme": "url"}],
    })
    table = f'Latest: **{latest["version"]}**. {len(records)} archived wheels.\n\n| Version | Released | Download |\n| --- | --- | --- |\n'
    table += "".join(f'| {r["version"]} | {r["released"][:10]} | [Wheel](wheels/{r["filename"]}) |\n' for r in reversed(records))
    readme_path = root / "README.md"
    readme = readme_path.read_text(encoding="utf-8-sig")
    pattern = r"<!-- WHEELS:START -->.*?<!-- WHEELS:END -->"
    if len(re.findall(pattern, readme, re.S)) != 1:
        raise ValueError("README must have exactly one WHEELS section")
    readme_path.write_text(re.sub(pattern, f"<!-- WHEELS:START -->\n{table}\n<!-- WHEELS:END -->", readme, flags=re.S), encoding="utf-8")
    print(f"Validated {len(records)} wheels; latest version {latest['version']}")


def publish(root=ROOT):
    latest = catalogue(root)[-1]
    tag = "v" + latest["version"]
    def run(*args, **kwargs):
        return subprocess.run(args, cwd=root, text=True, check=True, **kwargs)
    existing = subprocess.run(["gh", "release", "view", tag, "--json", "isDraft"], cwd=root, text=True, capture_output=True)
    if existing.returncode == 0 and not json.loads(existing.stdout)["isDraft"]:
        print(f"Release {tag} is already published; leaving it unchanged.")
        return
    notes = (f"Amyloid Modifier for ChimeraX {latest['version']}.\n\n"
             "Download the wheel below for installation. The source-code ZIP/TAR snapshot also contains "
             "all historical wheels in `wheels/`, their Python source, the MIT license, and checksums.\n\n"
             f"{MESSAGE}\n\nZenodo automatically archives this repository snapshot after publication.\n")
    assets = [str(root / "wheels" / latest["filename"]), str(root / "LICENSE"), str(root / "SHA256SUMS.txt")]
    if existing.returncode != 0:
        commit = run("git", "rev-parse", "HEAD", capture_output=True).stdout.strip()
        with tempfile.TemporaryDirectory() as temp:
            notes_file = Path(temp) / "release-notes.md"
            notes_file.write_text(notes, encoding="utf-8")
            run("gh", "release", "create", tag, "--draft", "--target", commit,
                "--title", f"Amyloid Modifier {latest['version']}", "--notes-file", str(notes_file))
    run("gh", "release", "upload", tag, *assets, "--clobber")
    run("gh", "release", "edit", tag, "--draft=false", "--latest")
    print(f"Published {tag}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("command", choices=["check", "prepare", "publish"])
    args = parser.parse_args()
    if args.command == "check":
        records = catalogue()
        print(f"Verified {len(records)} wheels; latest {records[-1]['version']}")
    elif args.command == "prepare":
        prepare()
    else:
        publish()
