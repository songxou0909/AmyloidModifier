import os
import math
from io import StringIO
from copy import deepcopy
from html import escape

import json
import tempfile
import hashlib
import base64
from collections import Counter, defaultdict


class StructureEditError(ValueError):
    pass


def _structure_dependencies():
    try:
        import gemmi
        import numpy as np
    except (ImportError, OSError) as exc:
        raise StructureEditError(
            "Amyloid Modifier's structure libraries could not be loaded. "
            "Install the latest Amyloid Modifier update through Tools > "
            "More Tools with dependency installation enabled, then restart "
            "ChimeraX. The installer supplies Gemmi automatically; no commands "
            "are needed. NumPy is included with ChimeraX.\n\n"
            f"Details: {exc}"
        ) from exc
    return gemmi, np


def _cif_rows(block, category):
    data = block.get_mmcif_category(category)
    return [dict(zip(data, values)) for values in zip(*data.values())]


def _set_cif_rows(block, category, rows):
    block.find_mmcif_category(category).erase()
    if rows:
        keys = list(dict.fromkeys(k for row in rows for k in row))
        block.set_mmcif_category(category, {k: [r.get(k) for r in rows] for k in keys})


def _present(value):
    return value is not None and value is not False and value not in ("", ".", "?")


def _null_id(value):
    return str(value) if _present(value) else ""


def _cif_annotation_records(categories):
    if not categories:
        return []
    prefix = 'REMARK 999 AMYLOID_CIF '
    width = 80 - len(prefix)
    encoded = base64.b64encode(json.dumps(categories, ensure_ascii=False).encode('utf-8')).decode('ascii')
    return [prefix + encoded[i:i + width] for i in range(0, len(encoded), width)]


def structure_chain_id(index, max_length=None):
    # Bijective base 62 avoids truncating or reusing chain IDs.
    if not isinstance(index, int) or index < 0:
        raise StructureEditError("Chain index must be a nonnegative integer.")
    alphabet = "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789"
    result = ""
    index += 1
    while index:
        index, digit = divmod(index - 1, len(alphabet))
        result = alphabet[digit] + result
    if max_length is not None and len(result) > max_length:
        raise StructureEditError("PDB chain namespace exhausted; export mmCIF.")
    return result


def _atom_key(row, model=True):
    fields = ("auth_asym_id", "auth_seq_id", "pdbx_PDB_ins_code", "label_comp_id",
              "label_atom_id", "label_alt_id")
    key = tuple(_null_id(row.get(k)) for k in fields)
    return ((_null_id(row.get("pdbx_PDB_model_num")) or "1",) + key) if model else key


def _residue_key(row):
    return tuple(_null_id(row.get(k)) for k in
                 ("auth_asym_id", "auth_seq_id", "pdbx_PDB_ins_code", "label_comp_id"))


def _endpoint_tests(row, prefix):
    # Unspecified CIF endpoint fields are wildcards.
    aliases = {
        "label_asym_id": (prefix + "label_asym_id",), "auth_asym_id": (prefix + "auth_asym_id",),
        "label_seq_id": (prefix + "label_seq_id",), "auth_seq_id": (prefix + "auth_seq_id",),
        "label_comp_id": (prefix + "label_comp_id",), "label_atom_id": (prefix + "label_atom_id",),
        "pdbx_PDB_ins_code": ("pdbx_" + prefix + "PDB_ins_code", prefix + "PDB_ins_code"),
        "label_alt_id": ("pdbx_" + prefix + "label_alt_id",),
    }
    return [(target, str(row[k])) for target, keys in aliases.items() for k in keys if _present(row.get(k))]


class _AtomIndex:

    def __init__(self, atoms):
        self.atoms = atoms
        self.indexes = {}

    def _candidates(self, tests):
        values = dict(tests)
        fields = next((fields for fields in (
            ('auth_asym_id', 'auth_seq_id', 'label_atom_id'),
            ('label_asym_id', 'label_seq_id', 'label_atom_id'),
            ('auth_asym_id', 'auth_seq_id'), ('label_asym_id', 'label_seq_id'),
            ('auth_asym_id',), ('label_asym_id',), ('label_atom_id',),
        ) if all(k in values for k in fields)), None)
        if fields is None:
            return self.atoms
        if fields not in self.indexes:
            index = defaultdict(list)
            for atom in self.atoms:
                index[tuple(str(atom.get(k)) for k in fields)].append(atom)
            self.indexes[fields] = index
        return self.indexes[fields].get(tuple(values[k] for k in fields), ())

    def matches(self, row, prefix):
        tests = _endpoint_tests(row, prefix)
        return [a for a in self._candidates(tests)
                if all(str(a.get(k)) == v for k, v in tests)] if tests else []

    def exists(self, row, prefix):
        tests = _endpoint_tests(row, prefix)
        return bool(tests) and any(all(str(a.get(k)) == v for k, v in tests)
                                   for a in self._candidates(tests))


class StructureBuffer:
    # Coordinates and edit provenance stay in memory until explicitly exported.

    def __init__(self, name, text='', report=None):
        self.name = os.path.basename(name)
        self.text = text
        self.report = deepcopy(report)
        self.format = 'cif' if self.name.lower().endswith(('.cif', '.mmcif')) else 'pdb'

    def copy_from(self, other):
        self.text = other.text
        self.report = deepcopy(other.report)
        self.format = other.format


class StructureEditor:
    # PDB supports two-character IDs by default; strict mode uses one character.
    # Tolerate source annotation problems, keeping their contents in the report.

    # Categories whose meaning is independent of chain instances/coordinates.
    SAFE_CATEGORIES = set("""
        entry audit_conform audit_author citation citation_author citation_editor
        struct struct_keywords exptl exptl_crystal exptl_crystal_grow exptl_crystal_grow_comp
        cell symmetry space_group space_group_symop atom_type chem_comp chem_comp_atom
        chem_comp_bond chem_comp_angle chem_comp_tor chem_comp_chir chem_comp_plane
        chem_comp_plane_atom entity entity_poly entity_poly_seq entity_src_gen
        entity_src_nat pdbx_entity_src_syn pdbx_entity_nonpoly entity_name_com
        pdbx_entity_name_com pdbx_entity_keywords struct_ref struct_conn_type struct_conf_type
        pdbx_nmr_sample_details pdbx_nmr_exptl_sample pdbx_nmr_exptl_sample_conditions
        pdbx_nmr_exptl pdbx_nmr_refine pdbx_nmr_software pdbx_nmr_spectrometer
        em_experiment em_sample_preparation em_vitrification em_imaging em_detector
        em_image_recording em_3d_reconstruction em_helical_entity em_entity_assembly
        em_software software audit em_buffer em_sample_support em_ctf_correction
        em_entity_assembly_naturalsource em_image_processing em_specimen
        pdbx_contact_author pdbx_audit_support pdbx_nmr_details
        pdbx_initial_refinement_model em_image_scans em_buffer_component
        em_entity_assembly_molwt em_entity_assembly_recombinant em_imaging_optics
        em_particle_selection
    """.split())
    CHAIN_CATEGORIES = set("""
        struct_ref_seq struct_ref_seq_dif pdbx_struct_mod_residue
        pdbx_unobs_or_zero_occ_residues pdbx_unobs_or_zero_occ_atoms
        pdbx_poly_seq_scheme pdbx_nonpoly_scheme pdbx_branch_scheme
        struct_conf struct_mon_prot_cis struct_conn struct_site_gen
    """.split())
    SPECIAL_CATEGORIES = {'atom_site', 'atom_site_anisotrop', 'struct_asym', 'struct_sheet',
                          'struct_sheet_range', 'struct_sheet_order', 'pdbx_struct_sheet_hbond',
                          'struct_site', 'audit_syntax', 'amyloid_unrecognized_pdb'}
    INVALIDATED_PREFIXES = (
        "pdbx_validate_", "pdbx_vrpt_", "refine", "pdbx_refine", "struct_ncs",
        "pdbx_struct_assembly", "pdbx_struct_oper", "pdbx_audit_revision",
        "pdbx_database_", "database_", "pdbx_nmr_ensemble", "pdbx_nmr_representative",
        "struct_biol", "atom_sites", "pdbx_struct_special_symmetry", "pdbx_helical_symmetry", "em_3d_fitting",
        "pdbx_entry_details", "em_admin",
    )

    def __init__(self, path, *, metadata_policy="preserve", pdb_mode="extended"):
        gemmi, np = _structure_dependencies()
        if metadata_policy not in ("error", "audit", "preserve") or pdb_mode not in ("extended", "strict"):
            raise StructureEditError("Invalid metadata policy or PDB mode.")
        buffer = path if isinstance(path, StructureBuffer) else None
        in_memory = buffer is not None or hasattr(path, 'read')
        if buffer is not None:
            self.path = 'memory:' + buffer.name
            text = buffer.text
            source_hash = hashlib.sha256(text.encode('utf-8')).hexdigest()
        elif in_memory:
            self.path = str(getattr(path, 'name', '<memory>'))
            text = path.read()
            source_hash = hashlib.sha256(text.encode('utf-8')).hexdigest()
        else:
            self.path = os.path.abspath(os.fspath(path))
            source_hash = self._hash(self.path)
            with open(self.path, encoding='utf-8-sig') as handle:
                text = handle.read()
        self.metadata_policy = metadata_policy
        self.pdb_mode = pdb_mode
        self.report = {"source": self.path, "source_sha256": source_hash,
                       "engine": "AmyloidModifier trial structured I/O", "gemmi": gemmi.__version__,
                       "pdb_mode": pdb_mode, "warnings": [], "metadata": {}, "chain_copies": [],
                       "validation_scope": "syntax, identities, references and coordinate round-trip; not wwPDB deposition validation"}
        previous_path = self.path + '.audit.json'
        previous = deepcopy(buffer.report) if buffer is not None else None
        if previous is not None or (not in_memory and os.path.isfile(previous_path)):
            try:
                if previous is None:
                    with open(previous_path, encoding='utf-8') as handle:
                        previous = json.load(handle)
                if previous.get('output_sha256') == self.report['source_sha256']:
                    fields = ('operation','source','source_sha256','output','output_sha256','chain_copies','geometry','inferred_covalent_links','ensemble_selection')
                    self.report['history'] = previous.get('history', []) + [{k:previous[k] for k in fields if k in previous}]
                    archive = dict(previous.get('inherited_metadata_archive', {}))
                    archive.update({k:v for k,v in previous.get('metadata',{}).items() if 'source_content' in v})
                    if previous.get('atom_coordinate_metadata_archive'):
                        archive['atom_coordinate_metadata:' + previous.get('source_sha256', '')] = {
                            'source_content': previous['atom_coordinate_metadata_archive']}
                    self.report['inherited_metadata_archive'] = archive
                    self.report['source_issues'] = deepcopy(previous.get('source_issues', []))
                    self.report['warnings'].extend(issue['warning'] for issue in self.report['source_issues'])
                    if previous.get('source_pdb_metadata'):
                        self.report['source_pdb_metadata'] = previous['source_pdb_metadata']
                else:
                    self.report['warnings'].append('Previous audit hash differs from input; previous provenance was not inherited.')
            except (OSError, ValueError, TypeError, AttributeError):
                self.report['warnings'].append('Previous audit could not be read; no provenance was inferred from it.')
        first_token_line = next((line.lstrip() for line in text.splitlines()
                                 if line.strip() and not line.lstrip().startswith('#')), '')
        self.is_cif = first_token_line.lower().startswith('data_')
        self.pdb_lines = [] if self.is_cif else text.splitlines()
        self.report['source_record_counts'] = dict(Counter(l[:6].strip() for l in self.pdb_lines))
        known_pdb = set('HEADER OBSLTE TITLE SPLIT CAVEAT COMPND SOURCE KEYWDS EXPDTA NUMMDL MDLTYP AUTHOR REVDAT SPRSDE JRNL REMARK DBREF DBREF1 DBREF2 SEQADV SEQRES MODRES HET HETNAM HETSYN FORMUL HELIX SHEET SSBOND LINK CISPEP SITE CRYST1 ORIGX1 ORIGX2 ORIGX3 SCALE1 SCALE2 SCALE3 MTRIX1 MTRIX2 MTRIX3 MODEL ATOM HETATM ANISOU TER ENDMDL CONECT MASTER END'.split())
        unknown_pdb = set(self.report['source_record_counts']) - known_pdb - {''}
        self.unrecognized_pdb_records = [line for line in self.pdb_lines if line[:6].strip() in unknown_pdb]
        if unknown_pdb and self.metadata_policy == 'error':
            raise StructureEditError('PDB annotations without an edit handler: ' + ', '.join(sorted(unknown_pdb)) + '. Use audit policy to archive them explicitly.')
        if unknown_pdb:
            self._source_warning('unreviewed_pdb_records',
                'PDB annotations without an edit handler: ' + ', '.join(sorted(unknown_pdb)) +
                ('. These entries are preserved unchanged when saving; their references may describe the original structure.'
                 if self.metadata_policy == 'preserve' else
                 '. These annotations are archived in memory and omitted from the working copy; atom coordinates are retained.'),
                records=[{'line': i, 'text': line} for i, line in enumerate(self.pdb_lines, 1)
                         if line[:6].strip() in unknown_pdb])
        if self.is_cif:
            doc = gemmi.cif.read_string(text)
            if len(doc) != 1:
                raise StructureEditError("Select one mmCIF data block before editing; multiple blocks are ambiguous.")
            self.block = doc.sole_block()
            # Keep all categories, quoted values and multiline text in the CIF DOM.
            self.document = doc
            st = gemmi.make_structure_from_block(self.block)
        else:
            st = gemmi.read_pdb_string(text)
            self._check_pdb_atoms(st)
            st.setup_entities()
            st.assign_label_seq_id()
            self.document = st.make_mmcif_document(gemmi.MmcifOutputGroups(True, auth_all=True))
            self.block = self.document.sole_block()
            embedded = ''.join(line[len('REMARK 999 AMYLOID_CIF '):].strip()
                               for line in self.pdb_lines if line.startswith('REMARK 999 AMYLOID_CIF '))
            if embedded:
                try:
                    categories = json.loads(base64.b64decode(embedded, validate=True))
                    for category, rows in categories.items():
                        if not category.startswith('_') or not category.endswith('.') or not isinstance(rows, list):
                            raise ValueError('Invalid category envelope')
                        if self.block.find_mmcif_category(category):
                            raise ValueError('Embedded metadata conflicts with parsed coordinates or metadata')
                        _set_cif_rows(self.block, category, rows)
                except (ValueError, TypeError, AttributeError) as exc:
                    raise StructureEditError('Cannot restore embedded mmCIF annotations: ' + str(exc)) from exc
        if self.metadata_policy == 'preserve':
            if self.unrecognized_pdb_records:
                _set_cif_rows(self.block, '_amyloid_unrecognized_pdb.',
                              [{'id': str(i), 'text': line} for i, line in enumerate(self.unrecognized_pdb_records, 1)])
            else:
                self.unrecognized_pdb_records = [row['text'] for row in _cif_rows(self.block, '_amyloid_unrecognized_pdb.')]
        self.atoms = _cif_rows(self.block, "_atom_site.")
        if not self.atoms:
            raise StructureEditError("No atom_site coordinates found.")
        required = {"id", "label_asym_id", "label_comp_id", "label_atom_id", "Cartn_x", "Cartn_y", "Cartn_z"}
        if not required <= self.atoms[0].keys():
            raise StructureEditError("Missing required atom_site columns: " + ", ".join(sorted(required - self.atoms[0].keys())))
        for row in self.atoms:
            for auth, label in (("auth_asym_id", "label_asym_id"), ("auth_seq_id", "label_seq_id"),
                                ("auth_atom_id", "label_atom_id"), ("auth_comp_id", "label_comp_id")):
                if not _present(row.get(auth)):
                    row[auth] = row.get(label)
            row.setdefault("pdbx_PDB_model_num", "1")
        self.models = list(dict.fromkeys(str(r["pdbx_PDB_model_num"]) for r in self.atoms))
        self.chains = list(dict.fromkeys(_null_id(r["auth_asym_id"]) for r in self.atoms))
        self._validate_atoms(self.atoms)
        self.poly_labels = {r["id"] for r in _cif_rows(self.block, "_struct_asym.")
                            if r.get("entity_id") in {e["entity_id"] for e in _cif_rows(self.block, "_entity_poly.")}}
        if not self.poly_labels:
            self.poly_labels = {r["label_asym_id"] for r in self.atoms if _present(r.get("label_seq_id"))}
        self.pdb_bonds = self._read_conect(st) if not self.is_cif else []
        self._check_source_connections()
        self.groups = []
        self.output_atoms = []
        self.output_block = None

    @classmethod
    def from_string(cls, text, *, source_name='<memory>', **options):
        with StringIO(text) as source:
            source.name = source_name
            return cls(source, **options)

    @staticmethod
    def _hash(path):
        with open(path, "rb") as handle:
            return hashlib.sha256(handle.read()).hexdigest()

    def _source_warning(self, kind, message, **details):
        prefix = ('Source annotation notice: ' if kind in ('unreviewed_pdb_records', 'unreviewed_cif_category')
                  else 'Non-standard source: ')
        warning = prefix + message
        issue = dict(kind=kind, source=self.path, source_sha256=self.report['source_sha256'],
                     warning=warning, **details)
        issues = self.report.setdefault('source_issues', [])
        if not any(old['kind'] == kind and old['warning'] == warning for old in issues):
            issues.append(issue)
        if warning not in self.report['warnings']:
            self.report['warnings'].append(warning)

    def unrecognized_cif_categories(self):
        return {category: _cif_rows(self.block, category)
                for category in self.block.get_mmcif_category_names()
                if category.strip('_.') not in self.SAFE_CATEGORIES | self.CHAIN_CATEGORIES | self.SPECIAL_CATEGORIES
                and not category.strip('_.').startswith(self.INVALIDATED_PREFIXES)}

    def _check_source_connections(self):
        rows = _cif_rows(self.block, '_struct_conn.')
        index = _AtomIndex(self.atoms)
        kept, missing = [], []
        for row in rows:
            (kept if all(index.exists(row, prefix) for prefix in ('ptnr1_', 'ptnr2_'))
             else missing).append(row)
        if missing:
            if self.metadata_policy == 'error':
                raise StructureEditError('Dangling struct_conn atom endpoint in source.')
            self._source_warning('dangling_struct_conn',
                f'{len(missing)} connection annotation(s) refer to missing atoms. '
                'Ignored these annotations in the working copy; atom coordinates and valid connections are retained.',
                rows=missing)
            _set_cif_rows(self.block, '_struct_conn.', kept)

    def _check_pdb_atoms(self, st):
        seen, anisou, model = set(), [], "1"
        serial_lines, identities = defaultdict(list), set()
        atom_count = 0
        for line_number, line in enumerate(self.pdb_lines, 1):
            rec = line[:6].strip()
            if rec == "MODEL":
                model = line[10:14].strip()
            elif rec in ("ATOM", "HETATM", "ANISOU"):
                if len(line) < 54:
                    raise StructureEditError("Truncated PDB coordinate/ANISOU record.")
                if rec == "ANISOU":
                    anisou.append((model, line[6:11], line[12:27]))
                    continue
                key = (model, line[6:11])
                if key in seen and self.metadata_policy == 'error':
                    raise StructureEditError("Duplicate atom serial within PDB MODEL: " + str(key))
                seen.add(key)
                serial_lines[key].append(line_number)
                identities.add((model, line[6:11], line[12:27]))
                atom_count += 1
                for start in (30, 38, 46):
                    if not math.isfinite(float(line[start:start + 8])):
                        raise StructureEditError("Nonfinite PDB coordinates.")
        if atom_count != sum(m.count_atom_sites() for m in st):
            raise StructureEditError("PDB parser did not retain every atom; resolve input ambiguity first.")
        if any(identity not in identities for identity in anisou):
            raise StructureEditError("ANISOU does not match its atom serial and identity.")
        duplicates = [{'model': m, 'serial': s.strip(), 'lines': lines}
                      for (m, s), lines in serial_lines.items() if len(lines) > 1]
        if duplicates:
            self._source_warning('duplicate_pdb_serials',
                f"{sum(len(item['lines']) for item in duplicates)} atoms share repeated PDB serial numbers. "
                'Atom identities and coordinates are retained; serials are regenerated in the working copy. '
                'Ambiguous CONECT references are ignored.', serials=duplicates)

    def _read_conect(self, st):
        # Read fixed-width fields, including adjacent five-digit serials.
        pairs = Counter()
        pair_records = defaultdict(list)
        for line_number, line in enumerate(self.pdb_lines, 1):
            if not line.startswith("CONECT"):
                continue
            try:
                nums = [int(line[i:i + 5]) for i in range(6, min(len(line), 31), 5) if line[i:i + 5].strip()]
            except ValueError as exc:
                raise StructureEditError("CONECT with nondecimal serials is not supported; convert to mmCIF first.") from exc
            if len(nums) > 1:
                pairs.update((nums[0], target) for target in nums[1:])
                for target in set(nums[1:]):
                    pair_records[(nums[0], target)].append({'line': line_number, 'text': line})
        serial_to_id = {}
        for m in st:
            for c in m:
                for res in c:
                    for atom in res:
                        serial_to_id.setdefault(atom.serial, []).append(
                            (str(m.num), c.name, str(res.seqid.num), _null_id(res.seqid.icode.strip()),
                             res.name, atom.name, atom.altloc.replace("\x00", "")))
        ids = {_atom_key(r): r["id"] for r in self.atoms}
        bonds, missing, ambiguous = {}, [], []
        for (a, b), order in pairs.items():
            common_models = ({key[0] for key in serial_to_id.get(a, ())} &
                             {key[0] for key in serial_to_id.get(b, ())})
            if not common_models:
                if self.metadata_policy == 'error':
                    raise StructureEditError(f"Dangling CONECT reference {a} -> {b} in source.")
                missing.append({'from': a, 'to': b, 'order': order, 'records': pair_records[(a, b)]})
                continue
            ambiguous_models = {m for m in common_models
                                if sum(key[0] == m for key in serial_to_id[a]) != 1
                                or sum(key[0] == m for key in serial_to_id[b]) != 1}
            if ambiguous_models:
                ambiguous.append({'from': a, 'to': b, 'models': sorted(ambiguous_models),
                                  'order': order, 'records': pair_records[(a, b)]})
            for ka in serial_to_id[a]:
                if ka[0] in ambiguous_models:
                    continue
                for kb in serial_to_id[b]:
                    if ka[0] == kb[0]:
                        key = tuple(sorted((ids[ka], ids[kb]), key=int))
                        bonds[key] = max(bonds.get(key, 0), order)
        if missing:
            examples = ', '.join(f"{ref['from']} -> {ref['to']}" for ref in missing[:5])
            if len(missing) > 5:
                examples += ', ...'
            self._source_warning('dangling_conect',
                f'{len(missing)} dangling CONECT reference(s) ({examples}). '
                'Ignored these references in the working copy because the atoms do not exist in the same model; '
                'atom coordinates and valid bonds are retained.', references=missing)
        if ambiguous:
            self._source_warning('ambiguous_conect',
                f'{len(ambiguous)} CONECT reference(s) use repeated atom serials. '
                'Ignored these references in affected models instead of guessing bond endpoints.',
                references=ambiguous)
        return [(a, b, n) for (a, b), n in bonds.items()]

    @staticmethod
    def _validate_atoms(atoms):
        identities, serials = set(), set()
        label_entities, label_auth = {}, {}
        for row in atoms:
            key = _atom_key(row)
            if key in identities:
                raise StructureEditError("Duplicate atom identity (including model, insertion code and altloc): " + str(key))
            identities.add(key)
            if row["id"] in serials:
                raise StructureEditError("atom_site.id must be unique across the mmCIF block.")
            serials.add(row["id"])
            for field in ("Cartn_x", "Cartn_y", "Cartn_z", "occupancy", "B_iso_or_equiv"):
                if _present(row.get(field)) and not math.isfinite(float(row[field])):
                    raise StructureEditError("Nonfinite atom value in " + field)
            if _present(row.get("occupancy")) and not 0 <= float(row["occupancy"]) <= 1.00001:
                raise StructureEditError("Occupancy outside [0, 1].")
            label = row["label_asym_id"]
            if not _present(label):
                raise StructureEditError("Every atom needs a label_asym_id.")
            auth = row["auth_asym_id"]
            if label in label_auth and label_auth[label] != auth:
                raise StructureEditError("One label_asym_id refers to multiple author chains.")
            label_auth[label] = auth
            entity = row.get("label_entity_id")
            if label in label_entities and label_entities[label] != entity:
                raise StructureEditError("One label_asym_id refers to multiple entities.")
            label_entities[label] = entity

    def ca_chains(self):
        # Choose one CA conformer per residue by occupancy, using the first model.
        _, np = _structure_dependencies()
        candidates = defaultdict(dict)
        for row in self.atoms:
            if str(row["pdbx_PDB_model_num"]) != self.models[0]:
                continue
            if row["label_asym_id"] not in self.poly_labels or row["label_atom_id"] != "CA" or row.get("type_symbol", "C").upper() != "C":
                continue
            key = _residue_key(row)[1:]
            chain = _null_id(row["auth_asym_id"])
            old = candidates[chain].get(key)
            rank = (float(row.get("occupancy") or 0), _null_id(row.get("label_alt_id")) in ("", "A"))
            if old is None or rank > old[0]:
                candidates[chain][key] = (rank, row)
        return {c: [{"id": (k[0] + k[1]), "key": k,
                      "coord": np.array([float(r[f]) for f in ("Cartn_x", "Cartn_y", "Cartn_z")])}
                    for k, (_, r) in residues.items()] for c, residues in candidates.items()}

    def preview(self):
        waters, others = [], []
        for r in self.atoms:
            if str(r["pdbx_PDB_model_num"]) != self.models[0] or r["label_asym_id"] in self.poly_labels:
                continue
            xyz = tuple(float(r[f]) for f in ("Cartn_x", "Cartn_y", "Cartn_z"))
            if r["label_comp_id"] in ("HOH", "WAT", "DOD"):
                if r["type_symbol"].upper() == "O":
                    waters.append(xyz)
            else:
                others.append(xyz)
        return self.ca_chains(), waters, others

    def associated_residues(self, core_chains, axial_limit=4.0, axis=None):
        # Assign complete nonpolymer residues once, to the nearest protein chain.
        _, np = _structure_dependencies()
        axis = np.asarray(axis if axis is not None else (0, 0, 1), dtype=float)
        axis /= np.linalg.norm(axis)
        ca = self.ca_chains()
        centers = {c: np.mean([r["coord"] for r in ca[c]], axis=0) for c in core_chains if c in ca}
        residues = defaultdict(list)
        for r in self.atoms:
            if str(r["pdbx_PDB_model_num"]) == self.models[0] and r["label_asym_id"] not in self.poly_labels:
                residues[_residue_key(r)].append(r)
        owners = {}
        for key, rows in residues.items():
            if key[0] in core_chains:
                owners[key] = key[0]
                continue
            pos = np.mean([[float(r[f]) for f in ("Cartn_x", "Cartn_y", "Cartn_z")] for r in rows], axis=0)
            choices = [(float(np.linalg.norm(pos - cen)), c) for c, cen in centers.items()
                       if abs(float(np.dot(pos - cen, axis))) <= axial_limit]
            if choices:
                owners[key] = min(choices)[1]
        return owners

    def trim(self, keep_chains, *, rename=None, keep_residues=()):
        keep = set(keep_chains)
        unknown = keep - set(self.chains)
        if unknown or not keep:
            raise StructureEditError("Empty selection or unknown chains: " + ", ".join(sorted(unknown)))
        residues = set(keep_residues)
        rows = [r for r in self.atoms if _null_id(r["auth_asym_id"]) in keep or _residue_key(r) in residues]
        mapping = {c: (rename or {}).get(c, c) for c in self.chains if any(_null_id(r["auth_asym_id"]) == c for r in rows)}
        self.report["operation"] = "trim/rename"
        return self.apply_copies([{"chains": mapping, "atom_ids": {r["id"] for r in rows}}])

    def rename(self, mapping):
        if set(mapping) - set(self.chains):
            raise StructureEditError("Rename map contains unknown author chain IDs.")
        self.report["operation"] = "rename"
        return self.apply_copies([{"chains": {c: mapping.get(c, c) for c in self.chains}}])

    def select_model(self, model_num, output_model="1"):
        # Select a source ensemble member before editing or merging live coordinates.
        model_num, output_model = str(model_num), str(output_model)
        if self.output_block is not None or model_num not in self.models:
            raise StructureEditError('Select an existing source model before editing.')
        original_models = list(self.models)
        self.atoms = [r for r in self.atoms if str(r['pdbx_PDB_model_num']) == model_num]
        atom_ids = {r['id'] for r in self.atoms}
        self.pdb_bonds = [(a, b, order) for a, b, order in self.pdb_bonds if a in atom_ids and b in atom_ids]
        unrecognized = self.unrecognized_cif_categories() if self.metadata_policy == 'preserve' else {}
        for category in self.block.get_mmcif_category_names():
            if category in unrecognized:
                continue
            rows = _cif_rows(self.block, category)
            if not rows:
                continue
            fields = [k for k in rows[0] if k.lower() in ('pdbx_pdb_model_num', 'pdb_model_num')]
            if not fields and category != '_atom_site_anisotrop.':
                continue
            kept = []
            for row in rows:
                if category == '_atom_site_anisotrop.' and row['id'] not in atom_ids:
                    continue
                if any(_present(row.get(k)) and str(row[k]) != model_num for k in fields):
                    continue
                kept.append(dict(row, **{k: output_model for k in fields if _present(row.get(k))}))
            _set_cif_rows(self.block, category, kept)
        for row in self.atoms:
            row['pdbx_PDB_model_num'] = output_model
        self.models = [output_model]
        self.chains = list(dict.fromkeys(_null_id(r['auth_asym_id']) for r in self.atoms))
        self.report['ensemble_selection'] = {'source_models': original_models,
            'selected_source_model': model_num, 'working_model': output_model}
        self.report['warnings'].append(f'Working copy contains selected ensemble member {model_num} of {len(original_models)}; the other source members were not combined into layers.')
        return self

    def apply_copies(self, copies):
        gemmi, np = _structure_dependencies()
        self.report.setdefault("operation", "explicit rigid copies")
        self.groups, self.output_atoms = [], []
        used_auth = set()
        next_label = 0
        for number, spec in enumerate(copies):
            mapping = dict(spec["chains"])
            if set(mapping) - set(self.chains):
                raise StructureEditError("Copy refers to an unknown source chain.")
            if len(set(mapping.values())) != len(mapping) or used_auth.intersection(mapping.values()):
                raise StructureEditError("Output author chain IDs must be unique; merging chains is not a rename.")
            if any(not isinstance(c, str) or not c or not c.isascii() or not c.isalnum() for c in mapping.values()):
                raise StructureEditError("Output chain IDs must be nonempty ASCII alphanumeric strings.")
            used_auth.update(mapping.values())
            rotation = np.asarray(spec.get("matrix", np.eye(3)), dtype=float)
            translation = np.asarray(spec.get("vector", np.zeros(3)), dtype=float)
            if rotation.shape != (3, 3) or translation.shape != (3,) or not np.isfinite(rotation).all() or not np.isfinite(translation).all():
                raise StructureEditError("Rigid transform must contain finite 3x3 R and length-3 t.")
            if not np.allclose(rotation.T @ rotation, np.eye(3), atol=1e-6) or not np.isclose(np.linalg.det(rotation), 1, atol=1e-6):
                raise StructureEditError("Transform must be a proper rotation (det R = +1), without scale/reflection.")
            selection = spec.get("atom_ids")
            rows = [r for r in self.atoms if _null_id(r["auth_asym_id"]) in mapping and (selection is None or r["id"] in selection)]
            if not rows:
                raise StructureEditError("A requested copy contains no atoms.")
            labels, ids, output_rows = {}, {}, []
            rotated_coordinates = not np.allclose(rotation, np.eye(3), atol=1e-12, rtol=0)
            identity_coordinates = np.array_equal(rotation, np.eye(3)) and np.array_equal(translation, np.zeros(3))
            for row in rows:
                old = row["label_asym_id"]
                if old not in labels:
                    # Ligands sharing a polymer author ID still need independent label/entity IDs.
                    labels[old] = structure_chain_id(next_label)
                    next_label += 1
                new = dict(row)
                new["id"] = str(len(self.output_atoms) + 1)
                ids[row["id"]] = new["id"]
                new["auth_asym_id"] = mapping[_null_id(row["auth_asym_id"])]
                new["label_asym_id"] = labels[old]
                if not identity_coordinates:
                    xyz = rotation @ np.array([float(row[f]) for f in ("Cartn_x", "Cartn_y", "Cartn_z")]) + translation
                    for field, value in zip(("Cartn_x", "Cartn_y", "Cartn_z"), xyz):
                        new[field] = f"{value:.6f}"
                # Drop frame-dependent uncertainties that cannot be rotated without covariances.
                uncertainty_fields = {"Cartn_x_esd", "Cartn_y_esd", "Cartn_z_esd"}
                invalidated = {field: value for field, value in row.items()
                               if _present(value) and (field.startswith('fract_') or
                                   (rotated_coordinates and field in uncertainty_fields))}
                if invalidated:
                    for field in invalidated:
                        new[field] = None
                    self.report.setdefault('atom_coordinate_metadata_archive', {}).setdefault(row['id'], {}).update(invalidated)
                    self.report['warnings'].append('Fractional coordinate metadata and rotated Cartesian ESDs were invalidated where present; source values are archived by atom ID.')
                self.output_atoms.append(new)
                output_rows.append(new)
            group = {"chains": mapping, "labels": labels, "ids": ids, "rotation": rotation,
                     "translation": translation, "rows": rows, "index": number,
                     "output_rows": output_rows, "atom_index": _AtomIndex(output_rows)}
            self.groups.append(group)
            self.report["chain_copies"].append({"copy": number, "author_chains": mapping,
                "label_chains": labels, "matrix": rotation.tolist(), "vector": translation.tolist()})
        if not self.output_atoms:
            raise StructureEditError("Refusing to write an empty structure.")
        self._validate_atoms(self.output_atoms)
        self.output_index = _AtomIndex(self.output_atoms)
        self.output_document = gemmi.cif.Document()
        self.output_block = self.output_document.add_new_block("AmyloidModifier")
        self._metadata()
        _set_cif_rows(self.output_block, "_atom_site.", self.output_atoms)
        self._anisotrop()
        self._conect_to_connections()
        self._collect_declared_bonds()
        self._validate_references()
        return self

    def _collect_declared_bonds(self):
        by_model = defaultdict(list)
        by_residue = defaultdict(list)
        for atom in self.output_atoms:
            by_model[str(atom['pdbx_PDB_model_num'])].append(atom)
            if atom.get('group_PDB') == 'HETATM':
                by_residue[(str(atom['pdbx_PDB_model_num']), _residue_key(atom))].append(atom)
        model_indexes = [_AtomIndex(atoms) for atoms in by_model.values()]
        def order(row):
            return {'sing':1,'doub':2,'trip':3}.get(str(row.get('pdbx_value_order', row.get('value_order','sing'))).lower(),1)
        for connection in _cif_rows(self.output_block, '_struct_conn.'):
            if connection.get('conn_type_id') not in ('covale','disulf','modres','metalc'):
                continue
            for index in model_indexes:
                partners = [index.matches(connection, prefix) for prefix in ('ptnr1_','ptnr2_')]
                for a in partners[0]:
                    for b in partners[1]:
                        alt_a, alt_b = _null_id(a.get('label_alt_id')), _null_id(b.get('label_alt_id'))
                        if alt_a and alt_b and alt_a != alt_b:
                            continue
                        self.output_bonds.append((a['id'],b['id'],order(connection)))
        bonds = defaultdict(list)
        for row in _cif_rows(self.output_block, '_chem_comp_bond.'):
            bonds[row['comp_id']].append(row)
        for (_, residue), atoms in by_residue.items():
            for row in bonds.get(residue[-1], []):
                a_atoms = [a for a in atoms if a['label_atom_id'] == row['atom_id_1']]
                b_atoms = [a for a in atoms if a['label_atom_id'] == row['atom_id_2']]
                for a in a_atoms:
                    for b in b_atoms:
                        aa, ab = _null_id(a.get('label_alt_id')), _null_id(b.get('label_alt_id'))
                        if aa and ab and aa != ab:
                            continue
                        self.output_bonds.append((a['id'],b['id'],order(row)))
        unique = {}
        for a,b,n in self.output_bonds:
            key = tuple(sorted((a,b),key=int))
            unique[key] = max(n,unique.get(key,0))
        self.output_bonds = [(a,b,n) for (a,b),n in unique.items()]

    def _audit(self, category, action, reason="", content=None):
        item = {"action": action, "reason": reason}
        if content is not None:
            item["source_content"] = content
        self.report["metadata"][category] = item

    def _remap_row(self, row, group, *, extra_label=(), extra_auth=()):
        out = dict(row)
        for key, value in row.items():
            lower = key.lower()
            label = ("asym_id" in lower and "auth" not in lower) or key in extra_label
            auth = ("auth_asym_id" in lower) or key in extra_auth or lower in ("pdbx_strand_id", "pdb_strand_id", "pdbx_pdb_strand_id")
            if not (label or auth) or not _present(value):
                continue
            mapping = group["chains"] if auth else group["labels"]
            # Chain-list fields have list semantics; partners/endpoints do not.
            if lower in ("pdbx_strand_id", "pdb_strand_id", "pdbx_pdb_strand_id"):
                values = [mapping[c.strip()] for c in str(value).split(",") if c.strip() in mapping]
                if not values:
                    return None
                out[key] = ",".join(dict.fromkeys(values))
            elif value not in mapping:
                return None
            else:
                out[key] = mapping[value]
        return out

    @staticmethod
    def _endpoint_exists(row, atoms, prefix):
        if isinstance(atoms, _AtomIndex):
            return atoms.exists(row, prefix)
        tests = _endpoint_tests(row, prefix)
        return bool(tests) and any(all(str(a.get(k)) == v for k, v in tests) for a in atoms)

    def _metadata(self):
        out = self.output_block
        special = {"atom_site", "atom_site_anisotrop", "struct_asym", "struct_sheet", "struct_sheet_range",
                   "struct_sheet_order", "pdbx_struct_sheet_hbond", "struct_site"}
        live_entities = {r.get("label_entity_id") for r in self.output_atoms if _present(r.get("label_entity_id"))}
        align_ids, site_ids = {}, {}
        for group in self.groups:
            for row in _cif_rows(self.block, '_struct_ref_seq.'):
                if self._remap_row(row, group) is not None:
                    align_ids[group['index'], row['align_id']] = str(len(align_ids) + 1)
            for row in _cif_rows(self.block, '_struct_site_gen.'):
                if self._remap_row(row,group) is not None:
                    key = (group['index'],row['site_id'])
                    site_ids.setdefault(key, str(len(site_ids)+1))
        for category in self.block.get_mmcif_category_names():
            name = category.strip("_.")
            rows = _cif_rows(self.block, category)
            if not rows or name in special:
                continue
            if name == 'amyloid_unrecognized_pdb':
                if self.metadata_policy == 'preserve':
                    _set_cif_rows(out, category, rows)
                continue
            if name == 'audit_syntax':
                self._audit(category, 'invalidated',
                            'Source serializer column-layout hints do not describe the rewritten CIF.', rows)
                continue
            if name.startswith(self.INVALIDATED_PREFIXES):
                self._audit(category, "invalidated", "Describes the source assembly, frame, refinement, archive history or validation; recompute for this derivative.", rows)
                continue
            if name in self.CHAIN_CATEGORIES:
                mapped = []
                for group in self.groups:
                    group_atoms = group["output_rows"]
                    group_index = group["atom_index"]
                    for row in rows:
                        new = self._remap_row(row, group)
                        if new is None:
                            continue
                        if name in ('struct_ref_seq','struct_ref_seq_dif') and 'align_id' in new:
                            new_align = align_ids.get((group['index'], row['align_id']))
                            if new_align is None:
                                continue
                            new['align_id'] = new_align
                        if name == 'struct_site_gen':
                            new['site_id'] = site_ids[group['index'],row['site_id']]
                        if name in ('pdbx_nonpoly_scheme','pdbx_branch_scheme'):
                            seq = new.get('pdb_seq_num',new.get('auth_seq_num'))
                            component = new.get('mon_id',new.get('pdb_mon_id'))
                            if _present(seq) and not any(a['label_asym_id'] == new.get('asym_id') and
                                a['auth_seq_id'] == seq and a['label_comp_id'] == component and
                                _null_id(a.get('pdbx_PDB_ins_code')) == _null_id(new.get('pdb_ins_code')) for a in group_atoms):
                                continue
                        if name == "struct_conn":
                            if any(_present(row.get(p + "symmetry")) and row[p + "symmetry"] not in ("1_555", "1555") for p in ("ptnr1_", "ptnr2_")):
                                if row.get('conn_type_id') in ('disulf','covale','modres'):
                                    raise StructureEditError('A retained covalent connection uses crystal symmetry. Resolve the symmetry-related atoms explicitly before editing.')
                                self._audit(category + "symmetry", "invalidated", "Nonidentity crystal contacts cannot be copied to a finite fibril.", rows)
                                continue
                            if not all(group_index.exists(new, p) for p in ("ptnr1_", "ptnr2_")):
                                continue
                        if name == "struct_conf" and not all(group_index.exists(new, p) for p in ("beg_", "end_")):
                            continue
                        if name == "struct_site_gen" and not group_index.exists(new, ""):
                            continue
                        if name == "struct_mon_prot_cis" and not all(group_index.exists(new, p) for p in ("", "pdbx_")):
                            continue
                        for key in ("id", "pdbx_id"):
                            if key in new:
                                new[key] = str(len(mapped) + 1)
                        if name == 'struct_conf' and 'pdbx_PDB_helix_id' in new:
                            new['pdbx_PDB_helix_id'] = str(len(mapped)+1)
                        mapped.append(new)
                _set_cif_rows(out, category, mapped)
                self._audit(category, "remapped", "All partners must survive in the same rigid copy.")
            elif name in self.SAFE_CATEGORIES:
                kept = []
                for row in rows:
                    row = dict(row)
                    entity = row.get("entity_id", row.get("pdbx_entity_id", row.get("id") if name == "entity" else None))
                    if _present(entity) and name not in ("em_entity_assembly",) and entity not in live_entities:
                        continue
                    if name == "entity_poly" and "pdbx_strand_id" in row:
                        row["pdbx_strand_id"] = ",".join(dict.fromkeys(a["auth_asym_id"] for a in self.output_atoms if a.get("label_entity_id") == entity))
                    if name == "entity" and "pdbx_number_of_molecules" in row:
                        members = [a for a in self.output_atoms if a.get('label_entity_id') == entity and str(a['pdbx_PDB_model_num']) == self.models[0]]
                        poly_out = {g['labels'][old] for g in self.groups for old in self.poly_labels if old in g['labels']}
                        row['pdbx_number_of_molecules'] = str(len({a['label_asym_id'] if a['label_asym_id'] in poly_out else _residue_key(a) for a in members}))
                    if name == 'em_software' and _present(row.get('fitting_id')):
                        self._audit('_em_software.fitting_id', 'invalidated', 'Removed source-model fitting software rows with invalidated fitting records.', rows)
                        continue
                    if name == 'em_software':
                        row.pop('fitting_id', None)
                    kept.append(row)
                _set_cif_rows(out, category, kept)
                self._audit(category, "preserved/filtered", "Source descriptive/experimental metadata; entity membership/counts updated.")
            else:
                if self.metadata_policy == "error":
                    raise StructureEditError(f"No reviewed edit policy for {category}. Use metadata_policy='audit' to remove it with a full export report, or add a tested handler.")
                self._source_warning('unreviewed_cif_category',
                    (f'Unreviewed mmCIF category {category} is preserved unchanged when saving; '
                     'its references may describe the original structure.' if self.metadata_policy == 'preserve' else
                     f'Unreviewed mmCIF category {category} is archived in memory and omitted from the working copy; atom coordinates are retained.'),
                    category=category, rows=rows)
                if self.metadata_policy == 'preserve':
                    _set_cif_rows(out, category, rows)
                self._audit(category, 'preserved/unmapped' if self.metadata_policy == 'preserve' else 'removed',
                            'Unreviewed category; cannot promise its references remain valid.', rows)
        asym = []
        seen = set()
        for row in self.output_atoms:
            label = row["label_asym_id"]
            if label not in seen:
                asym.append({"id": label, "entity_id": row.get("label_entity_id")})
                seen.add(label)
        _set_cif_rows(out, "_struct_asym.", asym)
        self._audit('_struct_asym.', 'rebuilt', 'Independent component identifiers for every retained/generated copy.')
        self._sequence_scheme()
        self._sheets()
        sites = _cif_rows(out, "_struct_site_gen.")
        site_counts = Counter(r["site_id"] for r in sites)
        site_headers = []
        for (group_index, old_id), new_id in site_ids.items():
            group = self.groups[group_index]
            for row in _cif_rows(self.block, '_struct_site.'):
                if row['id'] != old_id or new_id not in site_counts:
                    continue
                mapped = self._remap_row(row, group)
                if mapped is None:
                    continue
                if _present(mapped.get('pdbx_auth_asym_id')) and not group['atom_index'].exists(mapped, 'pdbx_'):
                    continue
                site_headers.append(dict(mapped, id=new_id, pdbx_num_residues=str(site_counts[new_id])))
        retained_sites = {r['id'] for r in site_headers}
        _set_cif_rows(out, '_struct_site.', site_headers)
        _set_cif_rows(out, '_struct_site_gen.', [r for r in sites if r['site_id'] in retained_sites])
        self._audit('_struct_site.', 'remapped/filtered',
                    'Site-defining residue/ligand and member references must survive in the same copy.',
                    _cif_rows(self.block, '_struct_site.'))
        self.report['source_entry_id'] = self.block.find_value('_entry.id')
        for category in out.get_mmcif_category_names():
            if self.report['metadata'].get(category, {}).get('action') == 'preserved/unmapped':
                continue
            rows = _cif_rows(out, category)
            changed = False
            for row in rows:
                for key in row:
                    if key == 'entry_id' or (category == '_entry.' and key == 'id'):
                        row[key] = 'AMYLOID_MODIFIED'
                        changed = True
            if changed:
                _set_cif_rows(out, category, rows)

    def _sequence_scheme(self):
        if _cif_rows(self.output_block, '_pdbx_poly_seq_scheme.'):
            return
        sequence = defaultdict(list)
        for row in _cif_rows(self.output_block, '_entity_poly_seq.'):
            sequence[row['entity_id']].append(row)
        by_label = defaultdict(list)
        for row in self.output_atoms:
            if str(row['pdbx_PDB_model_num']) == self.models[0]:
                by_label[row['label_asym_id']].append(row)
        scheme = []
        for label, atoms in by_label.items():
            entity = atoms[0].get('label_entity_id')
            if entity not in sequence:
                continue
            observed = {}
            for atom in atoms:
                observed.setdefault((atom.get('label_seq_id'), atom['label_comp_id']), atom)
            for residue in sequence[entity]:
                atom = observed.get((residue['num'], residue['mon_id']))
                scheme.append({'asym_id':label, 'entity_id':entity, 'seq_id':residue['num'],
                    'mon_id':residue['mon_id'], 'ndb_seq_num':residue['num'],
                    'pdb_seq_num':atom.get('auth_seq_id') if atom else None,
                    'auth_seq_num':atom.get('auth_seq_id') if atom else None,
                    'pdb_mon_id':atom['label_comp_id'] if atom else None,
                    'auth_mon_id':atom.get('auth_comp_id', atom['label_comp_id']) if atom else None,
                    'pdb_strand_id':atoms[0]['auth_asym_id'],
                    'pdb_ins_code':(_null_id(atom.get('pdbx_PDB_ins_code')) or False) if atom else False,
                    'hetero':residue.get('hetero','n')})
        if scheme:
            _set_cif_rows(self.output_block, '_pdbx_poly_seq_scheme.', scheme)
            self._audit('_pdbx_poly_seq_scheme.', 'rebuilt', 'Full source sequence with observed author numbers/insertion codes; missing coordinates remain unknown.')

    def _sheets(self):
        ranges = _cif_rows(self.block, "_struct_sheet_range.")
        orders = _cif_rows(self.block, "_struct_sheet_order.")
        hbonds = _cif_rows(self.block, "_pdbx_struct_sheet_hbond.")
        sheets_out, ranges_out, order_out, hbond_out = [], [], [], []
        self.sheet_origins = {}
        for group in self.groups:
            atoms = group["atom_index"]
            sheet_ids = list(dict.fromkeys(r["sheet_id"] for r in ranges))
            for sheet in sheet_ids:
                surviving = []
                for r in ranges:
                    if r["sheet_id"] == sheet:
                        mapped = self._remap_row(r, group)
                        if mapped and all(self._endpoint_exists(mapped, atoms, p) for p in ("beg_", "end_")):
                            surviving.append(mapped)
                if not surviving:
                    continue
                new_sheet = str(len(sheets_out) + 1)
                self.sheet_origins[new_sheet] = sheet
                range_map = {r["id"]: str(i + 1) for i, r in enumerate(surviving)}
                sheets_out.append({"id": new_sheet, "number_strands": str(len(surviving))})
                for row in surviving:
                    ranges_out.append(dict(row, sheet_id=new_sheet, id=range_map[row["id"]]))
                for source, dest in ((orders, order_out), (hbonds, hbond_out)):
                    for row in source:
                        if row["sheet_id"] != sheet or any(row.get(k) not in range_map for k in ("range_id_1", "range_id_2")):
                            continue
                        mapped = self._remap_row(row, group)
                        if mapped is None:
                            continue
                        if source is hbonds and not all(self._endpoint_exists(mapped, atoms, p) for p in ("range_1_", "range_2_")):
                            continue
                        mapped.update(sheet_id=new_sheet, range_id_1=range_map[row["range_id_1"]], range_id_2=range_map[row["range_id_2"]])
                        dest.append(mapped)
        for cat, rows in (("struct_sheet", sheets_out), ("struct_sheet_range", ranges_out),
                          ("struct_sheet_order", order_out), ("pdbx_struct_sheet_hbond", hbond_out)):
            _set_cif_rows(self.output_block, "_" + cat + ".", rows)
            self._audit("_" + cat + ".", "remapped", "Surviving strand endpoints/registrations only; no new interlayer hydrogen bonds inferred.")

    def propagate_layer_sheets(self, sandwiches, instances, layer_step=1):
        # Extend only annotated sheet relationships; do not infer missing registrations.
        _, np = _structure_dependencies()
        locations = {c:(pf,l) for pf,s in enumerate(sandwiches) for l,c in enumerate(s)}
        ranges = _cif_rows(self.block, '_struct_sheet_range.')
        orders = _cif_rows(self.block, '_struct_sheet_order.')
        hbonds = _cif_rows(self.block, '_pdbx_struct_sheet_hbond.')
        first_atoms = [a for a in self.output_atoms if str(a['pdbx_PDB_model_num']) == self.models[0]]
        first_index = _AtomIndex(first_atoms)
        source_index = _AtomIndex([a for a in self.atoms if str(a['pdbx_PDB_model_num']) == self.models[0]])
        labels_by_chain = defaultdict(set)
        for a in first_atoms:
            if _present(a.get('label_seq_id')):
                labels_by_chain[a['auth_asym_id']].add(a['label_asym_id'])
        source_labels = defaultdict(set)
        for atom in self.atoms:
            if atom['label_asym_id'] in self.poly_labels:
                source_labels[atom['auth_asym_id']].add(atom['label_asym_id'])
        def mapped(row, chain_mapping):
            labels = {}
            for chain, destination in chain_mapping.items():
                if source_labels[chain]:
                    choices = labels_by_chain[destination]
                    if len(choices) != 1:
                        raise StructureEditError('Sheet propagation needs one polymer component per author chain.')
                    labels.update((label, next(iter(choices))) for label in source_labels[chain])
            return self._remap_row(row, {'chains':chain_mapping,'labels':labels})
        changed_sheets, new_components = set(), []
        copy_origins = {new:old for group in self.groups for old,new in group['chains'].items()}
        for sheet_id in dict.fromkeys(r['sheet_id'] for r in ranges):
            sr = [r for r in ranges if r['sheet_id'] == sheet_id]
            so = [r for r in orders if r['sheet_id'] == sheet_id]
            if not so or any(r.get('beg_auth_asym_id') != r.get('end_auth_asym_id') or r.get('beg_auth_asym_id') not in locations for r in sr):
                continue
            if len({locations[r['beg_auth_asym_id']][1] for r in sr}) < 2:
                continue
            # Match sheet intervals by ordinal position to allow differing strand endpoints.
            per_chain = defaultdict(list)
            for r in sr:
                per_chain[r['beg_auth_asym_id']].append(r)
            per_pf = defaultdict(set)
            for chain, chain_ranges in per_chain.items():
                per_pf[locations[chain][0]].add(len(chain_ranges))
            if any(len(counts) != 1 for counts in per_pf.values()):
                self.report['warnings'].append(f'Sheet {sheet_id} has inconsistent strand counts across layers; retained only directly mappable annotations.')
                continue
            slots = {}
            for chain, chain_ranges in per_chain.items():
                for slot,r in enumerate(sorted(chain_ranges,key=lambda r:(int(r['beg_auth_seq_id']),int(r['end_auth_seq_id'])))):
                    slots[r['id']] = slot
            nodes, source_nodes = {}, {}
            for row in sr:
                chain = row['beg_auth_asym_id']
                pf, layer = locations[chain]
                pattern = slots[row['id']]
                source_nodes[row['id']] = (pf, layer, pattern, chain)
                for (dest_pf, dest_layer), dest_chain in instances.items():
                    if dest_pf != pf or (dest_layer-layer) % layer_step or copy_origins.get(dest_chain) != chain:
                        continue
                    node = (pf, dest_layer, pattern)
                    out = mapped(row, {chain:dest_chain})
                    if out and all(first_index.exists(out, p) for p in ('beg_','end_')):
                        nodes.setdefault(node, out)
            edges, bonds = {}, {}
            for row in so:
                a, b = source_nodes.get(row['range_id_1']), source_nodes.get(row['range_id_2'])
                if a is None or b is None:
                    raise StructureEditError('Source sheet order has a missing range.')
                for target in list(nodes):
                    if target[0] != a[0] or target[2] != a[2]:
                        continue
                    shift = target[1]-a[1]
                    other = (b[0], b[1]+shift, b[2])
                    if shift % layer_step or other not in nodes:
                        continue
                    edge = (target, other)
                    edges[edge] = row
                    # Retain source registrations when copied sheet edges overlap.
                    for h in hbonds:
                        if h['sheet_id'] != sheet_id or h['range_id_1'] != row['range_id_1'] or h['range_id_2'] != row['range_id_2']:
                            continue
                        chain_map = {a[3]:instances[target[:2]], b[3]:instances[other[:2]]}
                        out = mapped(h,chain_map)
                        if not out:
                            continue
                        old_ends = [source_index.matches(h, p) for p in ('range_1_','range_2_')]
                        new_ends = [first_index.matches(out, p) for p in ('range_1_','range_2_')]
                        if any(not x for x in old_ends + new_ends):
                            continue
                        def length(pair):
                            return float(np.linalg.norm([float(pair[0][0][k])-float(pair[1][0][k]) for k in ('Cartn_x','Cartn_y','Cartn_z')]))
                        error = abs(length(old_ends)-length(new_ends))
                        if error <= .5 and (edge not in bonds or shift == 0):
                            bonds[edge] = out
            if not edges:
                continue
            changed_sheets.add(sheet_id)
            neighbors = {n:set() for n in nodes}
            for a,b in edges:
                neighbors[a].add(b); neighbors[b].add(a)
            remaining = set(nodes)
            while remaining:
                seed = min(remaining)
                component, pending = set(), [seed]
                while pending:
                    node = pending.pop()
                    if node in component:
                        continue
                    component.add(node); remaining.discard(node)
                    pending.extend(neighbors[node]-component)
                # A directed path determines strand ordering for PDB as well.
                incoming = {n:0 for n in component}
                for a,b in edges:
                    if a in component and b in component: incoming[b] += 1
                ordered, todo = [], sorted(n for n in component if incoming[n] == 0)
                while todo:
                    node = todo.pop(0); ordered.append(node)
                    for a,b in edges:
                        if a == node and b in component:
                            incoming[b] -= 1
                            if incoming[b] == 0: todo.append(b)
                if len(ordered) != len(component):
                    ordered = sorted(component)  # mmCIF can represent cyclic graphs.
                new_components.append((ordered,nodes,edges,bonds))
        if not changed_sheets:
            return self
        categories = ('_struct_sheet.','_struct_sheet_range.','_struct_sheet_order.','_pdbx_struct_sheet_hbond.')
        retained_ids = {new for new, old in self.sheet_origins.items() if old not in changed_sheets}
        output = {cat:[r for r in _cif_rows(self.output_block,cat) if r.get('sheet_id',r.get('id')) in retained_ids] for cat in categories}
        count = max([int(r['id']) for r in output['_struct_sheet.']] or [0])
        for ordered,nodes,edges,bonds in new_components:
            count += 1; new_id = str(count)
            ids = {node:str(i+1) for i,node in enumerate(ordered)}
            output['_struct_sheet.'].append({'id':new_id,'number_strands':str(len(ordered))})
            for node in ordered:
                output['_struct_sheet_range.'].append(dict(nodes[node],sheet_id=new_id,id=ids[node]))
            for (a,b), row in edges.items():
                if a not in ids or b not in ids: continue
                ids_update = dict(sheet_id=new_id,range_id_1=ids[a],range_id_2=ids[b])
                output['_struct_sheet_order.'].append(dict(row,**ids_update))
                if (a,b) in bonds:
                    output['_pdbx_struct_sheet_hbond.'].append(dict(bonds[a,b],**ids_update))
        for category, rows in output.items():
            _set_cif_rows(self.output_block,category,rows)
            self._audit(category,'remapped/extended','Translated explicitly annotated sheet-order patterns across layers; registration endpoints and distances checked.')
        self.report['periodic_sheets_extended'] = len(changed_sheets)
        self._validate_references()
        return self

    def _anisotrop(self):
        _, np = _structure_dependencies()
        rows = _cif_rows(self.block, "_atom_site_anisotrop.")
        output = []
        for group in self.groups:
            for row in rows:
                if row["id"] not in group["ids"]:
                    continue
                new = self._remap_row(row, group)
                if new is None:
                    raise StructureEditError("Anisotropic atom identity disagrees with atom_site.")
                new["id"] = group["ids"][row["id"]]
                for prefix in ("U", "B"):
                    keys = [f"{prefix}[{i}][{j}]" for i, j in ((1, 1), (2, 2), (3, 3), (1, 2), (1, 3), (2, 3))]
                    if all(_present(row.get(k)) for k in keys) and not np.array_equal(group["rotation"], np.eye(3)):
                        u11, u22, u33, u12, u13, u23 = (float(row[k]) for k in keys)
                        tensor = np.array([[u11, u12, u13], [u12, u22, u23], [u13, u23, u33]])
                        rotated = group["rotation"] @ tensor @ group["rotation"].T
                        for key, (i, j) in zip(keys, ((0, 0), (1, 1), (2, 2), (0, 1), (0, 2), (1, 2))):
                            new[key] = f"{rotated[i, j]:.8f}"
                if not np.allclose(group["rotation"], np.eye(3)):
                    for key in list(new):
                        if "esd" in key.lower():
                            new[key] = None
                            self.report["warnings"].append("Anisotropic ESDs invalidated: full covariance is required to rotate uncertainties.")
                output.append(new)
        _set_cif_rows(self.output_block, "_atom_site_anisotrop.", output)

    def _conect_to_connections(self):
        rows = _cif_rows(self.output_block, "_struct_conn.")
        by_id = {r["id"]: r for r in self.output_atoms}
        first_index = _AtomIndex([a for a in self.output_atoms if str(a['pdbx_PDB_model_num']) == self.models[0]])
        existing_pairs = set()
        for row in rows:
            for a in first_index.matches(row, 'ptnr1_'):
                for b in first_index.matches(row, 'ptnr2_'):
                    existing_pairs.add(frozenset((a['id'], b['id'])))
        self.output_bonds = []
        for group in self.groups:
            for a, b, order in self.pdb_bonds:
                if a not in group["ids"] or b not in group["ids"]:
                    continue
                a, b = group["ids"][a], group["ids"][b]
                self.output_bonds.append((a, b, order))
                # struct_conn is model-independent; equivalent ensemble copies are represented once while PDB CONECT is written per model.
                if str(by_id[a]["pdbx_PDB_model_num"]) != self.models[0]:
                    continue
                # Merge duplicate bond annotations while preserving order and connection type.
                if frozenset((a, b)) in existing_pairs:
                    continue
                new = {"id": f"conect{len(rows) + 1}", "conn_type_id": "covale",
                       "pdbx_value_order": {1: "sing", 2: "doub", 3: "trip"}.get(order, "sing"),
                       "details": "Explicit source PDB CONECT; chemistry not inferred"}
                for prefix, atom in (("ptnr1_", by_id[a]), ("ptnr2_", by_id[b])):
                    for field in ("label_asym_id", "label_comp_id", "label_seq_id", "label_atom_id", "auth_asym_id", "auth_seq_id"):
                        new[prefix + field] = atom.get(field)
                    new["pdbx_" + prefix + "label_alt_id"] = atom.get("label_alt_id")
                    new["pdbx_" + prefix + "PDB_ins_code"] = atom.get("pdbx_PDB_ins_code")
                    new[prefix + "symmetry"] = "1_555"
                rows.append(new)
                # Unspecified altlocs remain wildcard endpoints, just as in the source matching rules, including subsequent CONECT records.
                for left in first_index.matches(new, 'ptnr1_'):
                    for right in first_index.matches(new, 'ptnr2_'):
                        existing_pairs.add(frozenset((left['id'], right['id'])))
        _set_cif_rows(self.output_block, "_struct_conn.", rows)
        if rows:
            types = {r["id"]: r for r in _cif_rows(self.output_block, "_struct_conn_type.")}
            for row in rows:
                types.setdefault(row["conn_type_id"], {"id": row["conn_type_id"]})
            _set_cif_rows(self.output_block, "_struct_conn_type.", list(types.values()))

    def propagate_layer_bonds(self, sandwiches, instances, layer_step=1, tolerance=0.20):
        # Propagate annotated covalent topology across repeat boundaries.
        _, np = _structure_dependencies()
        locations = {chain: (pf, layer) for pf, strand in enumerate(sandwiches) for layer, chain in enumerate(strand)}
        depth = len(sandwiches[0])
        source_by_id = {r['id']: r for r in self.atoms}
        source_first = _AtomIndex([r for r in self.atoms if str(r['pdbx_PDB_model_num']) == self.models[0]])
        output_by_key = {_atom_key(r): r for r in self.output_atoms}
        rows = _cif_rows(self.output_block, '_struct_conn.')
        templates = []
        for row in _cif_rows(self.block, '_struct_conn.'):
            kind = row.get('conn_type_id', '')
            if kind not in ('disulf', 'covale', 'modres'):
                continue
            if any(_present(row.get(p + 'symmetry')) and row[p + 'symmetry'] not in ('1_555','1555') for p in ('ptnr1_', 'ptnr2_')):
                raise StructureEditError('A covalent connection uses crystal symmetry. Resolve its symmetry-related atoms explicitly before helical expansion.')
            endpoints = []
            for prefix in ('ptnr1_', 'ptnr2_'):
                matches = source_first.matches(row, prefix)
                endpoints.append(matches)
            if any(not m for m in endpoints):
                raise StructureEditError('Source covalent connection has a missing atom endpoint.')
            for a in endpoints[0]:
                for b in endpoints[1]:
                    aa, ab = _null_id(a.get('label_alt_id')), _null_id(b.get('label_alt_id'))
                    if aa and ab and aa != ab:
                        continue
                    templates.append((a, b, row))
        for a, b, order in self.pdb_bonds:
            if str(source_by_id[a]['pdbx_PDB_model_num']) == self.models[0]:
                templates.append((source_by_id[a], source_by_id[b],
                                  {'conn_type_id':'covale', 'pdbx_value_order':{1:'sing',2:'doub',3:'trip'}.get(order, 'sing')}))
        def endpoint_key(row, prefix):
            return tuple(_null_id(row.get(prefix + k)) for k in ('auth_asym_id','auth_seq_id','label_comp_id','label_atom_id')) + (
                _null_id(row.get('pdbx_' + prefix + 'PDB_ins_code')), _null_id(row.get('pdbx_' + prefix + 'label_alt_id')))
        def pair_key(row):
            return tuple(sorted((endpoint_key(row, 'ptnr1_'), endpoint_key(row, 'ptnr2_'))))
        seen = {pair_key(r) for r in rows}
        made = 0
        all_layers = sorted({layer for pf, layer in instances})
        for a, b, template in templates:
            if a['auth_asym_id'] not in locations or b['auth_asym_id'] not in locations:
                continue  # Nonpolymer copies already follow their assigned owner.
            pfa, la = locations[a['auth_asym_id']]
            pfb, lb = locations[b['auth_asym_id']]
            for dest_a in all_layers:
                shift = dest_a - la
                dest_b = lb + shift
                if shift % layer_step or (pfa, dest_a) not in instances or (pfb, dest_b) not in instances:
                    continue
                if 0 <= dest_a < depth and 0 <= dest_b < depth:
                    continue
                chain_a, chain_b = instances[pfa, dest_a], instances[pfb, dest_b]
                model_pairs = []
                for model in self.models:
                    ka, kb = list(_atom_key(a)), list(_atom_key(b))
                    ka[0], kb[0] = model, model
                    sa, sb = output_by_key.get(tuple(ka)), output_by_key.get(tuple(kb))
                    ka[1], kb[1] = chain_a, chain_b
                    oa, ob = output_by_key.get(tuple(ka)), output_by_key.get(tuple(kb))
                    if any(x is None for x in (sa, sb, oa, ob)):
                        raise StructureEditError('Cannot propagate covalent topology: a matching residue/atom/altloc is absent in an output layer/model.')
                    def distance(x, y):
                        return float(np.linalg.norm([float(x[k])-float(y[k]) for k in ('Cartn_x','Cartn_y','Cartn_z')]))
                    reference, actual = distance(sa, sb), distance(oa, ob)
                    if abs(actual - reference) > tolerance:
                        raise StructureEditError(f'Expanded covalent link {chain_a}/{a["auth_seq_id"]}/{a["label_atom_id"]} - {chain_b}/{b["auth_seq_id"]}/{b["label_atom_id"]} would be {actual:.3f} A (source {reference:.3f} A). Review layer assignment/twist/rise; no file written.')
                    model_pairs.append((oa, ob, actual))
                first_a, first_b, length = model_pairs[0]
                new = dict(template)
                new.update(id=f'layerbond{len(rows)+1}', pdbx_dist_value=f'{length:.3f}',
                           details='Covalent pattern translated from source layer/residue identities; all model bond lengths checked')
                for prefix, atom in (('ptnr1_', first_a), ('ptnr2_', first_b)):
                    for field in ('label_asym_id','label_comp_id','label_seq_id','label_atom_id','auth_asym_id','auth_seq_id'):
                        new[prefix+field] = atom.get(field)
                    new['pdbx_'+prefix+'PDB_ins_code'] = atom.get('pdbx_PDB_ins_code')
                    new['pdbx_'+prefix+'label_alt_id'] = atom.get('label_alt_id')
                    new[prefix+'symmetry'] = '1_555'
                if pair_key(new) in seen:
                    continue
                seen.add(pair_key(new))
                rows.append(new)
                order = {'sing':1,'doub':2,'trip':3}.get(new.get('pdbx_value_order'),1)
                for oa, ob, length in model_pairs:
                    self.output_bonds.append((oa['id'], ob['id'], order))
                made += 1
        _set_cif_rows(self.output_block, '_struct_conn.', rows)
        if rows:
            types = {r['id']: r for r in _cif_rows(self.output_block,'_struct_conn_type.')}
            for r in rows:
                types.setdefault(r['conn_type_id'], {'id':r['conn_type_id']})
            _set_cif_rows(self.output_block,'_struct_conn_type.', list(types.values()))
        self.report['inferred_covalent_links'] = made
        self.report['covalent_inference'] = 'Translation of annotated source patterns using layer offsets; 0.20 A maximum bond-length change in every model.'
        self._validate_references()
        return self

    def _validate_references(self):
        block = self.output_block
        asym = {r["id"]: r.get("entity_id") for r in _cif_rows(block, "_struct_asym.")}
        entities = {r["id"] for r in _cif_rows(block, "_entity.")}
        for row in self.output_atoms:
            if row["label_asym_id"] not in asym or (_present(row.get("label_entity_id")) and row["label_entity_id"] not in entities):
                raise StructureEditError("Dangling atom_site -> struct_asym/entity reference.")
        atom_ids = {r["id"] for r in self.output_atoms}
        for row in _cif_rows(block, "_atom_site_anisotrop."):
            if row["id"] not in atom_ids:
                raise StructureEditError("Dangling anisotropic atom reference.")
        auth = {r["auth_asym_id"] for r in self.output_atoms}
        for cat in block.get_mmcif_category_names():
            if self.report['metadata'].get(cat, {}).get('action') == 'preserved/unmapped':
                continue
            for row in _cif_rows(block, cat):
                for key, val in row.items():
                    if "asym_id" in key and _present(val):
                        known = auth if "auth" in key else asym
                        if val not in known:
                            raise StructureEditError(f"Dangling chain reference {cat}{key}={val}")
        for row in _cif_rows(block, "_struct_conn."):
            if not all(self.output_index.exists(row, p) for p in ("ptnr1_", "ptnr2_")):
                raise StructureEditError("Dangling struct_conn atom endpoint.")
        sulfur_partners = defaultdict(set)
        for row in _cif_rows(block, '_struct_conn.'):
            if row.get('conn_type_id') != 'disulf':
                continue
            partners = [self.output_index.matches(row, p) for p in ('ptnr1_','ptnr2_')]
            for a in partners[0]:
                for b in partners[1]:
                    if a['pdbx_PDB_model_num'] != b['pdbx_PDB_model_num']:
                        continue
                    aa, ab = _null_id(a.get('label_alt_id')), _null_id(b.get('label_alt_id'))
                    if aa and ab and aa != ab:
                        continue
                    # Alternate positions of the same partner count once.
                    sulfur_partners[_atom_key(a)].add(_atom_key(b)[:-1])
                    sulfur_partners[_atom_key(b)].add(_atom_key(a)[:-1])
        if any(len(partners)>1 for partners in sulfur_partners.values()):
            raise StructureEditError('A disulfide sulfur would have multiple residue partners. Review the alternating-layer pattern or provide explicit covalent topology; no file written.')

    def _pdb_text(self):
        gemmi, _ = _structure_dependencies()
        ranges_by_sheet = defaultdict(list)
        for row in _cif_rows(self.output_block, '_struct_sheet_range.'):
            ranges_by_sheet[row['sheet_id']].append(row['id'])
        for row in _cif_rows(self.output_block, '_struct_sheet_order.'):
            ids = ranges_by_sheet[row['sheet_id']]
            if row['range_id_1'] not in ids or row['range_id_2'] not in ids or ids.index(row['range_id_2']) != ids.index(row['range_id_1']) + 1:
                raise StructureEditError('This sheet graph is not a sequential PDB strand list; export mmCIF to preserve its topology.')
        max_chain = 1 if self.pdb_mode == "strict" else 2
        if any(len(r["auth_asym_id"]) > max_chain for r in self.output_atoms):
            raise StructureEditError(f"This PDB mode supports at most {max_chain} character(s) per chain ID; choose extended mode or mmCIF.")
        for row in self.output_atoms:
            for field, width in (("label_atom_id", 4), ("label_comp_id", 3)):
                if len(str(row[field])) > width:
                    raise StructureEditError(f"{field} overflows PDB; export mmCIF.")
            for field in ("label_alt_id", "pdbx_PDB_ins_code"):
                if len(_null_id(row.get(field))) > 1:
                    raise StructureEditError(f"{field} overflows PDB; export mmCIF.")
            try:
                seq = int(row["auth_seq_id"])
            except (ValueError, TypeError) as exc:
                raise StructureEditError("PDB requires integer author residue numbers; export mmCIF.") from exc
            if not -999 <= seq <= 9999:
                raise StructureEditError("PDB residue number overflow; export mmCIF.")
            for field in ("Cartn_x", "Cartn_y", "Cartn_z"):
                if len(f"{float(row[field]):8.3f}") != 8:
                    raise StructureEditError("PDB coordinate field overflow; export mmCIF.")
            for field in ("occupancy", "B_iso_or_equiv"):
                if not _present(row.get(field)) or len(f"{float(row[field]):6.2f}") != 6:
                    raise StructureEditError("PDB requires representable occupancy and B factor; export mmCIF.")
        st = gemmi.make_structure_from_block(self.output_block)
        st.setup_entities()
        st.assign_serial_numbers(numbered_ter=True)
        serial_by_key = {}
        for model in st:
            for chain in model:
                for residue in chain:
                    for atom in residue:
                        key = (str(model.num), chain.name, str(residue.seqid.num), residue.seqid.icode.strip(),
                               residue.name, atom.name, atom.altloc.replace("\x00", ""))
                        serial_by_key[key] = atom.serial
            if model.count_atom_sites() + len(model) > 99999:
                raise StructureEditError("PDB atom/TER serial capacity exceeded; export mmCIF.")
        st.clear_conect()
        by_id = {r["id"]: r for r in self.output_atoms}
        bond_serials = {}
        for a, b, order in self.output_bonds:
            pair = tuple(sorted((serial_by_key[_atom_key(by_id[a])], serial_by_key[_atom_key(by_id[b])])))
            bond_serials[pair] = max(order, bond_serials.get(pair, 0))
        for (a, b), order in bond_serials.items():
            st.add_conect(a, b, order)
        options = gemmi.PdbWriteOptions()
        options.preserve_serial = True
        options.numbered_ter = True
        options.conect_records = True
        options.use_linkr = False
        lines = st.make_pdb_string(options).splitlines()
        lines = [l for l in lines if l[:6].strip() not in ("MASTER", "END")]
        lines = self._pdb_supplement(lines)
        lines = [line[:62] + '    ' + line[66:] if line.startswith('HEADER') else line for line in lines]
        at = next((i for i, line in enumerate(lines) if line[:6].strip() in ('DBREF','SEQRES','MODRES','HET','HELIX','SHEET','CRYST1','ATOM','HETATM','MODEL')), len(lines))
        lines.insert(at, 'REMARK 999 DERIVED COORDINATES FROM AMYLOIDMODIFIER; SEE .AUDIT.JSON REPORT')
        for line in lines:
            if len(line) > 80 and line not in self.unrecognized_pdb_records:
                raise StructureEditError("PDB writer produced an overlong record; export mmCIF.")
            if line[:6].strip() in ("ATOM", "HETATM", "ANISOU", "TER"):
                if not line[6:11].strip().isdigit():
                    raise StructureEditError("PDB serial overflow; hybrid-36 is not enabled. Export mmCIF.")
        counts = Counter(line[:6].strip() for line in lines)
        values = (counts["REMARK"], 0, counts["HET"], counts["HELIX"], counts["SHEET"], 0,
                  counts["SITE"], sum(counts[r + str(i)] for r in ("ORIGX", "SCALE", "MTRIX") for i in (1, 2, 3)),
                  counts["ATOM"] + counts["HETATM"], counts["TER"], counts["CONECT"], counts["SEQRES"])
        if max(values) > 99999:
            raise StructureEditError("PDB MASTER count overflow; export mmCIF.")
        lines += ["MASTER    " + "".join(f"{v:5d}" for v in values), "END"]
        if max_chain == 2 and any(len(r["auth_asym_id"]) == 2 for r in self.output_atoms):
            self.report["warnings"].append("Uses the ChimeraX/Gemmi two-character PDB chain extension; not strict wwPDB v3.3.")
        return "\n".join(l if l in self.unrecognized_pdb_records else l.ljust(80) for l in lines) + "\n"

    def _pdb_supplement(self, generated):
        rebuilt = set("HEADER TITLE KEYWDS EXPDTA NUMMDL DBREF DBREF1 DBREF2 SEQRES MODRES HET HELIX SHEET SSBOND LINK CISPEP CRYST1 ORIGX1 ORIGX2 ORIGX3 SCALE1 SCALE2 SCALE3 MTRIX1 MTRIX2 MTRIX3 MODEL ATOM HETATM ANISOU TER ENDMDL CONECT MASTER END".split())
        descriptive = {"AUTHOR", "JRNL", "HETNAM", "HETSYN", "FORMUL", "MDLTYP"}
        invalidated = {"OBSLTE", "SPLIT", "CAVEAT", "REVDAT", "SPRSDE", "REMARK"}
        extra = []
        compound, source = self._pdb_compound_source()
        extra.extend(compound + source)
        extra.extend(self._pdb_sites())
        live_compounds = {r['label_comp_id'] for r in self.output_atoms}
        for rec in dict.fromkeys(l[:6].strip() for l in self.pdb_lines):
            rows = [l for l in self.pdb_lines if l[:6].strip() == rec]
            if rec in rebuilt or rec in ("COMPND", "SOURCE", "SITE") or not rec:
                continue
            if self.metadata_policy == 'preserve' and any(line[:6].strip() == rec for line in self.unrecognized_pdb_records):
                self._audit('PDB:' + rec, 'preserved/unmapped', 'Original unrecognized records retained verbatim.', rows)
                continue
            if rec in descriptive:
                if rec in ("HETNAM", "HETSYN"):
                    rows = [r for r in rows if r[11:14].strip() in live_compounds]
                elif rec == "FORMUL":
                    # Source copy-count multipliers cannot describe a new assembly.
                    self._audit("PDB:FORMUL", "invalidated", "Chemical formula copy counts must be regenerated from chemistry, not guessed.", rows)
                    continue
                extra.extend(rows)
                self._audit("PDB:" + rec, "preserved", "Source descriptive metadata.")
            elif rec == "SEQADV":
                for group in self.groups:
                    for row in rows:
                        old = row[15:17].strip()
                        if old in group["chains"]:
                            extra.append(row[:15] + f"{group['chains'][old]:>2}" + row[17:])
                self._audit("PDB:SEQADV", "remapped")
            elif rec in invalidated or self.metadata_policy in ("audit", "preserve"):
                self._audit("PDB:" + rec, "invalidated", "Unstructured source annotation cannot be asserted for the edited model; retained in report.", rows)
            else:
                raise StructureEditError(f"No reviewed policy for PDB record {rec}; use audit policy to archive it in the report.")
        # Header sections precede sequence/connectivity annotations and coordinates.
        order = "HEADER OBSLTE TITLE SPLIT CAVEAT COMPND SOURCE KEYWDS EXPDTA NUMMDL MDLTYP AUTHOR REVDAT SPRSDE JRNL REMARK DBREF DBREF1 DBREF2 SEQADV SEQRES MODRES HET HETNAM HETSYN FORMUL HELIX SHEET SSBOND LINK CISPEP SITE CRYST1 ORIGX1 ORIGX2 ORIGX3 SCALE1 SCALE2 SCALE3 MTRIX1 MTRIX2 MTRIX3".split()
        head, body = [], []
        for line in generated + extra:
            (head if line[:6].strip() in order else body).append(line)
        head.sort(key=lambda l: order.index(l[:6].strip()))
        if self.metadata_policy == 'preserve':
            head.extend(self.unrecognized_pdb_records)
            categories = {category: _cif_rows(self.output_block, category)
                          for category, audit in self.report['metadata'].items()
                          if category.startswith('_') and audit['action'] == 'preserved/unmapped'}
            if categories:
                head.extend(_cif_annotation_records(categories))
        return head + body

    def _pdb_compound_source(self):
        import textwrap
        def parse(rec):
            joined = ' '.join(l[10:80].strip() for l in self.pdb_lines if l[:6].strip() == rec)
            groups, current = {}, None
            for token in joined.split(';'):
                if not token.strip():
                    continue
                if ':' not in token:
                    raise StructureEditError(f'Cannot safely parse {rec} key/value continuation.')
                key, value = (s.strip() for s in token.split(':', 1))
                if key == 'MOL_ID':
                    current = value
                    groups[current] = []
                elif current is None:
                    raise StructureEditError(f'{rec} lacks a MOL_ID.')
                else:
                    groups[current].append((key, value))
            return groups
        compounds = parse('COMPND')
        sources = parse('SOURCE')
        live = {}
        for mol, fields in compounds.items():
            revised = []
            keep = True
            for key, value in fields:
                if key == 'CHAIN':
                    original = [v.strip() for v in value.split(',')]
                    mapped = list(dict.fromkeys(g['chains'][v] for g in self.groups for v in original if v in g['chains']))
                    keep = bool(mapped)
                    value = ', '.join(mapped)
                revised.append((key, value))
            if keep:
                live[mol] = revised
        def render(rec, groups):
            text = ' '.join('MOL_ID: ' + mol + '; ' + ' '.join(k + ': ' + v + ';' for k, v in fields)
                            for mol, fields in groups.items())
            return [f'{rec:<6}{str(i + 1) if i else "":>4}' + part
                    for i, part in enumerate(textwrap.wrap(text, 70, break_long_words=False, break_on_hyphens=False))]
        if compounds:
            self._audit('PDB:COMPND', 'remapped', 'Molecule chain lists expanded/filtered without changing source molecule identity.')
        if sources:
            self._audit('PDB:SOURCE', 'preserved/filtered', 'Retained MOL_IDs only.')
        return render('COMPND', live), render('SOURCE', {m: f for m, f in sources.items() if m in live})

    def _pdb_sites(self):
        sites = defaultdict(list)
        for line in self.pdb_lines:
            if line[:6].strip() == 'SITE':
                for start in (18, 29, 40, 51):
                    part = line.ljust(80)[start:start + 10]
                    if part.strip():
                        sites[line[11:14].strip()].append(part)
        output, site_number = [], 0
        for group in self.groups:
            live = {(r['auth_asym_id'], r['auth_seq_id'], _null_id(r.get('pdbx_PDB_ins_code')), r['label_comp_id'])
                    for r in group['rows']}
            for old, parts in sites.items():
                kept = []
                for part in parts:
                    chain, seq, ins, res = part[3:5].strip(), part[5:9].strip(), part[9].strip(), part[:3].strip()
                    if chain in group['chains'] and (chain, seq, ins, res) in live:
                        kept.append(part[:3] + f"{group['chains'][chain]:>2}" + part[5:])
                if not kept:
                    continue
                site_number += 1
                if site_number > 999 or len(kept) > 99:
                    raise StructureEditError('SITE field overflow; export mmCIF.')
                for i in range(0, len(kept), 4):
                    output.append(f'SITE   {i // 4 + 1:3d} {site_number:3d} {len(kept):2d} ' + ' '.join(kept[i:i + 4]))
        if sites:
            self._audit('PDB:SITE', 'remapped', 'Surviving complete residue references only; site IDs and counts rebuilt.')
        return output

    def _serialize_checked(self, format, output_name):
        gemmi, np = _structure_dependencies()
        if self.output_block is None:
            raise StructureEditError("Choose an edit before writing.")
        if format not in ('cif', 'mmcif', 'pdb'):
            raise StructureEditError("Specify PDB or mmCIF output format.")
        if format != 'pdb' and self.pdb_lines:
            # Keep untranslatable PDB annotations in the in-memory provenance report.
            coordinate_records = {'ATOM','HETATM','ANISOU','TER','MODEL','ENDMDL','CONECT','MASTER','END'}
            self.report['source_pdb_metadata'] = [l for l in self.pdb_lines if l[:6].strip() not in coordinate_records]
            self.report['warnings'].append('PDB-to-mmCIF conversion: unmapped free-text/source-only PDB metadata is archived in this report, not asserted as edited mmCIF categories.')
        text = self._pdb_text() if format == "pdb" else self.output_document.as_string()
        reopened = (gemmi.read_pdb_string(text) if format == 'pdb' else
                    gemmi.make_structure_from_block(gemmi.cif.read_string(text).sole_block()))
        actual = []
        atom_properties = {}
        for model in reopened:
            for chain in model:
                for residue in chain:
                    for atom in residue:
                        key = (str(model.num), chain.name, str(residue.seqid.num), residue.seqid.icode.strip(),
                               residue.name, atom.name, atom.altloc.replace("\x00", ""))
                        actual.append((key,(atom.pos.x, atom.pos.y, atom.pos.z)))
                        atom_properties[key] = atom
        expected = {_atom_key(r): [float(r[f]) for f in ("Cartn_x", "Cartn_y", "Cartn_z")] for r in self.output_atoms}
        if len(actual) != len(expected) or {a[0] for a in actual} != set(expected):
            raise StructureEditError("Read-back atom identities/count differ from the requested edit.")
        tolerance = 0.00051 if format == "pdb" else 0.0000011
        if not np.allclose([xyz for key, xyz in actual], [expected[key] for key, xyz in actual], atol=tolerance, rtol=0):
            raise StructureEditError("Read-back coordinate precision check failed.")
        for row in self.output_atoms:
            atom = atom_properties[_atom_key(row)]
            for field, value in (('occupancy',atom.occ),('B_iso_or_equiv',atom.b_iso)):
                if _present(row.get(field)) and abs(float(row[field])-value) > (.0051 if format == 'pdb' else 1e-5) + abs(value)*1e-7:
                    raise StructureEditError('Read-back occupancy/B-factor check failed.')
            if _present(row.get('type_symbol')) and atom.element.name.upper() != str(row['type_symbol']).upper():
                raise StructureEditError('Read-back element identity differs.')
            if _present(row.get('pdbx_formal_charge')) and atom.charge != int(row['pdbx_formal_charge']):
                raise StructureEditError('Read-back formal charge differs.')
        back = reopened.make_mmcif_block(gemmi.MmcifOutputGroups(True,auth_all=True))
        def endpoint(row,prefix):
            return tuple(_null_id(row.get(prefix+k)) for k in ('auth_asym_id','auth_seq_id','label_comp_id','label_atom_id')) + (
                _null_id(row.get('pdbx_'+prefix+'PDB_ins_code',row.get(prefix+'PDB_ins_code'))),
                _null_id(row.get('pdbx_'+prefix+'label_alt_id')))
        def edge_set(block):
            return {tuple(sorted((endpoint(r,'ptnr1_'),endpoint(r,'ptnr2_'))))
                    for r in _cif_rows(block,'_struct_conn.') if r.get('conn_type_id') in ('disulf','covale','modres','metalc')}
        if not edge_set(self.output_block) <= edge_set(back):
            raise StructureEditError('Read-back lost an explicitly annotated covalent/metal connection.')
        def strand_set(block):
            return {(endpoint(r,'beg_'),endpoint(r,'end_')) for r in _cif_rows(block,'_struct_sheet_range.')}
        if strand_set(back) != strand_set(self.output_block):
            raise StructureEditError('Read-back sheet strand endpoints differ.')
        self.report.update(output=output_name, format=format, models=len(reopened), atoms=len(actual),
                           author_chains=len({r["auth_asym_id"] for r in self.output_atoms}),
                           output_sha256=hashlib.sha256(text.encode('utf-8')).hexdigest(), round_trip="passed",
                           round_trip_checks=['atom identities/counts','coordinates','occupancy/B','elements/charges','covalent/metal endpoints','sheet strand endpoints'])
        self.report["warnings"] = list(dict.fromkeys(self.report["warnings"]))
        return text

    def write(self, path, *, format=None, report_path=None, final_path=None):
        if isinstance(path, StructureBuffer):
            format = format or path.format
            text = self._serialize_checked(format, 'memory:' + path.name)
            # A failed check leaves both the previous coordinates and audit intact.
            path.text, path.report = text, deepcopy(self.report)
            path.format = 'pdb' if format == 'pdb' else 'cif'
            return self.report
        path = os.path.abspath(os.fspath(path))
        if format is None:
            suffix = os.path.splitext(path)[1].lower()
            format = 'cif' if suffix in ('.cif', '.mmcif') else 'pdb' if suffix in ('.pdb', '.ent') else None
        report_path = os.path.abspath(os.fspath(report_path or path + '.audit.json'))
        if report_path in (path, self.path):
            raise StructureEditError("Audit report must not overwrite a coordinate file.")
        output_name = os.path.abspath(os.fspath(final_path)) if final_path else path
        text = self._serialize_checked(format, output_name)
        suffix = '.pdb' if format == 'pdb' else '.cif'
        fd, temporary = tempfile.mkstemp(prefix='.amyloid-', suffix=suffix, dir=os.path.dirname(path))
        try:
            with os.fdopen(fd, 'w', encoding='utf-8', newline='\n') as handle:
                handle.write(text)
                handle.flush()
                os.fsync(handle.fileno())
            report_fd, report_tmp = tempfile.mkstemp(prefix='.amyloid-report-', dir=os.path.dirname(report_path))
            try:
                with os.fdopen(report_fd, 'w', encoding='utf-8') as handle:
                    json.dump(self.report, handle, indent=2, ensure_ascii=False)
                os.replace(report_tmp, report_path)
            finally:
                if os.path.exists(report_tmp):
                    os.unlink(report_tmp)
            os.replace(temporary, path)
        finally:
            if os.path.exists(temporary):
                os.unlink(temporary)
        return self.report


def _fit_detection_core(p, q):
    # Fit a stable 60% core so flexible termini do not dominate the transform.
    import numpy as np
    n = len(p)
    keep = max(6, int(math.ceil(n * .60)))
    seeds = [np.arange(n)]
    seeds.extend(np.arange(i, min(i + 12, n)) for i in range(0, n - 5, 6))
    best = None
    for mask in seeds:
        for _ in range(10):
            pc, qc = p[mask].mean(0), q[mask].mean(0)
            if np.linalg.matrix_rank(p[mask] - pc) < 2:
                break
            u, _, vt = np.linalg.svd((p[mask] - pc).T @ (q[mask] - qc))
            fix = np.eye(3)
            fix[2, 2] = np.linalg.det(vt.T @ u.T)
            rotation = vt.T @ fix @ u.T
            error = np.linalg.norm((p - pc) @ rotation.T + qc - q, axis=1)
            new = np.argsort(error, kind='stable')[:keep]
            score = float(np.sqrt(np.mean(error[new] ** 2)))
            if best is None or score < best['rmsd']:
                best = dict(rmsd=score, rotation=rotation, mask=new, error=error)
            if np.array_equal(np.sort(mask), np.sort(new)):
                break
            mask = new
        if best is not None and best['rmsd'] < .15:
            break
    return best


def _internal_rung_evidence(points, axis, spacing=4.8, max_rungs=6):
    import numpy as np
    n = len(points)
    counts = []
    for rung in range(1, max_rungs):
        pairs = 0
        participants = set()
        # Blocks avoid allocating an entire long-chain N x N x 3 tensor.
        for start in range(0, n, 128):
            delta = points[start:start + 128, None, :] - points[None, :, :]
            axial = delta @ axis
            transverse2 = np.maximum(0., np.sum(delta * delta, axis=2) - axial * axial)
            ii = np.arange(start, min(start + 128, n))[:, None]
            jj = np.arange(n)[None, :]
            hit = ((np.abs(axial - rung * spacing) <= .9) &
                   (transverse2 <= 2.5 ** 2) & (np.abs(ii - jj) >= 4))
            rows, cols = np.where(hit)
            pairs += len(rows)
            participants.update((rows + start).tolist())
            participants.update(cols.tolist())
        counts.append(dict(separation=rung, pairs=pairs, residues=len(participants)))
    supported = [r['separation'] + 1 for r in counts
                 if r['pairs'] >= 6 and r['residues'] >= max(12, math.ceil(n * .12))]
    return max([1] + supported), counts


def _antiparallel_contacts(p, q, axis=None):
    return _backbone_stack_contacts(p, q, axis, antiparallel=True)


def _backbone_stack_contacts(p, q, axis=None, antiparallel=True):
    import numpy as np
    direction = -1 if antiparallel else 1
    required = max(6, math.ceil(.25 * min(len(p), len(q))))
    if min(len(p), len(q)) < required + 2:
        return None
    # Reject distant chains before allocating contact blocks.
    separation = np.maximum(0., np.maximum(p.min(0) - q.max(0), q.min(0) - p.max(0)))
    if np.linalg.norm(separation) > 7.2:
        return None

    def tangents(points):
        bonds = np.linalg.norm(np.diff(points, axis=0), axis=1)
        vectors = points[2:] - points[:-2]
        lengths = np.linalg.norm(vectors, axis=1)
        valid = ((bonds[:-1] >= 2.8) & (bonds[:-1] <= 4.3) &
                 (bonds[1:] >= 2.8) & (bonds[1:] <= 4.3) & (lengths >= 5.4))
        return vectors / np.maximum(lengths[:, None], 1e-12), valid

    u, valid_p = tangents(p)
    v, valid_q = tangents(q)
    contacts = []
    for start in range(0, len(u), 128):
        delta = q[None, 1:-1] - p[1:-1][start:start + 128, None]
        distance2 = np.sum(delta * delta, axis=2)
        hit = ((direction * (u[start:start + 128] @ v.T) > .6) &
               valid_p[start:start + 128, None] & valid_q[None, :] &
               (distance2 >= 3.5 ** 2) & (distance2 <= 7.2 ** 2))
        cost = distance2
        if axis is not None:
            axial = delta @ axis
            cost = np.maximum(0., distance2 - axial * axial)
            hit &= (np.abs(axial) >= 3.5) & (np.abs(axial) <= 6.5) & (cost <= 3.0 ** 2)
        ii, jj = np.where(hit)
        contacts.extend((float(cost[i, j]), start + int(i), int(j), delta[i, j])
                        for i, j in zip(ii, jj))
    # One residue cannot inflate the support by contacting multiple partners.
    pairs, used_p, used_q = [], set(), set()
    for cost, i, j, vector in sorted(contacts, key=lambda c: c[:3]):
        if i not in used_p and j not in used_q:
            pairs.append((i, j, vector))
            used_p.add(i)
            used_q.add(j)
    if len(pairs) < required:
        return None
    registers = {(i, j) for i, j, _ in pairs}
    if not any(all((i + step, j + direction * step) in registers for step in range(4))
               for i, j in registers):
        return None
    vectors = np.array([d for _, _, d in pairs])
    if axis is None:
        strand = np.array([u[i] + direction * v[j] for i, j, _ in pairs])
        strand /= np.linalg.norm(strand, axis=1)[:, None]
        transverse = vectors - np.sum(vectors * strand, axis=1)[:, None] * strand
        direction = transverse.mean(0)
        norm = float(np.linalg.norm(direction))
        if not 3.5 <= norm <= 6.5:
            return None
        axis = direction / norm
        # Recheck contacts against the inferred direction and common spacing.
        return _backbone_stack_contacts(p, q, axis, antiparallel)
    axial = vectors @ axis
    sign = 1 if np.median(axial) > 0 else -1
    if np.count_nonzero(sign * axial > 0) < .85 * len(pairs):
        return None
    rise = float(np.mean(sign * axial))
    if not 4.0 <= rise <= 5.7:
        return None
    # Contacts within a long, overlapping multi-rung chain must not merge stacks.
    center_shift = q.mean(0) - p.mean(0)
    center_rise = float(center_shift @ axis) * sign
    # Transverse alignment is checked on the contact residues above. A whole-chain lateral centroid shifts when termini are missing or flexible.
    if not 2.5 <= center_rise <= 7.5:
        return None
    return dict(vector=axis * sign * rise, distance=rise, rise=rise,
                support=len(pairs), contact_pairs=[(i + 1, j + 1) for i, j, _ in pairs],
                rmsd=None, orientation='antiparallel' if antiparallel else 'parallel',
                rungs=1, gap=1, spacing=rise)


def detect_layers(chains, axis_hint=None):
    # Detect whole chain units; a chain may span multiple physical layers.
    import numpy as np
    ids = sorted(chains)
    if not ids:
        raise StructureEditError('No protein C-alpha atoms available for detection.')
    maps, points = {}, {}
    for cid in ids:
        maps[cid] = {}
        for residue in chains[cid]:
            key = tuple(residue.get('key', (residue['id'],)))
            if key in maps[cid]:
                raise StructureEditError('Duplicate C-alpha residue identity in chain ' + cid)
            coord = np.asarray(residue['coord'], dtype=float)
            if coord.shape != (3,) or not np.all(np.isfinite(coord)):
                raise StructureEditError('Invalid C-alpha coordinates in chain ' + cid)
            maps[cid][key] = coord
        points[cid] = np.array(list(maps[cid].values()))
        if not len(points[cid]):
            raise StructureEditError('Empty C-alpha chain ' + cid)

    candidates = []
    contact_candidates = []
    for i, a in enumerate(ids):
        for b in ids[i + 1:]:
            anti = _antiparallel_contacts(points[a], points[b])
            if anti is not None:
                contact_candidates.append(dict(anti, a=a, b=b))
            common = [k for k in maps[a] if k in maps[b]]
            # Use local backbone contacts when neighboring proteins have different sequences.
            sequence_a = tuple(k[-1] for k in maps[a] if len(k) >= 3)
            sequence_b = tuple(k[-1] for k in maps[b] if len(k) >= 3)
            if sequence_a != sequence_b:
                contact = _backbone_stack_contacts(points[a], points[b], antiparallel=False)
                if contact is not None:
                    contact_candidates.append(dict(contact, a=a, b=b))
            if len(common) < max(6, math.ceil(.5 * min(len(maps[a]), len(maps[b])))):
                continue
            p = np.array([maps[a][k] for k in common])
            q = np.array([maps[b][k] for k in common])
            lengths = np.linalg.norm(q - p, axis=1)
            if np.median(lengths) > 36 or np.median(lengths) < 3.5:
                continue
            fit = _fit_detection_core(p, q)
            if fit is None or fit['rmsd'] > 1.0:
                continue
            angle = math.degrees(math.acos(float(np.clip((np.trace(fit['rotation']) - 1) / 2, -1, 1))))
            if angle > 25:
                continue
            mask = fit['mask']
            vector = (q - p)[mask].mean(0)
            distance = float(np.linalg.norm(vector))
            if not 3.8 <= distance <= 34.2:
                continue
            # Long repeats need internal rung evidence. Ordinary missing-chain gaps may be at most three rungs; they never increase rungs/unit.
            if distance > 17.1:
                ka, _ = _internal_rung_evidence(points[a], vector / distance)
                kb, _ = _internal_rung_evidence(points[b], vector / distance)
                if distance > 5.7 * min(ka, kb):
                    continue
            candidates.append(dict(a=a, b=b, vector=vector, distance=distance,
                                   rmsd=fit['rmsd'], support=len(mask),
                                   orientation='parallel',
                                   center_a=p[mask].mean(0), center_b=q[mask].mean(0)))

    warnings = []
    hint = None
    if axis_hint is not None:
        hint = np.asarray(axis_hint, dtype=float)
        if hint.shape != (3,) or not np.all(np.isfinite(hint)) or np.linalg.norm(hint) < 1e-8:
            raise StructureEditError('An axis hint must be a finite nonzero three-vector.')
        hint = hint / np.linalg.norm(hint)
    axis = None
    edges = []
    rung_evidence = {}
    multiplicity = {c: 1 for c in ids}
    axis_candidates = candidates or contact_candidates
    if axis_candidates:
        # Estimate direction from nearest repeats, independent of sign and lateral separation.
        seeds = {min((i for i, e in enumerate(axis_candidates) if c in (e['a'], e['b'])),
                     key=lambda i: (axis_candidates[i]['distance'], axis_candidates[i]['rmsd'] or 0.))
                 for c in ids if any(c in (e['a'], e['b']) for e in axis_candidates)}
        directions = np.array([axis_candidates[i]['vector'] / axis_candidates[i]['distance'] for i in sorted(seeds)])
        # Choose the strongest angular cluster before estimating its direction.
        agreement = np.abs(directions @ directions.T) >= math.cos(math.radians(25))
        selected = agreement[np.argmax(agreement.sum(1))]
        _, _, vt = np.linalg.svd(directions[selected], full_matrices=False)
        axis = vt[0]
        if axis[np.argmax(np.abs(axis))] < 0:
            axis = -axis
        for c in ids:
            multiplicity[c], rung_evidence[c] = _internal_rung_evidence(points[c], axis)
        for e in candidates:
            dz = float(e['vector'] @ axis)
            lateral = float(np.linalg.norm(e['vector'] - dz * axis))
            k = min(multiplicity[e['a']], multiplicity[e['b']])
            # A long linker may create extra contacts in one chain; multiplicity must also agree with the observed repeat, on both sides.
            repeat_rungs = max(1, round(abs(dz) / 4.8))
            k = min(k, repeat_rungs)
            gap = max(1, round(abs(dz) / (4.8 * k)))
            spacing = abs(dz) / (gap * k)
            if gap > 3 or not 4.0 <= spacing <= 5.7 or lateral > max(2.0, .30 * abs(dz)):
                continue
            if dz < 0:
                e = dict(e, a=e['b'], b=e['a'], center_a=e['center_b'], center_b=e['center_a'])
            edges.append(dict(e, rise=abs(dz), rungs=k, gap=gap, spacing=spacing))

        for candidate in contact_candidates:
            a, b = candidate['a'], candidate['b']
            if any({a, b} == {e['a'], e['b']} and e['gap'] == 1 for e in edges):
                continue  # Keep an already supported rigid repeat unchanged.
            if multiplicity[a] != 1 or multiplicity[b] != 1:
                continue
            e = _backbone_stack_contacts(points[a], points[b], axis,
                                         antiparallel=candidate['orientation'] == 'antiparallel')
            if e is not None:
                if float(e['vector'] @ axis) < 0:
                    a, b = b, a
                    e['contact_pairs'] = [(j, i) for i, j in e['contact_pairs']]
                edges.append(dict(e, a=a, b=b))

    # Use the same disjoint paths for grouping and ordering chains.
    up, down, chosen = {}, {}, []
    for e in sorted(edges, key=lambda e: (e['gap'], e['rise'], e['rmsd'] or 0., e['a'], e['b'])):
        a, b = e['a'], e['b']
        if a in up or b in down:
            continue
        cursor = b
        while cursor in up and cursor != a:
            cursor = up[cursor]['b']
        if cursor == a:
            continue
        up[a] = e
        down[b] = e
        chosen.append(e)

    if not chosen and hint is not None:
        axis = hint
        for c in ids:
            multiplicity[c], rung_evidence[c] = _internal_rung_evidence(points[c], axis)
    tracks = []
    for start in ids:
        if start in down:
            continue
        path, unit_positions = [start], [0]
        current = start
        while current in up:
            e = up[current]
            current = e['b']
            path.append(current)
            unit_positions.append(unit_positions[-1] + e['gap'])
        path_edges = [up[c] for c in path[:-1]]
        k = int(round(float(np.median([e['rungs'] for e in path_edges])))) if path_edges else (multiplicity[start] if hint is not None and not chosen else 1)
        rise = float(np.median([e['rise'] / e['gap'] for e in path_edges])) if path_edges else None
        if path_edges and any(e['rungs'] != k for e in path_edges):
            warnings.append('Inconsistent internal layer multiplicity in chain stack ' + ', '.join(path))
        orientations = {e['orientation'] for e in path_edges if e['gap'] == 1}
        orientation = ('mixed' if len(orientations) > 1 else next(iter(orientations))) if orientations else 'undetermined'
        tracks.append(dict(chains=path, units=len(path), layers_per_unit=k, orientation=orientation,
                           layers=len(path) * k, unit_positions=unit_positions,
                           axial_layer_span=(unit_positions[-1] + 1) * k,
                           unit_rise=rise, layer_spacing=rise / k if rise else None,
                           complete=unit_positions == list(range(len(path)))))

    tracks.sort(key=lambda t: (-t['units'], tuple(t['chains'])))
    if not chosen and hint is None:
        axis = None
        warnings.append('No inter-unit repeat detected: axis and layers per unit are undetermined; singleton chains are shown provisionally as one layer each.')
    elif not chosen:
        warnings.append('Axis supplied from the previously detected parent stack; internal contacts were rechecked. This single unit does not independently determine a repeat axis.')
    elif any(t['units'] == 1 for t in tracks):
        warnings.append('Some chains have no supported stacking neighbor; their layer multiplicity is undetermined.')
    if any(not t['complete'] for t in tracks):
        warnings.append('Missing chain units detected; occupied layers and axial layer span differ.')
    uniform = len({t['layers_per_unit'] for t in tracks}) == 1
    if not uniform:
        warnings.append('Protofilaments have different layers per chain unit.')
    # Keep every input chain exactly once, even in a partial/ambiguous assembly.
    assert sorted(c for t in tracks for c in t['chains']) == ids
    orientations = {t['orientation'] for t in tracks} - {'undetermined'}
    orientation = ('mixed' if len(orientations) > 1 else next(iter(orientations))) if orientations else 'undetermined'
    return dict(version='3.1', protofilaments=tracks, orientation=orientation,
                antiparallel=any(e['orientation'] == 'antiparallel' for e in chosen),
                sandwiches=[t['chains'] for t in tracks],
                axis=axis.tolist() if axis is not None else None,
                axis_source=('neighbor repeats' if candidates else
                             ('antiparallel backbone contacts' if orientation == 'antiparallel' else 'backbone contacts')) if chosen else ('supplied parent axis' if hint is not None else 'undetermined'),
                detected_units=max(t['units'] for t in tracks),
                detected_layers=max(t['layers'] for t in tracks),
                layers_per_unit=tracks[0]['layers_per_unit'] if uniform else None,
                complete=uniform and all(t['complete'] for t in tracks) and
                         len({t['units'] for t in tracks}) == 1,
                warnings=warnings, rung_evidence=rung_evidence,
                accepted_edges=[{k: e[k] for k in ('a', 'b', 'rise', 'rungs', 'gap', 'rmsd', 'support', 'orientation', 'contact_pairs') if k in e} for e in chosen])


def _expansion_repeat_pattern(chains, sandwiches):
    import numpy as np
    tags, contact_vectors, pure_anti = [], [], False
    for stack in sandwiches:
        points = [np.array([r['coord'] for r in chains[c]]) for c in stack]
        polarity, flips, stack_tags = 1, [], []
        for i, cid in enumerate(stack):
            if i:
                contact = _antiparallel_contacts(points[i - 1], points[i])
                flips.append(contact is not None)
                if contact is not None:
                    polarity *= -1
                    contact_vectors.append(contact['vector'])
            identity = frozenset(tuple(r.get('key', (r['id'],))) for r in chains[cid])
            stack_tags.append((identity, polarity))
        pure_anti |= bool(flips) and all(flips)
        tags.append(stack_tags)
    depth = len(sandwiches[0])
    def compatible(period):
        return (not (pure_anti and period % 2) and
                all(all(t[i] == t[i + period] for i in range(len(t) - period)) for t in tags))
    period = next((p for p in range(1, depth + 1) if compatible(p)), depth + 1)
    axis = None
    if contact_vectors:
        _, _, vt = np.linalg.svd(np.array(contact_vectors), full_matrices=False)
        axis = vt[0]
        if np.dot(axis, np.mean(contact_vectors, axis=0)) < 0:
            axis = -axis
    return dict(period=period, compatible=compatible, axis=axis,
                antiparallel=bool(contact_vectors))


def expand_structure_layers(editor, sandwiches, layers_to_add, use_auto, manual_twist,
                            manual_rise, water_z_limit=4.0, use_alt=False,
                            use_computed_axis=True, label_generator=None,
                            repeat_units=None, random_seed=None):
    # Keep every input atom fixed; use one shared transform only for added copies.
    # repeat_units selects the cycle length; zero or negative values sample compatible templates.
    # Twist and rise are per chain unit, even for chains spanning several layers.
    _, np = _structure_dependencies()
    if not isinstance(layers_to_add, int) or layers_to_add <= 0:
        raise StructureEditError("Number of added layers must be a positive integer.")
    if not sandwiches or len({len(s) for s in sandwiches}) != 1 or not sandwiches[0]:
        raise StructureEditError("Expansion requires equally populated, ordered protofilaments.")
    core = [c for s in sandwiches for c in s]
    if len(set(core)) != len(core):
        raise StructureEditError("A chain occurs in more than one protofilament/layer.")
    ca = editor.ca_chains()
    if set(core) - set(ca):
        raise StructureEditError("Some selected chains have no polymer C-alpha atoms.")
    depth = len(sandwiches[0])
    pattern = _expansion_repeat_pattern(ca, sandwiches)
    if repeat_units is None:
        repeat_units = 2 if use_alt else pattern['period']
    if isinstance(repeat_units, bool) or not isinstance(repeat_units, int):
        raise StructureEditError('Alternating period must be a whole number; zero or any negative number selects random.')
    randomize = repeat_units <= 0
    if randomize:
        repeat_units = -1  # Canonical audit value for every random-mode alias.
    step = pattern['period'] if randomize else repeat_units
    if step > depth:
        raise StructureEditError(f'A repeat of {step} chain units requires at least {step} source units; only {depth} are available.')
    if not pattern['compatible'](step):
        raise StructureEditError('The repeat period would exchange different protein identities or opposite backbone orientations. '
                                 f'Use a compatible cycle (suggested: {pattern["period"]} chain units).')
    rng = None
    if randomize:
        if random_seed is None:
            random_seed = int(np.random.SeedSequence().entropy)
        if isinstance(random_seed, bool) or not isinstance(random_seed, int) or random_seed < 0:
            raise StructureEditError('Random seed must be a nonnegative integer.')
        rng = np.random.default_rng(random_seed)
    centers = np.array([np.mean([r["coord"] for s in sandwiches for r in ca[s[i]]], axis=0) for i in range(depth)])
    center = centers.mean(axis=0)
    axis = np.array([0., 0., 1.])
    if use_computed_axis and depth >= 2:
        # Same-phase displacements avoid the transverse stagger between A/B.
        repeated = centers[step:] - centers[:-step]
        axis_points = repeated if step > 1 and len(repeated) else centers - center
        _, singular, vt = np.linalg.svd(axis_points)
        if singular[0] < 1e-8:
            raise StructureEditError("Layer centers coincide; helical axis is indeterminate.")
        axis = vt[0]
        if not len(repeated) and pattern['axis'] is not None:
            axis = pattern['axis']
        if np.dot(axis, centers[-1] - centers[0]) < 0:
            axis = -axis
    rmsd = None
    if use_auto:
        if depth <= step:
            raise StructureEditError(f'Automatic fitting for a {step}-unit repeat needs at least {step + 1} source chain units. '
                                     'Load another repeat or enter manual twist/rise.')
        source, target, pair_ranges, repeat_directions = [], [], [], []
        for s in sandwiches:
            for i in range(depth - step):
                a = {r["key"]: r["coord"] for r in ca[s[i]]}
                b = {r["key"]: r["coord"] for r in ca[s[i + step]]}
                if set(a) != set(b):
                    raise StructureEditError("C-alpha residue identities differ between equivalent repeat positions; resolve alignment before automatic fitting.")
                pair_p, pair_q = np.array(list(a.values())), np.array([b[k] for k in a])
                pair_ranges.append((len(source), len(source) + len(a)))
                source.extend(pair_p)
                target.extend(pair_q)
                if use_computed_axis and (step > 1 or randomize) and len(a) >= 6:
                    fit = _fit_detection_core(pair_p, pair_q)
                    if fit is None or fit['rmsd'] > 1.5:
                        raise StructureEditError('No stable core between equivalent repeat positions.')
                    repeat_directions.append((pair_q - pair_p)[fit['mask']].mean(0))
        p, q = np.asarray(source), np.asarray(target)
        full_p, full_q = p, q
        if repeat_directions:
            axis = np.mean(repeat_directions, axis=0)
            norm = np.linalg.norm(axis)
            if norm < 1e-8:
                raise StructureEditError('Same-phase repeat directions cancel; check the ordered stacks.')
            axis /= norm
        total_fit_atoms = len(p)
        # Fit the stable core; copy complete chains, including their flexible regions.
        if use_computed_axis and len(p) >= 6:
            core_fit = _fit_detection_core(p, q)
            if core_fit is None or core_fit['rmsd'] > 1.0:
                raise StructureEditError('No consistent stable core for automatic expansion.')
            p, q = p[core_fit['mask']], q[core_fit['mask']]
        pc, qc = p.mean(axis=0), q.mean(axis=0)
        if len(p) < 3 or np.linalg.matrix_rank(p - pc) < 2:
            raise StructureEditError("At least three noncollinear matched C-alpha atoms are needed.")
        if use_computed_axis and (step > 1 or randomize):
            # Constrain periodic fits to the repeat direction to avoid unstable NMR screw axes.
            pp, qq = p - pc, q - qc
            pp -= (pp @ axis)[:, None] * axis
            qq -= (qq @ axis)[:, None] * axis
            angle = math.atan2(float(np.sum(np.cross(pp, qq) @ axis)), float(np.sum(pp * qq)))
            k = np.array([[0, -axis[2], axis[1]], [axis[2], 0, -axis[0]], [-axis[1], axis[0], 0]])
            rotation = np.eye(3) + math.sin(angle) * k + (1 - math.cos(angle)) * (k @ k)
        elif use_computed_axis:
            u, _, vt = np.linalg.svd((p - pc).T @ (q - qc))
            fix = np.eye(3)
            fix[2, 2] = np.linalg.det(vt.T @ u.T)
            rotation = vt.T @ fix @ u.T
            vals, vecs = np.linalg.eig(rotation)
            if not np.allclose(rotation, np.eye(3), atol=1e-7):
                axis = np.real(vecs[:, np.argmin(abs(vals - 1))])
                axis /= np.linalg.norm(axis)
                if np.dot(axis, centers[-1] - centers[0]) < 0:
                    axis = -axis
        else:
            pp, qq = p - pc, q - qc
            angle = math.atan2(np.sum(pp[:, 0] * qq[:, 1] - pp[:, 1] * qq[:, 0]),
                               np.sum(pp[:, 0] * qq[:, 0] + pp[:, 1] * qq[:, 1]))
            rotation = np.array([[math.cos(angle), -math.sin(angle), 0], [math.sin(angle), math.cos(angle), 0], [0, 0, 1]])
        translation = qc - rotation @ pc
        skew = np.array([rotation[2, 1] - rotation[1, 2], rotation[0, 2] - rotation[2, 0], rotation[1, 0] - rotation[0, 1]]) / 2
        twist = math.degrees(math.atan2(float(np.dot(skew, axis)), float(np.clip((np.trace(rotation) - 1) / 2, -1, 1)))) / step
        rise = float(np.dot(translation, axis)) / step
        rmsd = float(np.sqrt(np.mean(np.sum((p @ rotation.T + translation - q) ** 2, axis=1))))
        if step > 1 or randomize:
            # A pooled trimmed fit must not discard an entire protein/phase.
            errors = np.linalg.norm(full_p @ rotation.T + translation - full_q, axis=1)
            phase_errors = []
            for first, last in pair_ranges:
                keep = max(3, math.ceil(.6 * (last - first)))
                core_error = float(np.sqrt(np.mean(np.sort(errors[first:last])[:keep] ** 2)))
                phase_errors.append(core_error)
                if core_error > 1.5:
                    raise StructureEditError('The source units do not share a consistent repeat transform in every phase; choose a different period or manual parameters.')
        editor.report['automatic_fit'] = {'method': 'stable-core proper rigid fit' if use_computed_axis else 'Z-axis fit',
                                        'matched_ca_atoms': total_fit_atoms, 'fitted_ca_atoms': len(p),
                                        'core_rmsd': rmsd, 'period_units': step,
                                        'axis_constrained_to_repeats': bool(use_computed_axis and (step > 1 or randomize))}
        if step > 1 or randomize:
            editor.report['automatic_fit']['pair_core_rmsds'] = phase_errors
    else:
        twist, rise = float(manual_twist), float(manual_rise)
        if not math.isfinite(twist) or not math.isfinite(rise) or abs(rise) < 1e-8:
            raise StructureEditError("Twist/rise must be finite, with nonzero rise.")
        angle = math.radians(twist * step)
        k = np.array([[0, -axis[2], axis[1]], [axis[2], 0, -axis[0]], [-axis[1], axis[0], 0]])
        rotation = np.eye(3) + math.sin(angle) * k + (1 - math.cos(angle)) * (k @ k)
        translation = center - rotation @ center + rise * step * axis
    if abs(rise) < 1e-8:
        raise StructureEditError("Estimated helical rise is zero; refusing overlapping copies.")
    owners = editor.associated_residues(set(core), water_z_limit, axis)
    used = set(editor.chains)
    next_index = 0
    def allocate():
        nonlocal next_index
        while True:
            label = (label_generator or structure_chain_id)(next_index)
            next_index += 1
            if label not in used:
                used.add(label)
                return label
    # The identity copy retains all existing atoms; only new chains are transformed.
    copies = [{"chains": {c: c for c in editor.chains}}]
    template_choices = []
    instances = {(pf, layer): chain for pf, strand in enumerate(sandwiches) for layer, chain in enumerate(strand)}
    top = layers_to_add // 2
    bottom = layers_to_add - top
    alignment_cache = {}
    compatibility_cache = {}
    source_at = {i: i for i in range(depth)}
    # Sample outwards, checking compatibility with each placed same-phase neighbor.
    lower = list(range(-1, -bottom - 1, -1)) if randomize else list(range(-bottom, 0))
    for dest in lower + list(range(depth, depth + top)):
        anchor = dest % step if dest < 0 else depth - step + (dest - depth) % step
        alignment = np.eye(4)
        if randomize:
            neighbor = source_at[dest + step if dest < 0 else dest - step]
            eligible = []
            for candidate in range(dest % step, depth, step):
                key = tuple(sorted((candidate, neighbor)))
                if key not in compatibility_cache:
                    compatible = True
                    for stack in sandwiches:
                        if candidate == neighbor:
                            continue
                        a = {r['key']: r['coord'] for r in ca[stack[candidate]]}
                        b = {r['key']: r['coord'] for r in ca[stack[neighbor]]}
                        fit = _fit_detection_core(np.array(list(a.values())), np.array([b[k] for k in a])) if len(a) >= 6 else None
                        if fit is None or fit['rmsd'] > 1.0:
                            compatible = False
                            break
                    compatibility_cache[key] = compatible
                if compatibility_cache[key]:
                    eligible.append(candidate)
            ref = eligible[int(rng.integers(len(eligible)))]
            # Align sampled copies to the boundary core to avoid accumulating donor-position errors.
            if ref != anchor:
                key = ref, anchor
                if key not in alignment_cache:
                    donor, reference = [], []
                    for stack in sandwiches:
                        a = {r['key']: r['coord'] for r in ca[stack[ref]]}
                        b = {r['key']: r['coord'] for r in ca[stack[anchor]]}
                        donor.extend(a.values())
                        reference.extend(b[k] for k in a)
                    p, q = np.array(donor), np.array(reference)
                    fit = _fit_detection_core(p, q) if len(p) >= 6 else None
                    if fit is None or fit['rmsd'] > 1.5:
                        raise StructureEditError('A sampled unit has no consistent stable core with its repeat position.')
                    p, q = p[fit['mask']], q[fit['mask']]
                    pc, qc = p.mean(0), q.mean(0)
                    u, _, vt = np.linalg.svd((p - pc).T @ (q - qc))
                    fix = np.eye(3)
                    fix[2, 2] = np.linalg.det(vt.T @ u.T)
                    r = vt.T @ fix @ u.T
                    transform = np.eye(4)
                    transform[:3, :3], transform[:3, 3] = r, qc - r @ pc
                    alignment_cache[key] = transform
                alignment = alignment_cache[key]
        else:
            ref = anchor
        source_at[dest] = ref
        power = (dest - anchor) // step
        homogeneous = np.eye(4)
        homogeneous[:3, :3], homogeneous[:3, 3] = rotation, translation
        transformed = np.linalg.matrix_power(homogeneous, power) @ alignment
        source_chains = {s[ref] for s in sandwiches}
        selected = [r for r in editor.atoms if _null_id(r["auth_asym_id"]) in source_chains or owners.get(_residue_key(r)) in source_chains]
        chain_order = list(dict.fromkeys(_null_id(r["auth_asym_id"]) for r in selected))
        mapping = {c: allocate() for c in chain_order}
        for pf, strand in enumerate(sandwiches):
            instances[pf, dest] = mapping[strand[ref]]
        template_choices.append(dict(destination_unit=dest, source_unit=ref, anchor_unit=anchor,
                                     eligible_templates=len(eligible) if randomize else 1,
                                     source_chains=[s[ref] for s in sandwiches],
                                     destination_chains=[instances[pf, dest] for pf in range(len(sandwiches))]))
        copies.append({"chains": mapping, "atom_ids": {r["id"] for r in selected},
                       "matrix": transformed[:3, :3], "vector": transformed[:3, 3]})
    editor.report["operation"] = "expand layers"
    editor.report["geometry"] = {"twist_degrees": twist, "rise_angstrom": rise, "axis": axis.tolist(),
                                 "fit_rmsd_angstrom": rmsd, "alternating": step > 1,
                                 "repeat_units": repeat_units, "transform_period_units": step,
                                 "random_seed": random_seed if randomize else None,
                                 "template_choices": template_choices,
                                 "cycle_matrix": rotation.tolist(), "cycle_vector": translation.tolist(),
                                 "manual_axis_point": center.tolist() if not use_auto else None}
    editor.apply_copies(copies)
    editor.propagate_layer_bonds(sandwiches, instances, step)
    editor.propagate_layer_sheets(sandwiches, instances, step)
    source_atoms = {_atom_key(row): row for row in editor.atoms}
    retained = {_atom_key(row): row for row in editor.output_atoms if _atom_key(row) in source_atoms}
    if retained.keys() != source_atoms.keys() or any(
            retained[key][field] != row[field]
            for key, row in source_atoms.items() for field in ('Cartn_x', 'Cartn_y', 'Cartn_z')):
        raise StructureEditError('Expansion changed an existing atom; the edit was rejected.')
    editor.report['source_preservation'] = {'atoms': len(source_atoms), 'coordinates': 'unchanged'}
    tilt = math.degrees(math.acos(float(np.clip(abs(axis[2]), -1, 1))))
    return twist, rise, axis, tilt


def select_live_ensemble_member(editor, live, scene_position=None):
    # ChimeraX model IDs are session IDs, not source MODEL numbers.
    if len(editor.models) <= 1 or len(live.models) != 1:
        return False
    _, np = _structure_dependencies()
    live_keys = [_atom_key(r, model=False) for r in live.atoms]
    xyz = np.array([[float(r[f]) for f in ('Cartn_x', 'Cartn_y', 'Cartn_z')] for r in live.atoms])
    if scene_position is not None:
        xyz = scene_position.inverse().transform_points(xyz)
    by_model = defaultdict(dict)
    for row in editor.atoms:
        by_model[str(row['pdbx_PDB_model_num'])][_atom_key(row, model=False)] = row
    matches = []
    for number, atoms in by_model.items():
        if any(key not in atoms for key in live_keys):
            continue
        source_xyz = [[float(atoms[key][f]) for f in ('Cartn_x', 'Cartn_y', 'Cartn_z')] for key in live_keys]
        if np.allclose(xyz, source_xyz, atol=.0011, rtol=0):
            matches.append(number)
    if len(matches) != 1:
        raise StructureEditError('Cannot uniquely identify the selected ensemble member in the original file after coordinate edits. Save that member as a separate PDB/mmCIF, open it, and load it; source conformer identity will not be guessed.')
    editor.select_model(matches[0], live.models[0])
    return True


def merge_structure_coordinates(editor, live_path):
    # Import live coordinates only after checking source atom identities.
    live = live_path if isinstance(live_path, StructureEditor) else StructureEditor(live_path)
    source_keys = {_atom_key(r): r for r in editor.atoms}
    live_keys = {_atom_key(r): r for r in live.atoms}
    if set(live.models) != set(editor.models):
        raise StructureEditError('The live selection omits or renumbers ensemble models. Select/export the complete ensemble, or explicitly load a single model.')
    if set(live_keys) - set(source_keys):
        raise StructureEditError('Live atom identities changed. Reload that model explicitly before editing; metadata cannot be mapped by atom order.')
    fields = ('Cartn_x', 'Cartn_y', 'Cartn_z', 'occupancy', 'B_iso_or_equiv')
    changed = set(live_keys) != set(source_keys)
    for key, row in live_keys.items():
        source = source_keys[key]
        for field in fields:
            if _present(row.get(field)):
                if not _present(source.get(field)) or abs(float(source[field]) - float(row[field])) > 1e-6:
                    changed = True
                source[field] = row[field]
    if not changed:
        return False
    selected = {source_keys[k]['id'] for k in live_keys}
    # Live tensors are expressed in the same frame as the exported coordinates.
    live_anis = {r['id']: r for r in _cif_rows(live.block, '_atom_site_anisotrop.')}
    anis = []
    for key, row in live_keys.items():
        if row['id'] in live_anis:
            tensor = dict(live_anis[row['id']])
            tensor['id'] = source_keys[key]['id']
            # Retain just tensor values, avoiding live label/auth namespace drift.
            anis.append({k: v for k, v in tensor.items() if k == 'id' or k.startswith(('U[', 'B['))})
    _set_cif_rows(editor.block, '_atom_site_anisotrop.', anis)
    editor.report['operation'] = 'synchronize checked ChimeraX coordinates'
    chains = {r['auth_asym_id'] for r in editor.atoms if r['id'] in selected}
    editor.apply_copies([{'chains': {c: c for c in chains}, 'atom_ids': selected}])
    editor.report['warnings'].append('Live coordinates changed: secondary structure and map fit require review; automatic DSSP changes are not imported as sheet topology.')
    return True


def sync_structure_from_chimerax(session, model_id, path):
    model = _live_structure(session, model_id)
    live = structure_from_chimerax(session, model)
    editor = StructureEditor(path)
    if merge_structure_coordinates(editor, live):
        report = editor.write(path, format='cif' if editor.is_cif else 'pdb')
        if isinstance(path, StructureBuffer):
            model._amyloid_structure = StructureBuffer(path.name, path.text, path.report)
        session.logger.info(f"Synchronized {report['atoms']} atoms with the working structure.")


def _live_structure(session, model_id):
    from chimerax.atomic import AtomicStructure
    model = next((m for m in session.models.list(type=AtomicStructure)
                  if m.id_string == model_id), None)
    if model is None:
        raise StructureEditError('The working model is no longer open in ChimeraX.')
    return model


def structure_from_chimerax(session, model):
    from chimerax.mmcif import mmcif_write
    mmcif_write._set_standard_residues()
    position = model.scene_position
    transform = None if position.is_identity() else position
    with StringIO() as stream:
        mmcif_write.save_structure(session, stream, [model], [transform], set(),
                                  False, False, True, False, False, False)
        return StructureEditor.from_string(stream.getvalue(),
                                           source_name=f'ChimeraX model #{model.id_string}')


def create_structure_working_copy(session, model, path, format):
    source = getattr(model, '_amyloid_structure', None)
    if source is None:
        source = getattr(model, 'filename', None)
    live = structure_from_chimerax(session, model)
    if isinstance(source, StructureBuffer) or (isinstance(source, str) and os.path.isfile(source)):
        editor = StructureEditor(source)
        # Validate against current live state before copying source metadata.
        selected = select_live_ensemble_member(editor, live, getattr(model, 'scene_position', None))
        changed = merge_structure_coordinates(editor, live)
        if not changed:
            editor.rename({})
        report = editor.write(path, format='cif' if format == 'mmcif' else 'pdb')
        if selected:
            selection = report['ensemble_selection']
            session.logger.info(f"Loaded selected ensemble member {selection['selected_source_model']} of {len(selection['source_models'])} as one working structure ({report['atoms']} atoms). Other ensemble members are unchanged.")
    else:
        # Without a source document, only the live model metadata is available.
        editor = live
        editor.rename({})
        editor.report['warnings'].append('No accessible original source; metadata already discarded by ChimeraX cannot be recovered.')
        report = editor.write(path, format='cif' if format == 'mmcif' else 'pdb')
    for issue in report.get('source_issues', []):
        session.logger.warning(issue['warning'])
    return report


def update_secondary_structure_from_chimerax(session, model_id, path):
    from chimerax.core.commands import run
    model = _live_structure(session, model_id)
    run(session, f'dssp #{model_id}')
    editor, live = StructureEditor(path), structure_from_chimerax(session, model)
    source = {_atom_key(a):a for a in editor.atoms}
    labels = {}
    for atom in live.atoms:
        key = _atom_key(atom)
        if key not in source:
            raise StructureEditError('DSSP snapshot atom identities differ from the working file.')
        labels[atom['label_asym_id']] = source[key]['label_asym_id']
    group = {'labels':labels,'chains':{c:c for c in editor.chains}}
    for category in ('_struct_conf.','_struct_conf_type.','_struct_sheet.','_struct_sheet_range.','_struct_sheet_order.','_pdbx_struct_sheet_hbond.'):
        rows = [editor._remap_row(row,group) for row in _cif_rows(live.block,category)]
        if any(row is None for row in rows):
            raise StructureEditError('DSSP annotation contains an unmappable chain reference.')
        _set_cif_rows(editor.block,category,rows)
    editor.rename({})
    editor.report['operation'] = 'explicit ChimeraX DSSP secondary structure assignment'
    editor.report['warnings'].append('Secondary structure assigned by ChimeraX DSSP; this is computational annotation, not experimental validation.')
    editor.write(path,format='cif' if editor.is_cif else 'pdb')


def _attach_preserved_annotations(model, source):
    editor = StructureEditor(source)
    categories = editor.unrecognized_cif_categories()
    records = editor.unrecognized_pdb_records
    # The native PDB saver carries arbitrary model metadata entries through.
    by_record = defaultdict(list)
    for line in records:
        by_record[line[:6].strip()].append(line)
    for record, lines in by_record.items():
        model.set_metadata_entry(record, lines)
    if categories:
        remarks = [line for line in model.metadata.get('REMARK', [])
                   if not line.startswith('REMARK 999 AMYLOID_CIF ')]
        model.set_metadata_entry('REMARK', remarks + _cif_annotation_records(categories))
    if records:
        categories['_amyloid_unrecognized_pdb.'] = [
            {'id': str(i), 'text': line} for i, line in enumerate(records, 1)]
    model._amyloid_unrecognized_categories = categories
    if not categories:
        return
    from chimerax.mmcif import mmcif_write
    native = mmcif_write.save_structure
    if getattr(native, '_amyloid_preserves_annotations', False):
        return

    def save_with_annotations(session, file, models, *args, **kwargs):
        result = native(session, file, models, *args, **kwargs)
        preserved = {}
        for saved_model in models:
            for category, rows in getattr(saved_model, '_amyloid_unrecognized_categories', {}).items():
                if category not in preserved:
                    preserved[category] = deepcopy(rows)
                elif preserved[category] != rows:
                    preserved[category].extend(deepcopy(rows))
        if preserved:
            gemmi, _ = _structure_dependencies()
            block = gemmi.cif.Block('annotations')
            for category, rows in preserved.items():
                _set_cif_rows(block, category, rows)
            # Append categories to the current data block, before the next model.
            file.write('\n# Original annotations retained without remapping.\n' +
                       block.as_string().split('\n', 1)[1])
        return result

    save_with_annotations._amyloid_preserves_annotations = True
    mmcif_write.save_structure = save_with_annotations


def open_structure_buffer(session, source):
    if source.format == 'pdb':
        from chimerax.pdb import open_pdb
        with StringIO(source.text) as stream:
            models, _ = open_pdb(session, stream, file_name=source.name, log_info=False)
    else:
        from chimerax.atomic import AtomicStructure
        from chimerax.mmcif import mmcif, _mmcif
        if not mmcif._initialized:
            mmcif._initialize(session)
        # The native buffer parser requires the final ignore_styling argument.
        data = source.text.encode('utf-8')
        try:
            pointers = _mmcif.parse_mmCIF_buffer(data, mmcif._additional_categories,
                                                session.logger, False, True, False)
        except _mmcif.error as exc:
            if 'PDBx/mmCIF styling lost' not in str(exc):
                raise
            # Retry non-fixed-column CIF with the general in-memory parser.
            pointers = _mmcif.parse_mmCIF_buffer(data, mmcif._additional_categories,
                                                session.logger, False, True, True)
        models = [AtomicStructure(session, name=source.name, c_pointer=p, log_info=False)
                  for p in pointers]
        for model in models:
            model.is_mmcif = True
            if '_missing_poly_seq' in model.metadata:
                model.set_metadata_entry('_missing_poly_seq', None)
            combine = getattr(model, 'combine_sym_atoms', None)
            if combine is not None:
                combine()
    if not models:
        raise StructureEditError('No atomic structure could be read from the working copy.')
    try:
        for model in models:
            # Attach in-memory provenance so later tool instances can reuse source metadata.
            model._amyloid_structure = StructureBuffer(source.name, source.text, source.report)
            _attach_preserved_annotations(model, source)
        session.models.add(models)
    except Exception:
        for model in models:
            model.delete()
        raise
    return models


try:
    from PyQt6.QtWidgets import (QApplication, QWidget, QLabel, QLineEdit, 
                                 QPushButton, QTextEdit, QMessageBox, QFileDialog, 
                                 QVBoxLayout, QHBoxLayout, QSizePolicy, QSlider,
                                 QSplitter, QSpinBox, QGroupBox, QFormLayout, 
                                 QCheckBox, QTableWidget, QHeaderView, QTableWidgetItem, 
                                 QTextBrowser, QDialog, QComboBox, QTabWidget)
    from PyQt6.QtCore import Qt, QTimer
    from chimerax.core.tools import ToolInstance
    from chimerax.core.commands import run
    
except ModuleNotFoundError as exc:
    if not (exc.name.startswith('PyQt6') or exc.name.startswith('chimerax')):
        raise
    ToolInstance = object

toplevel_windows = []

def show_error_message(message):
    msg = QMessageBox()
    msg.setIcon(QMessageBox.Icon.Critical)
    msg.setWindowTitle("Error")
    msg.setText(message)
    msg.exec()

def open_chain_modifier():
    try:
        import numpy as np
        import matplotlib
        matplotlib.use('qtagg')
        from PyQt6.QtGui import QDoubleValidator, QIntValidator
        from matplotlib.backends.backend_qtagg import FigureCanvasQTAgg as FigureCanvas
        from matplotlib.figure import Figure
    except ImportError as e:
        print(f"Error loading Matplotlib: {e}")
        return

    def generate_chain_id(index):
        return structure_chain_id(index)

    class TrimDialog(QDialog):
        def __init__(self, layers, full_cif_data, original_filename, parent=None):
            super().__init__(parent)
            self.setWindowTitle("Trim Layers")
            self.resize(800, 600)
            self.layers = layers  
            self.full_cif_data = full_cif_data 
            self.original_filename = original_filename
            self.initUI()
            qt_app = QApplication.instance()
            if qt_app:
                self.setStyleSheet(qt_app.styleSheet())

        def initUI(self):
            layout = QVBoxLayout()
            layout.addWidget(QLabel("Select complete chain-unit groups to KEEP\nA chain may span multiple physical layers. Chains remain intact and are renamed in the new model."))
            self.table = QTableWidget()
            self.table.setColumnCount(4)
            self.table.setHorizontalHeaderLabels(["Keep", "Physical Layers", "Original Chains", "Residues"])
            self.table.horizontalHeader().setSectionResizeMode(QHeaderView.ResizeMode.Stretch)
            self.table.horizontalHeader().setSectionResizeMode(0, QHeaderView.ResizeMode.ResizeToContents)
            self.table.horizontalHeader().setSectionResizeMode(1, QHeaderView.ResizeMode.ResizeToContents)
            self.populate_table()
            layout.addWidget(self.table)
            btn_layout = QHBoxLayout(); btn_layout.addStretch()
            self.btn_save = QPushButton("Create Trimmed Model")
            self.btn_save.clicked.connect(self.save_cif)
            btn_layout.addWidget(self.btn_save)
            layout.addLayout(btn_layout)
            self.setLayout(layout)

        def update_data(self, layers, full_cif_data, filename):
            self.layers = layers
            self.full_cif_data = full_cif_data
            self.original_filename = filename
            self.populate_table()

        def populate_table(self):
            self.table.setRowCount(len(self.layers))
            self.table.setFocusPolicy(Qt.FocusPolicy.NoFocus) 
            self.table.setSelectionMode(QTableWidget.SelectionMode.NoSelection) 

            for i, layer in enumerate(self.layers):
                widget = QWidget(); chk = QCheckBox()
                layout = QHBoxLayout(widget); layout.addWidget(chk)
                layout.setAlignment(Qt.AlignmentFlag.AlignCenter); layout.setContentsMargins(0, 0, 0, 0)
                self.table.setCellWidget(i, 0, widget)

                rank_item = QTableWidgetItem(layer.get('layer_label', str(i + 1)))
                rank_item.setTextAlignment(Qt.AlignmentFlag.AlignCenter)
                rank_item.setFlags(Qt.ItemFlag.ItemIsEnabled) 
                self.table.setItem(i, 1, rank_item)

                chains_item = QTableWidgetItem(", ".join(layer['chains']))
                chains_item.setFlags(Qt.ItemFlag.ItemIsEnabled)
                self.table.setItem(i, 2, chains_item)

                res_item = QTableWidgetItem(str(layer['residues']))
                res_item.setFlags(Qt.ItemFlag.ItemIsEnabled)
                self.table.setItem(i, 3, res_item)

        def save_cif(self):
            selected = []
            for i in range(self.table.rowCount()):
                widget = self.table.cellWidget(i, 0)
                checkbox = widget.findChild(QCheckBox)
                if checkbox and checkbox.isChecked():
                    for cid in self.layers[i]['chains']: 
                        selected.append(cid)
            
            if not selected:
                QMessageBox.warning(self, "Warning", "No layers selected.")
                return
            
            mapping = {old: generate_chain_id(idx) for idx, old in enumerate(selected)}
            name_part, ext = os.path.splitext(os.path.basename(self.original_filename))
            
            path = StructureBuffer(f"{name_part}_trimmed.cif")
            
            try:
                self.write_cif(path, mapping)
                
                parent_widget = self.parent()
                if parent_widget:
                    open_structure_buffer(parent_widget.session, path)
                    
                QMessageBox.information(self, "Success", f"Trimmed model loaded into ChimeraX!\n(Renamed {len(mapping)} chains)")
                self.accept()
            except Exception as e: 
                QMessageBox.critical(self, "Error", str(e))

        def write_cif(self, path, mapping):
            editor = StructureEditor(self.full_cif_data['source'])
            editor.trim(set(mapping), rename=mapping)
            editor.write(path, format='cif')


    class Mol3DCanvas(FigureCanvas):
        def __init__(self, parent=None, width=3, height=4, dpi=100):
            self.fig = Figure(figsize=(width, height), dpi=dpi)
            self.fig.set_facecolor('#2C2C2C')
            self.axes = self.fig.add_subplot(111, projection='3d')
            self.axes.set_facecolor('#2C2C2C')
            self.axes.axis('off') 
            self.fig.subplots_adjust(left=0, right=1, bottom=0, top=1)
            super().__init__(self.fig)
            self.setParent(parent)
            self.setSizePolicy(QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Expanding)
            self.setMinimumSize(200, 200)
            self.updateGeometry()
            self.global_center = (0, 0, 0)
            self.max_range = 10.0
            self.zoom_level = 50 
            self.mpl_connect('scroll_event', self.on_scroll)
            self.slider_callback = None 

        def on_scroll(self, event):
            if event.inaxes != self.axes: return
            step = 5 if event.button == 'up' else -5
            new_val = max(1, min(100, self.zoom_level + step))
            self.zoom_level = new_val
            if self.slider_callback:
                self.slider_callback(new_val)
            self.apply_zoom()

        def set_zoom_from_slider(self, value):
            self.zoom_level = value
            self.apply_zoom()

        def apply_zoom(self):
            if self.max_range == 0: return
            factor = 2.0 - ((self.zoom_level / 100.0) * 1.9)
            if factor < 0.05: factor = 0.05
            radius = self.max_range * factor
            gx, gy, gz = self.global_center
            self.axes.set_xlim(gx - radius, gx + radius)
            self.axes.set_ylim(gy - radius, gy + radius)
            self.axes.set_zlim(gz - radius, gz + radius)
            self.draw()

        def plot_chains(self, chains_data, label_ids=None, waters=None, ions=None):
            self.axes.clear()
            self.axes.axis('off')
            if not chains_data:
                self.draw()
                return
                
            if label_ids is None: label_ids = set()
            chain_ids = sorted(list(chains_data.keys()))
            
            all_points = []
            calc_targets = label_ids if label_ids else chain_ids
            for cid in calc_targets:
                if cid in chains_data: all_points.extend(chains_data[cid])
            
            if not all_points: return

            xs_all = [p[0] for p in all_points]
            ys_all = [p[1] for p in all_points]
            zs_all = [p[2] for p in all_points]
            gx, gy, gz = sum(xs_all)/len(xs_all), sum(ys_all)/len(ys_all), sum(zs_all)/len(zs_all)
            self.global_center = (gx, gy, gz)

            max_dist = max([math.sqrt((p[0]-gx)**2 + (p[1]-gy)**2 + (p[2]-gz)**2) for p in all_points]) if all_points else 0
            self.max_range = max_dist if max_dist > 0 else 10.0

            for cid in chain_ids:
                coords = chains_data[cid]
                if not coords: continue
                xs, ys, zs = zip(*coords)
                color = "#FFFFFF" 
                self.axes.plot(xs, ys, zs, c=color, linewidth=1.5, alpha=0.6)
                if cid in label_ids and len(coords) > 0:
                    lx, ly, lz = coords[0]
                    self.axes.text(lx, ly, lz, f"{cid}", color='black', fontsize=11, weight='bold', horizontalalignment='center', verticalalignment='center', bbox=dict(facecolor='white', alpha=0.9, edgecolor='#333', boxstyle='round,pad=0.3'))

            if waters:
                wx, wy, wz = zip(*waters)
                self.axes.scatter(wx, wy, wz, c='lightblue', s=10, alpha=0.8, edgecolors='none')
            if ions:
                ix, iy, iz = zip(*ions)
                self.axes.scatter(ix, iy, iz, c='lightpink', s=15, alpha=0.9, edgecolors='none')

            self.apply_zoom()

    class CIFMol3DCanvas(FigureCanvas):
        def __init__(self, parent=None, width=3, height=4, dpi=100):
            self.fig = Figure(figsize=(width, height), dpi=dpi)
            self.fig.set_facecolor('#2C2C2C')
            self.axes = self.fig.add_subplot(111, projection='3d')
            self.axes.set_facecolor('#2C2C2C')
            self.axes.axis('off')
            self.fig.subplots_adjust(left=0, right=1, bottom=0, top=1)
            super().__init__(self.fig)
            self.setParent(parent)
            self.setSizePolicy(QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Expanding)
            self.setMinimumSize(200, 200)
            self.updateGeometry()
            self.global_center = (0,0,0); self.max_range = 10.0; self.zoom_level = 50
            self.mpl_connect('scroll_event', self.on_scroll)
            self.slider_callback = None

        def on_scroll(self, event):
            if event.inaxes != self.axes: return
            step = 5 if event.button == 'up' else -5
            self.zoom_level = max(1, min(100, self.zoom_level + step))
            if self.slider_callback: self.slider_callback(self.zoom_level)
            self.apply_zoom()

        def set_zoom_from_slider(self, value):
            self.zoom_level = value
            self.apply_zoom()

        def apply_zoom(self):
            if self.max_range == 0: return
            factor = 2.0 - ((self.zoom_level / 100.0) * 1.9)
            if factor < 0.05: factor = 0.05
            
            gx, gy, gz = self.global_center
            
            if hasattr(self, 'data_ranges'):
                dx, dy, dz = self.data_ranges
                self.axes.set_xlim(gx - (dx/2)*factor, gx + (dx/2)*factor)
                self.axes.set_ylim(gy - (dy/2)*factor, gy + (dy/2)*factor)
                self.axes.set_zlim(gz - (dz/2)*factor, gz + (dz/2)*factor)
            else:
                r = self.max_range * factor
                self.axes.set_xlim(gx-r, gx+r); self.axes.set_ylim(gy-r, gy+r); self.axes.set_zlim(gz-r, gz+r)
            self.draw()

        def plot_chains(self, chains_data, label_ids=None, preserve_view=False, waters=None, ions=None):
            elev, azim = self.axes.elev, self.axes.azim
            self.axes.clear()
            self.axes.axis('off')
            
            if not chains_data: 
                self.draw()
                return
                
            if label_ids is None: label_ids = set()
            cids = sorted(chains_data.keys())
            
            if not preserve_view or self.max_range == 10.0:
                all_p = [p for cid in cids if cid in chains_data for p in chains_data[cid]]
                if all_p:
                    xs, ys, zs = zip(*all_p)
                    gx, gy, gz = sum(xs)/len(xs), sum(ys)/len(ys), sum(zs)/len(zs)
                    self.global_center = (gx, gy, gz)
                    
                    dx = max(xs) - min(xs)
                    dy = max(ys) - min(ys)
                    dz = max(zs) - min(zs)
                    
                    self.data_ranges = (max(dx, 1)*1.1, max(dy, 1)*1.1, max(dz, 1)*1.1)
                    
                    try:
                        self.axes.set_box_aspect(self.data_ranges)
                    except AttributeError:
                        pass
                    
                    self.max_range = max([math.sqrt((p[0]-gx)**2+(p[1]-gy)**2+(p[2]-gz)**2) for p in all_p]) or 10.0
            
            for cid in cids:
                coords = chains_data[cid]
                if not coords: continue
                xs, ys, zs = zip(*coords)
                col, w, a = ("#98c379", 3.0, 0.9) if cid in label_ids else ("#FFFFFF", 1.0, 0.3)
                self.axes.plot(xs, ys, zs, c=col, linewidth=w, alpha=a)
                if cid in label_ids:
                    self.axes.text(coords[0][0], coords[0][1], coords[0][2], f"{cid}", color='black', fontsize=11, weight='bold', bbox=dict(facecolor='white', alpha=0.9, edgecolor='#333', boxstyle='round,pad=0.3'))
            
            if waters:
                wx, wy, wz = zip(*waters)
                self.axes.scatter(wx, wy, wz, c='lightblue', s=10, alpha=0.8, edgecolors='none')
            if ions:
                ix, iy, iz = zip(*ions)
                self.axes.scatter(ix, iy, iz, c='lightpink', s=15, alpha=0.9, edgecolors='none')
                
            if preserve_view: self.axes.view_init(elev=elev, azim=azim)
            self.apply_zoom()

    class CIFLayerIdentifier(QWidget):
        def __init__(self, tool_instance, session):
            super().__init__()
            self.tool_instance = tool_instance
            self.session = session
            self._last_model_ids = set()
            
            
            self.working_structure = None; self.layers_result = []; self.full_cif_data = {}
            self.initUI()
            qt_app = QApplication.instance()
            if qt_app: 
                self.setStyleSheet(qt_app.styleSheet())

        def initUI(self):
            main_l = QHBoxLayout()
            splitter = QSplitter(Qt.Orientation.Horizontal)
            left = QWidget(); l_lay = QVBoxLayout(left); l_lay.setSpacing(2); l_lay.setContentsMargins(4, 4, 4, 4)
            
            hb_load = QHBoxLayout()
            self.cmb_models = QComboBox()
            self.cmb_models.wheelEvent = lambda event: event.ignore()
            self.cmb_models.setSizePolicy(QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Fixed)
            self.populate_models()
            hb_load.addWidget(self.cmb_models)
            
            self.auto_refresh_timer = QTimer(self)
            self.auto_refresh_timer.timeout.connect(self.populate_models)
            self.auto_refresh_timer.timeout.connect(self.check_live_edits)
            self.auto_refresh_timer.start(1000)
            
            self.btn_load_model = QPushButton('Load Model')
            self.btn_load_model.clicked.connect(self.load_from_chimerax)
            hb_load.addWidget(self.btn_load_model)
            l_lay.addLayout(hb_load)
            
            self.lbl_status = QLabel("No file loaded"); self.lbl_status.setStyleSheet("color: gray; font-style: italic;")
            l_lay.addWidget(self.lbl_status)
            
            self.grp_ops = QGroupBox("Operations"); self.grp_ops.setEnabled(False)
            ops_l = QVBoxLayout(); ops_l.setSpacing(2); ops_l.setContentsMargins(4, 10, 4, 4)
            self.btn_trim = QPushButton("Trim Chain Units")
            self.btn_trim.setToolTip("Keep selected layers in a new model; use ChimeraX Save to save it")
            self.btn_trim.clicked.connect(self.open_trim_dialog)
            ops_l.addWidget(self.btn_trim)
            self.grp_ops.setLayout(ops_l); l_lay.addWidget(self.grp_ops)
            
            grp_p = QGroupBox("Detection Parameters"); form = QFormLayout(); form.setContentsMargins(4, 10, 4, 4); form.setVerticalSpacing(2)
            self.ed_chn = QLineEdit("4"); self.ed_chn.setValidator(QIntValidator())
            form.addRow("Chains per Layer:", self.ed_chn)
            self.ed_zmn = QLineEdit("0.0"); self.ed_zmn.setValidator(QDoubleValidator())
            form.addRow("Z Shift Min (Å):", self.ed_zmn)
            self.ed_zmx = QLineEdit("4.0"); self.ed_zmx.setValidator(QDoubleValidator())
            form.addRow("Z Shift Max (Å):", self.ed_zmx)
            self.ed_chn.editingFinished.connect(self.recalc); self.ed_zmn.editingFinished.connect(self.recalc); self.ed_zmx.editingFinished.connect(self.recalc)
            for field in (self.ed_chn, self.ed_zmn, self.ed_zmx):
                field.hide()
                form.labelForField(field).hide()
            grp_p.setLayout(form); l_lay.addWidget(grp_p)
            
            self.txt = QTextBrowser(); self.txt.setOpenLinks(False); self.txt.anchorClicked.connect(self.link_clk)
            self.txt.setStyleSheet("font-family: Consolas, monospace;")
            self.txt.setSizePolicy(QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Expanding)
            l_lay.addWidget(self.txt, 1)
            
            right = QWidget(); r_lay = QHBoxLayout(right); r_lay.setContentsMargins(0,0,0,0)
            self.cvs = CIFMol3DCanvas(self, dpi=100); r_lay.addWidget(self.cvs, 1)
            sl = QSlider(Qt.Orientation.Vertical); sl.setRange(1, 100); sl.setValue(50); sl.valueChanged.connect(lambda v: self.cvs.set_zoom_from_slider(v))
            self.cvs.slider_callback = lambda v: sl.setValue(v)
            r_lay.addWidget(sl)
            
            splitter.addWidget(left); splitter.addWidget(right)
             
            splitter.setStretchFactor(0, 5); splitter.setStretchFactor(1, 5)
            main_l.addWidget(splitter); self.setLayout(main_l)

        def sync_viewer_to_structure(self):
            if self.working_model_id and self.working_structure:
                try:
                    sync_structure_from_chimerax(self.session, self.working_model_id, self.working_structure)
                except Exception as exc:
                    self.session.logger.error('Coordinate synchronization stopped: ' + str(exc))
                    return False
            return True


        def check_live_edits(self):
            if getattr(self, '_load_in_progress', False) or getattr(self, '_reload_in_progress', False): return
            if not getattr(self, 'working_model_id', None): return
            try:
                from chimerax.atomic import AtomicStructure
                found = False
                for m in self.session.models.list(type=AtomicStructure):
                    if m.id_string == self.working_model_id:
                        found = True
                        current_atoms = len(m.atoms)
                        if getattr(self, '_last_live_atom_count', -1) != current_atoms:
                            self._last_live_atom_count = current_atoms
                            if not self.sync_viewer_to_structure():
                                return
                            self.process_cif(self.working_structure)
                        break
                
                if not found:
                    self.working_model_id = None
                    self.working_structure = None
                    self.grp_ops.setEnabled(False)
                    self.cvs.axes.clear()
                    self.cvs.axes.axis('off')
                    self.cvs.draw()
                    self.txt.append("<br><b>>>> Working model was closed in ChimeraX. Tool cleared.</b><br>")
            except Exception:
                pass

        def populate_models(self):
            from chimerax.atomic import AtomicStructure
            current_models = self.session.models.list(type=AtomicStructure)
            current_ids = {(m.id_string, id(m)) for m in current_models}
            if current_ids != self._last_model_ids:
                current_sel = self.cmb_models.currentData()
                try:
                    current_sel_id = current_sel.id_string if current_sel and not getattr(current_sel, 'deleted', False) else None
                except Exception:
                    current_sel_id = None
                    
                if not current_sel_id and getattr(self, 'working_model_id', None):
                    current_sel_id = self.working_model_id
                    
                self.cmb_models.blockSignals(True)
                self.cmb_models.clear()
                for model in current_models:
                    item_text = f"#{model.id_string} {model.name}"
                    self.cmb_models.addItem(item_text, userData=model)
                    self.cmb_models.setItemData(self.cmb_models.count() - 1, item_text, Qt.ItemDataRole.ToolTipRole)
                if current_sel_id:
                    for i in range(self.cmb_models.count()):
                        m = self.cmb_models.itemData(i)
                        if m and getattr(m, 'id_string', '') == current_sel_id:
                            self.cmb_models.setCurrentIndex(i)
                            break
                self.cmb_models.blockSignals(False)
                self._last_model_ids = current_ids

        def reload_working_model(self):
            target_id = getattr(self, 'working_model_id', None)
            
            captured_state = False
            was_cartoon_visible = False
            was_atom_visible = True
            atom_style = None
            
            try:
                from chimerax.atomic import AtomicStructure
                for m in self.session.models.list(type=AtomicStructure):
                    if m.id_string == target_id:
                        captured_state = True
                        if hasattr(m, 'residues') and len(m.residues) > 0:
                            was_cartoon_visible = any(getattr(r, 'ribbon_display', False) for r in m.residues)
                        if hasattr(m, 'atoms') and len(m.atoms) > 0:
                            was_atom_visible = any(getattr(a, 'display', True) for a in m.atoms)
                            for a in m.atoms:
                                if getattr(a, 'display', False):
                                    mode_val = getattr(a, 'draw_mode', 2)
                                    mode_str = str(mode_val).lower()
                                    if mode_val == 0 or "sphere" in mode_str: atom_style = "sphere"
                                    elif mode_val == 1 or "ball" in mode_str: atom_style = "ball"
                                    elif mode_val == 3 or "wire" in mode_str: atom_style = "wire"
                                    else: atom_style = "stick"
                                    break
                        break
            except Exception:
                pass
            
            try:
                from chimerax.core.commands import run
                
                try: run(self.session, "view name chimerax_modifier_locked_view")
                except: pass
                
                models = open_structure_buffer(self.session, self.working_structure)
                if target_id:
                    try: run(self.session, f"close #{target_id}")
                    except: pass
                    
                if models:
                    first_item = models[0]
                    first_model = first_item[0] if isinstance(first_item, (list, tuple)) else first_item
                    if hasattr(first_model, 'id_string'):
                        self.working_model_id = first_model.id_string
                        self._last_live_atom_count = len(first_model.atoms) if hasattr(first_model, 'atoms') else 0
                        run(self.session, f"color #{self.working_model_id} bypolymer")
                        
                        if captured_state:
                            if was_cartoon_visible: run(self.session, f"show #{self.working_model_id} cartoons")
                            else: run(self.session, f"hide #{self.working_model_id} cartoons")
                            
                            if was_atom_visible:
                                run(self.session, f"show #{self.working_model_id} atoms")
                                if atom_style: run(self.session, f"style #{self.working_model_id} {atom_style}")
                            else: run(self.session, f"hide #{self.working_model_id} atoms")
                
                try: run(self.session, "view chimerax_modifier_locked_view")
                except: pass
                
            except Exception as e:
                self.txt.append(f"\nModel display warning: {e}")

        def load_from_chimerax(self):
            model = self.cmb_models.currentData()
            if not model: return
            self._axis_reference = None
            self._load_in_progress = True
            try:
                self.loaded_filename = model.name
                
                name_part, _ = os.path.splitext(self.loaded_filename)
                working = StructureBuffer(f"{name_part}_modified.cif")
                
                create_structure_working_copy(self.session, model, working, "mmcif")
                self.working_structure = working
                self.grp_ops.setEnabled(True)
                
                self.working_model_id = None
                self.reload_working_model()
                
                self.lbl_status.setText(f"Loaded: {name_part}_modified.cif")
                self.txt.clear(); self.txt.append(f"Created CIF working copy from {model.name}")
                self.process_cif(working)
                
                self.populate_models()
                for i in range(self.cmb_models.count()):
                    m = self.cmb_models.itemData(i)
                    if m and getattr(m, 'id_string', '') == getattr(self, 'working_model_id', ''):
                        self.cmb_models.setCurrentIndex(i)
                        break
            except Exception as e:
                self.txt.append(f"Load Error: {e}")
            finally:
                self._load_in_progress = False

        def recalc(self):
            if self.working_structure: 
                if not self.sync_viewer_to_structure():
                    return
                QTimer.singleShot(10, lambda: self.process_cif(self.working_structure))

        def link_clk(self, url):
            try:
                idx = int(url.toString())
                if 0 <= idx < len(self.layers_result):
                    self.cvs.plot_chains(self.chains_data_plot, label_ids=set(self.layers_result[idx]['chains']), preserve_view=True, waters=getattr(self, 'current_waters', []), ions=getattr(self, 'current_ions', []))
            except: pass

        def open_trim_dialog(self):
            if not self.layers_result: return
            if hasattr(self, 'trim_dialog') and self.trim_dialog and self.trim_dialog.isVisible():
                self.trim_dialog.raise_()
                self.trim_dialog.activateWindow()
                return
            self.trim_dialog = TrimDialog(self.layers_result, self.full_cif_data, self.loaded_filename, self)
            self.trim_dialog.show()

        def parse_cif(self, path):
            editor = StructureEditor(path)
            chains, waters, ions = editor.preview()
            return {c: [r['coord'] for r in residues] for c, residues in chains.items()}, waters, ions, {'source': path}


        def process_cif(self, path):
            try:
                editor = StructureEditor(path)
                chains, waters, ions = editor.preview()
                result = detect_layers(chains)
                self.detection_result = result
                self.full_cif_data = {'source': path}
                coords = {c: [r['coord'] for r in rows] for c, rows in chains.items()}
                self.chains_data_plot = {c: [tuple(p) for p in coords[c]] for c in coords}
                self.current_waters = waters
                self.current_ions = ions
                self.layers_result = []
                groups = defaultdict(list)
                for track in result['protofilaments']:
                    k = track['layers_per_unit']
                    for position, cid in zip(track['unit_positions'], track['chains']):
                        groups[(position * k + 1, (position + 1) * k)].append(cid)
                for (first, last), cur in sorted(groups.items()):
                    self.layers_result.append({'chains': cur, 'count': len(cur),
                        'residues': sum(len(coords[c]) for c in cur),
                        'layer_label': str(first) if first == last else f'{first}–{last}'})
                h_lines = [result['orientation'].capitalize() + '<br>',
                           f"Protofilaments: {len(result['sandwiches'])}; maximum chain units: {result['detected_units']}; physical layers: {result['detected_layers']}<br>"]
                h_lines.extend('Warning: ' + escape(issue['warning']) + '<br>'
                               for issue in editor.report.get('source_issues', []))
                h_lines.extend('Note: ' + w + '<br>' for w in result['warnings'])
                high_ids = set()
                for i, lay in enumerate(self.layers_result):
                    lnk = f"<a href='{i}' style='color: #97c379; font-weight: bold; text-decoration: none;'>Layers {lay['layer_label']}</a>"
                    h_lines.append(f"{lnk}: {', '.join(lay['chains'])} [Chains, Residues]:{lay['count']}, {lay['residues']}<br>")
                    if i==0: high_ids.update(lay['chains'])
                self.txt.setHtml("<br>".join(h_lines))
                self.cvs.plot_chains(self.chains_data_plot, label_ids=high_ids, waters=waters, ions=ions)
                
                if hasattr(self, 'trim_dialog') and self.trim_dialog and self.trim_dialog.isVisible():
                    self.trim_dialog.update_data(self.layers_result, self.full_cif_data, self.loaded_filename)

            except Exception as e:
                self.layers_result = []
                self.detection_result = None
                self.txt.append(str(e))

    class PDBLayerIdentifier(QWidget):
        

        def __init__(self, tool_instance, session):
            super().__init__()
            self.tool_instance = tool_instance
            self.session = session
            self.working_model_id = None
            self._reload_in_progress = False
            self.final_sandwiches = [] 
            self.detected_layers = 0 
            self.detected_units = 0
            self.layers_per_unit = 1
            self.detection_result = None
            self._last_model_ids = set() 
            

            self.working_structure = None 
            self.original_filename_display = ""
            self.initUI()

        def initUI(self):
            main_layout = QHBoxLayout()
            left_widget = QWidget()
            layout = QVBoxLayout(left_widget)
            layout.setSpacing(2)
            layout.setContentsMargins(4, 4, 4, 4)

            hb_load = QHBoxLayout()
            self.cmb_models = QComboBox()
            self.cmb_models.wheelEvent = lambda event: event.ignore()
            self.cmb_models.setSizePolicy(QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Fixed)
            self.populate_models()
            hb_load.addWidget(self.cmb_models)
            
            self.auto_refresh_timer = QTimer(self)
            self.auto_refresh_timer.timeout.connect(self.populate_models)
            self.auto_refresh_timer.timeout.connect(self.check_live_edits)
            self.auto_refresh_timer.start(1000)
            
            self.btn_load_model = QPushButton('Load Model')
            self.btn_load_model.clicked.connect(self.load_from_chimerax)
            hb_load.addWidget(self.btn_load_model)

            self.lbl_orientation = QLabel("Undetermined")
            self.lbl_orientation.setToolTip("Parallel and antiparallel stacking are detected from coordinates.")
            hb_load.addWidget(self.lbl_orientation)
            layout.addLayout(hb_load)

            self.grp_ops = QGroupBox("Operations")
            self.grp_ops.setEnabled(False) 
            v_ops = QVBoxLayout()
            v_ops.setSpacing(2)
            v_ops.setContentsMargins(4, 10, 4, 4)

            hb_trim = QHBoxLayout()
            self.btn_trim = QPushButton("Trim/Expand Layers")
            self.btn_trim.setToolTip("Trim or expand the structure to the target number of layers")
            self.btn_trim.clicked.connect(self.action_trim)
            hb_trim.addWidget(self.btn_trim)
            hb_trim.addWidget(QLabel("Target Layers:"))
            self.edit_trim = QLineEdit("5")
            self.edit_trim.setValidator(QIntValidator(1, 500))
            self.edit_trim.setFixedWidth(50)
            hb_trim.addWidget(self.edit_trim)

            self.chk_fit_map = QCheckBox("Fit to Map")
            self.chk_fit_map.setToolTip("Fit to Map after trimming only. Expansion keeps existing coordinates and the shared helical transform unchanged.")
            self.chk_fit_map.toggled.connect(self.on_fit_map_toggled)
            self.fit_map_model_id = None
            self.fit_map_mode = None
            hb_trim.addWidget(self.chk_fit_map)

            hb_trim.addStretch()
            v_ops.addLayout(hb_trim)

            hb_pf_trim = QHBoxLayout()
            
            self.btn_remove_pf = QPushButton("Trim")
            self.btn_remove_pf.setFixedWidth(70)
            self.btn_remove_pf.setToolTip("Trim away the selected protofilament entirely")
            self.btn_remove_pf.clicked.connect(self.action_remove_pf)
            hb_pf_trim.addWidget(self.btn_remove_pf)

            self.btn_select_pf = QPushButton("Select")
            self.btn_select_pf.setFixedWidth(70)
            self.btn_select_pf.setToolTip("Left-click to select this protofilament; right-click to deselect it.")
            self.btn_select_pf.clicked.connect(self.action_select_pf)
            self.btn_select_pf.setContextMenuPolicy(Qt.ContextMenuPolicy.CustomContextMenu)
            self.btn_select_pf.customContextMenuRequested.connect(lambda _pos: self.action_select_pf(deselect=True))
            hb_pf_trim.addWidget(self.btn_select_pf)

            self.cmb_protofilaments = QComboBox()
            self.cmb_protofilaments.setSizePolicy(QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Fixed)
            hb_pf_trim.addWidget(self.cmb_protofilaments, 1)
            
            v_ops.addLayout(hb_pf_trim)

            hb_ren = QHBoxLayout()
            self.btn_rename = QPushButton("Rename Chains")
            self.btn_rename.clicked.connect(self.action_rename)
            hb_ren.addWidget(self.btn_rename)
            hb_ren.addWidget(QLabel("(Rename Chain ID from the Middle Out)"))
            v_ops.addLayout(hb_ren)

            self.chk_numbers_first = QCheckBox("Use Numbers (0-9) after Z")
            self.chk_numbers_first.setToolTip("Use only one letter/number to rename chains to avoid problems loading PDB in old softwares")
            v_ops.addWidget(self.chk_numbers_first)

            self.chk_only_numbers = QCheckBox("Use only numbers")
            self.chk_only_numbers.setToolTip("Use only numbers to rename chains")
            v_ops.addWidget(self.chk_only_numbers)
            self.chk_numbers_first.toggled.connect(lambda state: state and self.chk_only_numbers.setChecked(False))
            self.chk_only_numbers.toggled.connect(lambda state: state and self.chk_numbers_first.setChecked(False))

            self.btn_beta = QPushButton("Re-evaluate β-sheet structure")
            self.btn_beta.clicked.connect(self.action_evaluate_beta)
            self.btn_beta.setEnabled(False) 
            v_ops.addWidget(self.btn_beta)
            

            hb_res = QHBoxLayout()
            self.btn_restrain = QPushButton("Generate Restrain File")
            self.btn_restrain.clicked.connect(self.action_restrain)
            hb_res.addWidget(self.btn_restrain)
            hb_res.addWidget(QLabel("(Export .cxc File)"))
            v_ops.addLayout(hb_res)

            self.grp_ops.setLayout(v_ops)
            layout.addWidget(self.grp_ops)

            self.grp_params = QGroupBox("Detection Parameters")
            params_layout = QFormLayout()
            self.detection_form = params_layout
            params_layout.setContentsMargins(4, 10, 4, 4)
            params_layout.setVerticalSpacing(2)

            def create_param_input(default_val, tooltip_msg, is_int=False):
                line_edit = QLineEdit()
                line_edit.setText(str(default_val))
                if is_int: line_edit.setValidator(QIntValidator())
                else: line_edit.setValidator(QDoubleValidator())
                line_edit.editingFinished.connect(self.on_param_change) 
                line_edit.setToolTip(tooltip_msg) 
                return line_edit

            self.edit_z_min = create_param_input(4.6, "Minimum vertical distance to be considered as a layer stack")
            params_layout.addRow("Z Shift Min (Å):", self.edit_z_min)
            self.edit_z_max = create_param_input(5.0, "Maximum vertical distance to be considered as a layer stack")
            params_layout.addRow("Z Shift Max (Å):", self.edit_z_max)
            self.edit_xy_limit = create_param_input(3.0, "Maximum horizontal drift allowed between stacked subunits, increase this number if your amino acids extend too long from fibril axis")
            params_layout.addRow("XY Shift Limit (Å):", self.edit_xy_limit)
            self.edit_neighbor_count = create_param_input(6, "How many closest chains to search for stacking each time", is_int=True)
            params_layout.addRow("Neighbor Search Count:", self.edit_neighbor_count)
            for field in (self.edit_z_min, self.edit_z_max, self.edit_xy_limit, self.edit_neighbor_count):
                field.hide()
                params_layout.labelForField(field).hide()
            self.edit_water_z = create_param_input(4.0, "Maximum vertical distance from the protein layer centroids for a water/ion to be included")
            params_layout.addRow("Water/Ion Axial Limit:", self.edit_water_z)

            self.grp_params.setLayout(params_layout)

            self.text_area = QTextEdit()
            self.text_area.setReadOnly(True)
            self.text_area.setSizePolicy(QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Expanding)
            font = self.text_area.font()
            font.setFamily("Courier New")
            font.setPointSize(7)
            self.text_area.setFont(font)
            layout.addWidget(self.text_area, 1)

            right_widget = QWidget()
            r_layout = QVBoxLayout(right_widget)
            r_layout.setContentsMargins(0, 0, 0, 0)
            r_layout.setSpacing(2)
            r_layout.addWidget(self.grp_params)

            viewer_layout = QHBoxLayout()

            self.mol_canvas = Mol3DCanvas(self, width=3, height=4, dpi=100)
            viewer_layout.addWidget(self.mol_canvas, 1)

            self.zoom_slider = QSlider(Qt.Orientation.Vertical)
            self.zoom_slider.setRange(1, 100)
            self.zoom_slider.setValue(50)
            self.zoom_slider.setTickPosition(QSlider.TickPosition.TicksRight)
            self.zoom_slider.setTickInterval(10)
            self.zoom_slider.valueChanged.connect(lambda v: self.mol_canvas.set_zoom_from_slider(v))
            self.mol_canvas.slider_callback = lambda v: self.zoom_slider.setValue(v)
            viewer_layout.addWidget(self.zoom_slider)

            r_layout.addLayout(viewer_layout)

            splitter = QSplitter(Qt.Orientation.Horizontal)
            splitter.addWidget(left_widget)
            splitter.addWidget(right_widget)
            
            splitter.setStretchFactor(0, 4)
            splitter.setStretchFactor(1, 6)

            main_layout.addWidget(splitter)
            self.setLayout(main_layout)

        def on_fit_map_toggled(self, checked):
            if not checked:
                self.fit_map_model_id = None
                self.fit_map_mode = None
                return

            vols = [m for m in self.session.models.list() if hasattr(m, 'data')]
            if not vols:
                QMessageBox.warning(self, "No Maps", "No map/volume models currently loaded in ChimeraX.")
                self.chk_fit_map.blockSignals(True)
                self.chk_fit_map.setChecked(False)
                self.chk_fit_map.blockSignals(False)
                return

            dialog = QDialog(self)
            dialog.setWindowTitle("Fit to Map Options")
            layout = QVBoxLayout(dialog)
            
            layout.addWidget(QLabel("Select the target map (.mrc):"))
            combo = QComboBox()
            for v in vols:
                combo.addItem(f"#{v.id_string} {v.name}", userData=v.id_string)
            layout.addWidget(combo)
            
            layout.addWidget(QLabel("Select fitting mode:"))
            from PyQt6.QtWidgets import QRadioButton
            radio_whole = QRadioButton("Fit as a whole")
            radio_whole.setToolTip("Performs a global Fit-to-Map after trimming. Automatic map fitting is skipped during expansion to preserve the existing structure.")
            radio_chain = QRadioButton("Fit as per chain")
            radio_chain.setToolTip("Fits individual chains after trimming. Automatic map fitting is skipped during expansion to preserve the existing structure.")
            radio_whole.setChecked(True)
            
            layout.addWidget(radio_whole)
            layout.addWidget(radio_chain)
            
            btn_box = QHBoxLayout()
            btn_ok = QPushButton("OK")
            btn_ok.clicked.connect(dialog.accept)
            btn_cancel = QPushButton("Cancel")
            btn_cancel.clicked.connect(dialog.reject)
            btn_box.addWidget(btn_ok)
            btn_box.addWidget(btn_cancel)
            layout.addLayout(btn_box)

            if dialog.exec() == QDialog.DialogCode.Accepted:
                self.fit_map_model_id = combo.currentData()
                self.fit_map_mode = "whole" if radio_whole.isChecked() else "chain"
            else:
                self.chk_fit_map.blockSignals(True)
                self.chk_fit_map.setChecked(False)
                self.chk_fit_map.blockSignals(False)
        
        def sync_viewer_to_structure(self):
            if self.working_model_id and self.working_structure:
                try:
                    sync_structure_from_chimerax(self.session, self.working_model_id, self.working_structure)
                except Exception as exc:
                    self.session.logger.error('Coordinate synchronization stopped: ' + str(exc))
                    return False
            return True


        def on_param_change(self):
            if self.working_structure:
                if not self.sync_viewer_to_structure():
                    return
                self.text_area.append("\n--- Parameters changed: Re-calculating ---")
                QTimer.singleShot(10, lambda: self.process_pdb(self.working_structure))

        def get_param(self, widget, default):
            try: return float(widget.text())
            except ValueError: return default
        
        def get_kept_heteroatoms(self, keep_layer_indices, xy_limit, water_z_limit):
            if not self.working_structure or not self.final_sandwiches:
                return set()
            editor = StructureEditor(self.working_structure)
            core = {c for s in self.final_sandwiches for c in s}
            kept = {s[i] for s in self.final_sandwiches for i in keep_layer_indices if 0 <= i < len(s)}
            axis = (getattr(self, 'detection_result', None) or {}).get('axis')
            owners = editor.associated_residues(core, water_z_limit, axis)
            return {residue for residue, owner in owners.items() if owner in kept}


        def check_live_edits(self):
            if getattr(self, '_load_in_progress', False) or getattr(self, '_reload_in_progress', False): return
            if not getattr(self, 'working_model_id', None): return
            try:
                from chimerax.atomic import AtomicStructure
                found = False
                for m in self.session.models.list(type=AtomicStructure):
                    if m.id_string == self.working_model_id:
                        found = True
                        current_atoms = len(m.atoms)
                        if getattr(self, '_last_live_atom_count', -1) != current_atoms:
                            self._last_live_atom_count = current_atoms
                            if not self.sync_viewer_to_structure():
                                return
                            self.process_pdb(self.working_structure)
                        break
                
                if not found:
                    self.working_model_id = None
                    self.working_structure = None
                    self.grp_ops.setEnabled(False)
                    self.mol_canvas.axes.clear()
                    self.mol_canvas.axes.axis('off')
                    self.mol_canvas.draw()
                    self.text_area.append("\n>>> Working model was closed in ChimeraX. Tool cleared.")
            except Exception:
                pass

        def populate_models(self):
            from chimerax.atomic import AtomicStructure
            current_models = self.session.models.list(type=AtomicStructure)
            current_ids = {(m.id_string, id(m)) for m in current_models}
            
            if current_ids != self._last_model_ids:
                current_sel = self.cmb_models.currentData()
                try:
                    current_sel_id = current_sel.id_string if current_sel and not getattr(current_sel, 'deleted', False) else None
                except Exception:
                    current_sel_id = None
                
                if not current_sel_id and getattr(self, 'working_model_id', None):
                    current_sel_id = self.working_model_id
                
                self.cmb_models.blockSignals(True)
                self.cmb_models.clear()
                
                for model in current_models:
                    item_text = f"#{model.id_string} {model.name}"
                    self.cmb_models.addItem(item_text, userData=model)
                    self.cmb_models.setItemData(self.cmb_models.count() - 1, item_text, Qt.ItemDataRole.ToolTipRole)
                
                if current_sel_id:
                    for i in range(self.cmb_models.count()):
                        m = self.cmb_models.itemData(i)
                        if m and getattr(m, 'id_string', '') == current_sel_id:
                            self.cmb_models.setCurrentIndex(i)
                            break
                            
                self.cmb_models.blockSignals(False)
                self._last_model_ids = current_ids

        def load_from_chimerax(self):
            model = self.cmb_models.currentData()
            if not model: return
            self._load_in_progress = True
            try:
                self._axis_reference = None
                self.original_filename_display = model.name
                
                name_part, ext = os.path.splitext(self.original_filename_display)
                is_cif = ext.lower() in ('.cif', '.mmcif')
                save_ext = ".cif" if is_cif else ".pdb"
                save_format = "mmcif" if is_cif else "pdb"
                
                working = StructureBuffer(f"{name_part}_modified{save_ext}")
                
                create_structure_working_copy(self.session, model, working, save_format)
                self.working_structure = working
                self.grp_ops.setEnabled(True)
                
                self.working_model_id = None 
                
                self.reload_working_model()
                self.text_area.clear()
                self.text_area.append(f"Created {'mmCIF' if is_cif else 'PDB'} working copy from {model.name}")
                self.check_spatial_and_id_duplicates(self.working_structure)
                self.process_pdb(self.working_structure)
                
                self.populate_models()
                for i in range(self.cmb_models.count()):
                    m = self.cmb_models.itemData(i)
                    if m and getattr(m, 'id_string', '') == self.working_model_id:
                        self.cmb_models.setCurrentIndex(i)
                        break
                        
            except Exception as e:
                self.text_area.append(f"Load Error: {e}")
            finally:
                self._load_in_progress = False

        def reload_working_model(self):
            if not self.working_structure:
                return False

            target_id = self.working_model_id
            
            captured_state = False
            was_cartoon_visible = False
            was_atom_visible = True
            atom_style = None
            
            try:
                from chimerax.atomic import AtomicStructure
                for m in self.session.models.list(type=AtomicStructure):
                    if m.id_string == target_id:
                        captured_state = True
                        if hasattr(m, 'residues') and len(m.residues) > 0:
                            was_cartoon_visible = any(getattr(r, 'ribbon_display', False) for r in m.residues)
                        if hasattr(m, 'atoms') and len(m.atoms) > 0:
                            was_atom_visible = any(getattr(a, 'display', True) for a in m.atoms)
                            for a in m.atoms:
                                if getattr(a, 'display', False):
                                    mode_val = getattr(a, 'draw_mode', 2)
                                    mode_str = str(mode_val).lower()
                                    if mode_val == 0 or "sphere" in mode_str: atom_style = "sphere"
                                    elif mode_val == 1 or "ball" in mode_str: atom_style = "ball"
                                    elif mode_val == 3 or "wire" in mode_str: atom_style = "wire"
                                    else: atom_style = "stick"
                                    break
                        break
            except Exception:
                pass
            
            from chimerax.core.commands import run
            from chimerax.atomic import AtomicStructure

            self._reload_in_progress = True
            try:
                try: run(self.session, "view name chimerax_modifier_locked_view")
                except: pass

                # Validate the replacement before closing the current model.
                try:
                    models = open_structure_buffer(self.session, self.working_structure)
                except Exception as e:
                    self.text_area.append(f"\nModel reload failed; previous model kept: {e}")
                    return False

                def first_atomic_model(value):
                    if isinstance(value, AtomicStructure):
                        return value
                    if isinstance(value, (list, tuple)):
                        for item in value:
                            model = first_atomic_model(item)
                            if model is not None:
                                return model
                    return None

                first_model = first_atomic_model(models)
                if first_model is None or not hasattr(first_model, 'id_string'):
                    self.text_area.append("\nModel reload failed; previous model kept: no atomic structure was opened.")
                    return False

                new_model_id = first_model.id_string
                self.working_model_id = new_model_id
                self._last_live_atom_count = len(first_model.atoms) if hasattr(first_model, 'atoms') else 0

                if target_id and target_id != new_model_id:
                    try: run(self.session, f"close #{target_id}")
                    except Exception as e:
                        self.text_area.append(f"\nPrevious working model could not be closed: {e}")

                try:
                    run(self.session, f"color #{new_model_id} bypolymer")

                    if captured_state:
                        if was_cartoon_visible: run(self.session, f"show #{new_model_id} cartoons")
                        else: run(self.session, f"hide #{new_model_id} cartoons")

                        if was_atom_visible:
                            run(self.session, f"show #{new_model_id} atoms")
                            if atom_style: run(self.session, f"style #{new_model_id} {atom_style}")
                        else: run(self.session, f"hide #{new_model_id} atoms")
                except Exception as e:
                    self.text_area.append(f"\nModel display warning: {e}")

                return True
            finally:
                try: run(self.session, "view chimerax_modifier_locked_view")
                except: pass
                self._reload_in_progress = False

        def action_remove_pf(self):
            if not self.working_structure or not self.final_sandwiches: return
            if self.cmb_protofilaments.count() == 0: return
            
            if not self.sync_viewer_to_structure():
                return
            try:
                pf_index = self.cmb_protofilaments.currentData()
                if pf_index is None or pf_index < 0 or pf_index >= len(self.final_sandwiches): return
                
                total_layers = self.detected_units
                keep_chains = set()
                
                for i, sandwich in enumerate(self.final_sandwiches):
                    if i != pf_index:
                        for chain_id in sandwich:
                            keep_chains.add(chain_id)
                            
                if not keep_chains:
                    QMessageBox.warning(self, "Warning", "Cannot remove the last remaining protofilament.")
                    return
                
                mid_chain = self.final_sandwiches[pf_index][len(self.final_sandwiches[pf_index])//2]
                self.text_area.append(f"Removing Protofilament with middle chain {mid_chain}...")
                
                edited = StructureBuffer(self.working_structure.name)
                valid_layer_indices = range(total_layers)
                xy_lim = self.get_param(self.edit_xy_limit, 3.0)
                water_z_limit = self.get_param(self.edit_water_z, 4.0)
                
                keep_het_lines = self.get_kept_heteroatoms(valid_layer_indices, xy_lim, water_z_limit)
                
                if self.working_structure.format == 'cif':
                    self.write_trimmed_cif(self.working_structure, edited, keep_chains, keep_het_lines)
                else:
                    self.write_trimmed_pdb(self.working_structure, edited, keep_chains, keep_het_lines)
                    
                self.working_structure.copy_from(edited)
                self.text_area.append("\n>>>> Applied Protofilament Trimming")
                
                if not self.reload_working_model():
                    self.text_area.append("\nThe edit was prepared, but the previous ChimeraX model was kept because the replacement could not be opened.")
                    return
                self.process_pdb(self.working_structure)
                
            except Exception as e:
                QMessageBox.critical(self, "Error", f"Action failed: {e}")

        def action_select_pf(self, checked=False, *, deselect=False):
            if not getattr(self, 'working_model_id', None) or not getattr(self, 'final_sandwiches', []): return
            if self.cmb_protofilaments.count() == 0: return
            
            try:
                pf_index = self.cmb_protofilaments.currentData()
                if pf_index is None or pf_index < 0 or pf_index >= len(self.final_sandwiches): return
                
                target_sandwich = self.final_sandwiches[pf_index]
                if not target_sandwich: return
                
                chain_str = ",".join(target_sandwich)
                from chimerax.core.commands import run
                
                if deselect:
                    run(self.session, f"select subtract #{self.working_model_id}/{chain_str}")
                else:
                    run(self.session, f"select #{self.working_model_id}/{chain_str}")
                verb = 'Deselected' if deselect else 'Selected'
                self.text_area.append(f"\n>>>> {verb} Protofilament with chains: {chain_str}")
                
            except Exception as e:
                QMessageBox.critical(self, "Error", f"Selection failed: {e}")

        def action_rename(self):
            if not self.working_structure: return
            if not self.sync_viewer_to_structure():
                return
            try:
                old_sandwiches = self.final_sandwiches
                mapping = self.create_renaming_mapping(old_sandwiches)
                
                edited = StructureBuffer(self.working_structure.name)
                if self.working_structure.format == 'cif':
                    self.write_renamed_cif(self.working_structure, edited, self.final_sandwiches)
                else:
                    self.write_renamed_pdb(self.working_structure, edited, self.final_sandwiches)
                self.working_structure.copy_from(edited)
                
                self.text_area.append("\n>>>> Applied Renaming")
                self.reload_working_model()
                reference = getattr(self, '_axis_reference', None)
                if reference:
                    reference['chains'] = {mapping.get(c, c): rows for c, rows in reference['chains'].items()}
                self.process_pdb(self.working_structure, old_sandwiches=old_sandwiches, rename_map=mapping)
            except Exception as e: QMessageBox.critical(self, "Error", f"Rename failed: {e}")

        class ExpandDialog(QDialog):
            def __init__(self, current_layers, parent=None):
                super().__init__(parent)
                self.setWindowTitle("Expand Layers")
                self.setModal(False)
                self.setWindowModality(Qt.WindowModality.NonModal)
                self.current_layers = current_layers
                self.unit_layers = getattr(parent, 'layers_per_unit', 1) or 1
                detection = getattr(parent, 'detection_result', None) or {}
                self.suggested_repeat = detection.get('suggested_repeat_units', 2 if detection.get('antiparallel') or detection.get('parent_antiparallel') else 1)
                self.compatible_repeats = detection.get('compatible_repeat_units')
                rises = [track['unit_rise'] for track in detection.get('protofilaments', [])
                         if track.get('unit_rise') is not None]
                self.manual_rise = float(np.median(rises)) if rises else 4.8 * self.unit_layers
                self.options = None
                self.initUI()

            def initUI(self):
                layout = QVBoxLayout(self)
                layout.setContentsMargins(12, 12, 12, 12)
                layout.setSpacing(10)
                self.chk_auto = QCheckBox("Auto Computing")
                self.chk_auto.setChecked(True)
                layout.addWidget(self.chk_auto)

                period_form = QFormLayout()
                self.edit_period = QLineEdit(str(self.suggested_repeat))
                self.edit_period.setMinimumWidth(150)
                self.edit_period.setToolTip("1: A-A-A; 2: A-B-A-B; 3: A-B-C-A-B-C. Zero or any negative whole number: random compatible source conformations.")
                period_form.addRow("Alternating period (chain units):", self.edit_period)
                layout.addLayout(period_form)

                self.chk_compute_axis = QCheckBox("Compute Helical Axis")
                self.chk_compute_axis.setToolTip("Use the source stacking direction. When unchecked, use the global Z axis; align the structure with Z first.")
                self.chk_compute_axis.setChecked(self.current_layers > 1)
                self.chk_compute_axis.stateChanged.connect(self.on_compute_axis_changed)
                layout.addWidget(self.chk_compute_axis)

                form_layout = QFormLayout()
                self.edit_twist = QLineEdit("0.000")
                self.edit_rise = QLineEdit(f"{self.manual_rise:.3f}")
                self.edit_twist.setPlaceholderText("Automatic")
                self.edit_rise.setPlaceholderText("Automatic")
                self.edit_rise.setToolTip("Manual starting value: detected chain-unit rise, or 4.8 Å per layer when no repeat was detected.")
                self._showing_auto = False
                form_layout.addRow("Twist per chain unit (°):", self.edit_twist)
                form_layout.addRow("Rise per chain unit (Å):", self.edit_rise)
                layout.addLayout(form_layout)

                self.parameter_note = QLabel()
                self.parameter_note.setWordWrap(True)
                self.parameter_note.setSizePolicy(QSizePolicy.Policy.Preferred, QSizePolicy.Policy.Minimum)
                layout.addWidget(self.parameter_note)
                self.availability_note = QLabel()
                self.availability_note.setWordWrap(True)
                self.availability_note.setSizePolicy(QSizePolicy.Policy.Preferred, QSizePolicy.Policy.Minimum)
                layout.addWidget(self.availability_note)

                self.btn_expand = QPushButton("Expand")
                self.btn_expand.clicked.connect(self.accept)
                layout.addWidget(self.btn_expand)
                self.setMinimumWidth(max(520, self.fontMetrics().horizontalAdvance("Alternating period (chain units):") + 210))
                self.chk_auto.stateChanged.connect(self.toggle_manual)
                self.edit_period.textChanged.connect(self.toggle_manual)
                self.toggle_manual()

            def alternating_period(self):
                try:
                    return int(self.edit_period.text().strip())
                except ValueError:
                    return None

            def on_compute_axis_changed(self, state):
                if self.chk_compute_axis.isChecked() and self.current_layers <= 1:
                    QMessageBox.warning(self, "Reminder", "Only one chain unit per stack is available; a repeat axis cannot be computed from neighboring units.")

            def toggle_manual(self):
                period = self.alternating_period()
                randomize = period is not None and period <= 0
                effective = self.suggested_repeat if randomize else period
                compatible = self.compatible_repeats is None or effective in self.compatible_repeats
                valid = effective is not None and 0 < effective <= self.current_layers and compatible
                auto_available = valid and self.current_layers > effective
                self.chk_auto.setEnabled(auto_available)
                if not auto_available and self.chk_auto.isChecked():
                    self.chk_auto.blockSignals(True)
                    self.chk_auto.setChecked(False)
                    self.chk_auto.blockSignals(False)
                self.btn_expand.setEnabled(valid)
                is_auto = self.chk_auto.isChecked()
                if is_auto and not self._showing_auto:
                    self._manual_values = (self.edit_twist.text(), self.edit_rise.text())
                    self.edit_twist.clear()
                    self.edit_rise.clear()
                elif not is_auto and self._showing_auto:
                    self.edit_twist.setText(self._manual_values[0])
                    self.edit_rise.setText(self._manual_values[1])
                self._showing_auto = is_auto
                self.edit_twist.setEnabled(not is_auto)
                self.edit_rise.setEnabled(not is_auto)
                layer_word = 'layer' if self.unit_layers == 1 else 'layers'
                if randomize:
                    pattern = 'Random compatible source chain units\nProtein identities and strand directions remain in order.'
                elif period is not None:
                    unit_word = 'unit' if period == 1 else 'units'
                    pattern = f'Repeat every {period} chain {unit_word} ({period * self.unit_layers} layers)'
                else:
                    pattern = 'Enter a whole number; 0 or any negative number selects random.'
                self.parameter_note.setText(
                    f'Current structure:\nOne chain unit spans {self.unit_layers} physical {layer_word}\n\n'
                    f'Alternating pattern:\n{pattern}')
                availability = ''
                if effective is not None and effective > self.current_layers:
                    availability = f'At least {effective} source units are needed to include every template.'
                elif effective is not None and not compatible:
                    availability = 'This period would exchange different proteins or opposite strand orientations. Choose a compatible alternating period.'
                elif not auto_available and valid:
                    availability = f'Automatic mode needs at least {effective + 1} source units. Enter manual parameters for this fragment.'
                self.availability_note.setText(availability)
                self.availability_note.setVisible(bool(availability))
                self.edit_twist.setStyleSheet("color: gray;" if is_auto else "")
                self.edit_rise.setStyleSheet("color: gray;" if is_auto else "")
                self._fit_text()
                QTimer.singleShot(0, self._fit_text)

            def _fit_text(self):
                # Reserve wrapped text height under the current ChimeraX theme.
                text_width = self.minimumWidth() - 24
                for label in (self.parameter_note, self.availability_note):
                    label.setMinimumHeight(max(0, label.heightForWidth(text_width)) if label.text() else 0)
                self.layout().activate()
                height = max(self.sizeHint().height(), self.minimumSizeHint().height())
                self.setMinimumHeight(height)
                self.resize(max(self.width(), self.minimumWidth()), max(self.height(), height))

            def accept(self):
                if not self.btn_expand.isEnabled():
                    return
                automatic = self.chk_auto.isChecked()
                twist = rise = 0.0
                if not automatic:
                    try:
                        twist, rise = float(self.edit_twist.text()), float(self.edit_rise.text())
                        if not math.isfinite(twist) or not math.isfinite(rise) or abs(rise) < 1e-8:
                            raise ValueError()
                    except ValueError:
                        QMessageBox.warning(self, "Invalid Input", "Enter finite twist/rise values, with a nonzero rise.")
                        return
                self.options = dict(use_auto=automatic, repeat_units=self.alternating_period(),
                                    use_computed_axis=self.chk_compute_axis.isChecked(),
                                    manual_twist=twist, manual_rise=rise)
                super().accept()

        def action_trim(self):
            pending = getattr(self, '_expand_dialog', None)
            if pending is not None and pending.isVisible():
                pending.raise_()
                pending.activateWindow()
                return
            self._change_layer_count()

        def _expansion_context(self):
            model = next((m for m in self.session.models.list()
                          if getattr(m, 'id_string', None) == self.working_model_id), None)
            return dict(path=self.working_structure, model=model,
                        units=self.detected_units, layers_per_unit=self.layers_per_unit,
                        stacks={frozenset(s) for s in self.final_sandwiches})

        def _finish_expansion_dialog(self, dialog):
            if getattr(self, '_expand_dialog', None) is dialog:
                self._expand_dialog = None
            dialog.deleteLater()

        def _change_layer_count(self, target_layers=None, expansion_options=None, expansion_context=None):
            if not self.working_structure: return
            if expansion_context is not None:
                current = self._expansion_context()
                if current['path'] != expansion_context['path'] or current['model'] is not expansion_context['model']:
                    QMessageBox.warning(self, "Source changed", "The working model changed while the expansion window was open. Reopen expansion for the current model.")
                    return
            if not self.sync_viewer_to_structure():
                return
            if expansion_context is not None:
                # Refresh live coordinates after scene movement and reject changed topology.
                self.process_pdb(self.working_structure)
                if self._expansion_context() != expansion_context:
                    QMessageBox.warning(self, "Source changed", "The chain units changed while the expansion window was open. Review the current assignment and reopen expansion.")
                    return
            try:
                target_layers = int(self.edit_trim.text()) if target_layers is None else target_layers
            except ValueError:
                QMessageBox.warning(self, "Invalid Input", "Please enter a valid integer for Target Layers.")
                return
                
            total_layers = self.detected_layers
            unit_layers = self.layers_per_unit or 1
            total_units = self.detected_units
            
            if target_layers <= 0:
                QMessageBox.warning(self, "Invalid Input", "Keep at least one layer.")
                return

            if target_layers == total_layers:
                self.text_area.append(f"\nYour model is already having {total_layers} layer(s)")
                return
            if self.detection_result and not self.detection_result['complete']:
                QMessageBox.warning(self, "Incomplete stacks", "The detected stacks have missing units, unequal lengths or different layer multiplicities. Whole-assembly trimming/expansion requires complete matching stacks.")
                return
            if target_layers % unit_layers:
                QMessageBox.warning(self, "Whole chain units required", f"Each chain forms {unit_layers} layers. Choose a multiple of {unit_layers}; this operation keeps chains intact.")
                return
            target_units = target_layers // unit_layers
            
            try:
                edited = StructureBuffer(self.working_structure.name)
                
                if target_layers < total_layers:
                    start_idx = (total_units - target_units) // 2
                    end_idx = start_idx + target_units
                    valid_layer_indices = range(start_idx, end_idx)
                    
                    keep_chains = set()
                    for sandwich in self.final_sandwiches:
                        for i, chain_id in enumerate(sandwich):
                            if i in valid_layer_indices: keep_chains.add(chain_id)
                    
                    self.text_area.append("Calculating water/ion retention based on layer centroids...")
                    xy_lim = self.get_param(self.edit_xy_limit, 3.0)
                    water_z_limit = self.get_param(self.edit_water_z, 4.0)
                    keep_het_lines = self.get_kept_heteroatoms(valid_layer_indices, xy_lim, water_z_limit)
                    
                    if self.working_structure.format == 'cif':
                        self.write_trimmed_cif(self.working_structure, edited, keep_chains, keep_het_lines)
                    else:
                        self.write_trimmed_pdb(self.working_structure, edited, keep_chains, keep_het_lines)
                    self.working_structure.copy_from(edited)
                    self.text_area.append(f"\n>>>> Applied Trimming (Kept middle {target_layers} layers + associated waters)")
                
                else:
                    if self.working_structure.format != 'cif':
                        current_atoms = getattr(self, '_last_live_atom_count', 0)
                        if current_atoms > 0 and total_layers > 0:
                            estimated_atoms = (current_atoms / total_layers) * target_layers
                            estimated_chains = (len(self.results) / total_layers) * target_layers
                            
                            if estimated_atoms > 99999 or estimated_chains > 3906:
                                QMessageBox.warning(
                                    self, 
                                    "PDB Limit Exceeded", 
                                    "The expanded model will exceed the maximum number of atoms (99,999) or chain IDs (at most 2 characters in extended PDB mode) a PDB file can hold. Please convert the current model to CIF format and proceed."
                                )
                                return

                    if expansion_options is None:
                        dialog = self.ExpandDialog(total_units, self)
                        dialog.setWindowTitle(f"Expand to {target_layers} Layers")
                        self._expand_dialog = dialog
                        context = self._expansion_context()
                        dialog.accepted.connect(lambda: self._change_layer_count(
                            target_layers, dialog.options, context))
                        dialog.finished.connect(lambda _result: self._finish_expansion_dialog(dialog))
                        dialog.show()
                        return
                    use_auto = expansion_options['use_auto']
                    repeat_units = expansion_options['repeat_units']
                    use_computed_axis = expansion_options['use_computed_axis']
                    manual_twist = expansion_options['manual_twist']
                    manual_rise = expansion_options['manual_rise']

                    self.text_area.append("Calculating expansion transformations...")
                    layers_to_add = target_units - total_units
                    water_z_limit = self.get_param(self.edit_water_z, 4.0)
                    if self.working_structure.format == 'cif':
                        applied_twist, applied_rise, applied_axis, tilt_deg = self.write_expanded_cif(self.working_structure, edited, layers_to_add, use_auto, manual_twist, manual_rise, water_z_limit, use_computed_axis=use_computed_axis, repeat_units=repeat_units)
                    else:
                        applied_twist, applied_rise, applied_axis, tilt_deg = self.write_expanded_pdb(self.working_structure, edited, layers_to_add, use_auto, manual_twist, manual_rise, water_z_limit, use_computed_axis=use_computed_axis, repeat_units=repeat_units)
                    self.working_structure.copy_from(edited)
                    pattern = 'random' if repeat_units <= 0 else f'{repeat_units} chain units'
                    self.text_area.append(f"\n>>>> Applied Expansion (Added {layers_to_add} chain units / {layers_to_add * unit_layers} layers per protofilament, Alternating period: {pattern})")
                    self.text_area.append("Existing atom coordinates preserved; the shared helical transform was applied only to added copies.")

                if not self.reload_working_model():
                    self.text_area.append("\nThe in-memory edit could not be opened; the previous ChimeraX model was kept.")
                    return
                
                if target_layers > total_layers and getattr(self, 'working_model_id', None):
                    try:
                        from chimerax.core.commands import run
                        run(self.session, f"dssp #{self.working_model_id}")
                        if not self.sync_viewer_to_structure():
                            return
                    except:
                        pass

                self.process_pdb(self.working_structure)
                
                if target_layers > total_layers:
                    self.text_area.append(f"\nApplied Twist: {applied_twist:.5f}°")
                    self.text_area.append(f"Applied Rise: {applied_rise:.5f} Å")
                    self.text_area.append(f"Helical Axis Direction: [{applied_axis[0]:.4f}, {applied_axis[1]:.4f}, {applied_axis[2]:.4f}], with {tilt_deg:.2f}° tilted from the Z-Axis")
                
                if target_layers > total_layers and self.chk_fit_map.isChecked():
                    self.text_area.append("Automatic Fit to Map skipped: expansion preserves existing coordinates and shared helical geometry.")
                if target_layers < total_layers and self.chk_fit_map.isChecked() and getattr(self, 'fit_map_model_id', None) and self.working_model_id:
                    mode = getattr(self, 'fit_map_mode', 'whole')
                    self.text_area.append(f"\n>>>> Fitting ({mode}) to Map #{self.fit_map_model_id}...")
                    from chimerax.core.commands import run
                    fit_count = 0
                    try:
                        if mode == "chain":
                            for sandwich in self.final_sandwiches:
                                for chain_id in sandwich:
                                    run(self.session, f"fitmap #{self.working_model_id}/{chain_id} inMap #{self.fit_map_model_id}")
                                    fit_count += 1
                            self.text_area.append(f"Successfully fitted {fit_count} chains to map.")
                        else:
                            run(self.session, f"fitmap #{self.working_model_id} inMap #{self.fit_map_model_id}")
                            self.text_area.append("Successfully fitted the entire model to map.")
                        
                        if not self.sync_viewer_to_structure():
                            return
                        self.process_pdb(self.working_structure)
                    except Exception as e:
                        self.text_area.append(f"Fit to map warning/error: {e}")
                
            except Exception as e: 
                QMessageBox.critical(self, "Error", f"Action failed: {e}")


        def write_expanded_pdb(self, input_path, output_path, layers_to_add, use_auto, manual_twist, manual_rise, water_z_limit=4.0, use_alt=False, use_computed_axis=True, repeat_units=None):
            editor = StructureEditor(input_path)
            result = expand_structure_layers(editor, self.final_sandwiches, layers_to_add,
                use_auto, manual_twist, manual_rise, water_z_limit, use_alt, use_computed_axis,
                label_generator=self.generate_label, repeat_units=repeat_units)
            report = editor.write(output_path, format='pdb')
            self.text_area.append('Validated expansion: {} atoms.'.format(report['atoms']))
            invalid = [k for k, v in report['metadata'].items() if v['action'] in ('invalidated', 'removed')]
            if invalid:
                self.text_area.append('Source metadata invalidated: ' + ', '.join(invalid))
            return result


        def action_restrain(self):
            if not self.final_sandwiches or self.detected_layers == 0: return
            file_path, _ = QFileDialog.getSaveFileName(self, "Save Restrain File", "restrain_torsion.cxc", "ChimeraX Command (*.cxc);;All Files (*)")
            if file_path:
                try:
                    mapping = self.create_renaming_mapping(self.final_sandwiches)
                    mid_layer_idx = self.detected_units // 2 
                    with open(file_path, 'w') as f:
                        for sandwich in self.final_sandwiches:
                            mid_chain_original = sandwich[mid_layer_idx]
                            mid_id_new = mapping[mid_chain_original]
                            for layer_k in range(len(sandwich)):
                                if layer_k == mid_layer_idx: continue 
                                target_chain_original = sandwich[layer_k]
                                target_id_new = mapping[target_chain_original]
                                f.write(f"isolde restrain torsions /{target_id_new} template /{mid_id_new} angleRange 180\n")
                    self.text_area.append(f"\n[Export] Restrain file generated: {os.path.basename(file_path)}")
                except Exception as e: QMessageBox.critical(self, "Error", f"Failed to save restrain file:\n{str(e)}")

        def action_evaluate_beta(self):
            if not self.working_structure or self.detected_layers <= 1: return
            if not self.working_model_id: return
            
            if not self.sync_viewer_to_structure():
                return
            try:
                self.text_area.append(f"\n>>> Re-evaluating secondary structure using ChimeraX dssp...")
                update_secondary_structure_from_chimerax(self.session, self.working_model_id, self.working_structure)
                self.reload_working_model()
                self.process_pdb(self.working_structure)
                self.text_area.append(">>> β-sheet evaluation completed\n>>> Updated working model")
                    
            except Exception as e:
                self.text_area.append(f"ChimeraX DSSP Error: {e}")

        def check_spatial_and_id_duplicates(self, file_path):
            # Reject duplicate identities rather than leave ambiguous metadata references.
            StructureEditor(file_path)


        def rewrite_split_chains(self, file_path, split_indices):
            raise StructureEditError('Ambiguous chain splitting requires an explicit residue mapping; no automatic renaming was performed.')


        def process_pdb(self, filename, old_sandwiches=None, rename_map=None):
            try:
                editor = StructureEditor(filename)
                chains, waters, ions = editor.preview()
                for issue in editor.report.get('source_issues', []):
                    self.text_area.append('Warning: ' + escape(issue['warning']))
                hint = None
                reference = getattr(self, '_axis_reference', None)
                if reference and reference['path'] == filename:
                    p, q = [], []
                    for cid, rows in chains.items():
                        old = reference['chains'].get(cid, {})
                        for row in rows:
                            if row['id'] in old:
                                p.append(old[row['id']]); q.append(row['coord'])
                    if len(p) >= 6:
                        fit = _fit_detection_core(np.asarray(p), np.asarray(q))
                        if fit is not None and fit['rmsd'] <= 1.0:
                            hint = fit['rotation'] @ np.asarray(reference['axis'])
                result = detect_layers(chains, axis_hint=hint)
                if hint is not None and reference.get('antiparallel', False) and result['orientation'] == 'undetermined':
                    result['parent_antiparallel'] = True
                    result['warnings'].append('The parent stack was antiparallel; this fragment has insufficient neighbors to recheck its orientation.')
                pattern = _expansion_repeat_pattern(chains, result['sandwiches'])
                result['suggested_repeat_units'] = pattern['period']
                result['compatible_repeat_units'] = [p for p in range(1, result['detected_units'] + 1)
                                                     if pattern['compatible'](p)]
                if result.get('parent_antiparallel'):
                    result['suggested_repeat_units'] = max(2, result['suggested_repeat_units'])
                    result['compatible_repeat_units'] = [p for p in result['compatible_repeat_units'] if p % 2 == 0]
                self.detection_result = result
                self.final_sandwiches = result['sandwiches']
                self.detected_units = result['detected_units']
                self.detected_layers = result['detected_layers']
                self.layers_per_unit = result['layers_per_unit']
                if result['axis'] is not None:
                    self._axis_reference = dict(path=filename, axis=result['axis'],
                        antiparallel=result['antiparallel'] or result.get('parent_antiparallel', False),
                        chains={c: {r['id']: r['coord'] for r in rows} for c, rows in chains.items()})
                self.results = {c: {'centroid': self.get_centroid(rows)} for c, rows in chains.items()}
                self.chains_data_plot = {c: [tuple(r['coord']) for r in rows] for c, rows in chains.items()}
                self.current_waters, self.current_ions = waters, ions
                self.lbl_orientation.setText(result['orientation'].capitalize())
                lines = ['--------------- Automatic Detection 3.1 ---------------',
                         result['orientation'].capitalize(),
                         f'Total Chains: {len(chains)}',
                         f'Protofilaments: {len(self.final_sandwiches)}',
                         f'Max Chain Units: {self.detected_units}',
                         f'Max Physical Layers: {self.detected_layers}']
                if result['axis'] is None:
                    lines.append('Stacking axis: undetermined (no inter-unit repeat)')
                else:
                    lines.append('Stacking axis (' + result['axis_source'] + '): [' + ', '.join(f'{v:.4f}' for v in result['axis']) + ']')
                for i, track in enumerate(result['protofilaments']):
                    lines.append(f"\nProtofilament #{i+1} ({track['orientation']}): {track['units']} units × {track['layers_per_unit']} layers/unit = {track['layers']} layers")
                    lines.append(', '.join(track['chains']))
                    if not track['complete']:
                        lines.append(f"Unit positions: {track['unit_positions']}; axial span: {track['axial_layer_span']} layers")
                if old_sandwiches and rename_map:
                    lines.append('\nRenamed chains: ' + ', '.join(f'{c} → {rename_map.get(c, c)}' for s in old_sandwiches for c in s))
                lines.extend('Note: ' + w for w in result['warnings'])
                self.text_area.append('\n'.join(lines))
                self.btn_beta.setEnabled(self.detected_layers > 1)
                self.cmb_protofilaments.blockSignals(True)
                self.cmb_protofilaments.clear()
                for i, pf in enumerate(self.final_sandwiches):
                    self.cmb_protofilaments.addItem(f'Protofilament with Chain {pf[len(pf)//2]}', userData=i)
                self.cmb_protofilaments.blockSignals(False)
                highlights = {s[len(s)//2] for s in self.final_sandwiches}
                self.mol_canvas.plot_chains(self.chains_data_plot, label_ids=highlights, waters=waters, ions=ions)
            except Exception as exc:
                # A failed new analysis must not leave a previous model's assignment active.
                self.detection_result = None
                self.final_sandwiches = []
                self.results = {}
                self.detected_units = self.detected_layers = 0
                self.layers_per_unit = 1
                self._axis_reference = None
                self.lbl_orientation.setText('Undetermined')
                self.cmb_protofilaments.clear()
                self.btn_beta.setEnabled(False)
                self.text_area.append(f'Analysis Error: {exc}')

        def parse_cif_manual_logic(self, filepath):
            return StructureEditor(filepath).preview()


        def parse_pdb_manual_logic(self, filepath):
            return StructureEditor(filepath).preview()


        def get_centroid(self, chain_data):
            if not chain_data: return np.array([0.0, 0.0, 0.0])
            coords = [res['coord'] for res in chain_data]
            return np.mean(coords, axis=0)

        def check_stacking(self, chain_a, chain_b, z_min, z_max, xy_limit):
            res_map_a = {r['id']: r['coord'] for r in chain_a}
            res_map_b = {r['id']: r['coord'] for r in chain_b}
            common_ids = set(res_map_a.keys()).intersection(set(res_map_b.keys()))
            if len(common_ids) < 3: return 0
            
            z_diffs = []
            for rid in common_ids:
                pos_a, pos_b = res_map_a[rid], res_map_b[rid]
                dx, dy, dz = pos_b[0] - pos_a[0], pos_b[1] - pos_a[1], pos_b[2] - pos_a[2]
                if math.sqrt(dx*dx + dy*dy) > xy_limit: return 0
                z_diffs.append(dz)

            avg_dz = np.mean(z_diffs)
            abs_dz = abs(avg_dz)
            
            if z_min <= abs_dz <= z_max:
                if np.std(z_diffs) > 1.5: return 0 
                return 1 if avg_dz > 0 else -1
            return 0

        def get_expansion_order(self, n_layers):
            mid = n_layers // 2; order = [mid]; offset = 1
            while True:
                added = False
                up = mid + offset
                if up < n_layers: order.append(up); added = True
                down = mid - offset
                if down >= 0: order.append(down); added = True
                if not added: break
                offset += 1
            return order

        def generate_label(self, index):
            if self.chk_only_numbers.isChecked(): return str(index)
            use_numbers = self.chk_numbers_first.isChecked()
            if use_numbers:
                if index < 26: return chr(ord('A') + index)
                elif index < 36: return str(index - 26)
                else:
                    temp_index = index - 10; label = ""; temp_index += 1
                    while temp_index > 0:
                        temp_index -= 1; label = chr(ord('A') + (temp_index % 26)) + label; temp_index //= 26
                    return label
            else:
                label = ""; index += 1
                while index > 0: 
                    index -= 1; label = chr(ord('A') + (index % 26)) + label; index //= 26
                return label

        def create_renaming_mapping(self, sandwiches):
            editor = StructureEditor(self.working_structure)
            core = {c for s in sandwiches for c in s}
            reserved = set(editor.chains) - core
            mapping, index = {}, 0
            for layer in self.get_expansion_order(self.detected_units):
                for sandwich in sandwiches:
                    if layer >= len(sandwich):
                        continue
                    while True:
                        label = self.generate_label(index)
                        index += 1
                        if label not in reserved and label not in mapping.values():
                            break
                    mapping[sandwich[layer]] = label
            return mapping

        
        def write_renamed_cif(self, input_filename, output_filename, sandwiches):
            editor = StructureEditor(input_filename)
            editor.rename(self.create_renaming_mapping(sandwiches))
            editor.write(output_filename, format='cif')


        def write_trimmed_cif(self, input_path, output_path, keep_chains, keep_het_lines):
            editor = StructureEditor(input_path)
            # Selections are complete residue identities, never physical text line numbers.
            residues = keep_het_lines or set()
            editor.trim(keep_chains, keep_residues=residues)
            report = editor.write(output_path, format='cif')
            self.text_area.append('Validated trim: {} atoms.'.format(report['atoms']))


        def write_expanded_cif(self, input_path, output_path, layers_to_add, use_auto, manual_twist, manual_rise, water_z_limit=4.0, use_alt=False, use_computed_axis=True, repeat_units=None):
            editor = StructureEditor(input_path)
            result = expand_structure_layers(editor, self.final_sandwiches, layers_to_add,
                use_auto, manual_twist, manual_rise, water_z_limit, use_alt, use_computed_axis,
                label_generator=self.generate_label, repeat_units=repeat_units)
            report = editor.write(output_path, format='cif')
            self.text_area.append('Validated expansion: {} atoms.'.format(report['atoms']))
            invalid = [k for k, v in report['metadata'].items() if v['action'] in ('invalidated', 'removed')]
            if invalid:
                self.text_area.append('Source metadata invalidated: ' + ', '.join(invalid))
            return result


        def write_renamed_pdb(self, input_filename, output_filename, sandwiches):
            editor = StructureEditor(input_filename)
            editor.rename(self.create_renaming_mapping(sandwiches))
            editor.write(output_filename, format='pdb')


        def write_trimmed_pdb(self, input_path, output_path, keep_chains, keep_atom_indices=None):
            editor = StructureEditor(input_path)
            # Selections are complete residue identities, never physical text line numbers.
            residues = keep_atom_indices or set()
            editor.trim(keep_chains, keep_residues=residues)
            report = editor.write(output_path, format='pdb')
            self.text_area.append('Validated trim: {} atoms.'.format(report['atoms']))



    
    class DebugWidget(QWidget):
        def __init__(self, pdb_widget_instance, parent=None):
            super().__init__(parent)
            self.pdb_widget = pdb_widget_instance
            self.initUI()

        def initUI(self):
            layout = QVBoxLayout()
            
            lbl = QLabel("<b>Debug Mode: Inspect Raw Files</b><br><br>"
                         "Save the raw working file <i>before</i> ChimeraX parses and auto-corrects it. "
                         "This allows you to inspect the exact output of your Python scripts and catch formatting/syntax "
                         "errors that ChimeraX silently deletes upon loading.<br><br>"
                         "<i>Note: This fetches the working file from the 'Modifier' tab.</i>")
            lbl.setWordWrap(True)
            layout.addWidget(lbl)
            
            self.btn_save_raw = QPushButton("Export Raw Working File")
            self.btn_save_raw.clicked.connect(self.save_raw_file)
            layout.addWidget(self.btn_save_raw)
            
            layout.addStretch()
            self.setLayout(layout)

        def save_raw_file(self):
            working_file = self.pdb_widget.working_structure
            if working_file is None:
                QMessageBox.warning(self, "No File", "No working file is currently loaded in the Modifier tab.")
                return
            
            ext = ".cif" if working_file.format == 'cif' else ".pdb"
            
            save_path, _ = QFileDialog.getSaveFileName(self, "Save Raw Working File", f"debug_raw{ext}", f"Structure Files (*{ext});;All Files (*)")
            if save_path:
                try:
                    with open(save_path, 'w', encoding='utf-8', newline='\n') as handle:
                        handle.write(working_file.text)
                    QMessageBox.information(self, "Success", f"Raw file successfully exported to:\n{save_path}")
                except Exception as e:
                    QMessageBox.critical(self, "Error", f"Failed to export raw file:\n{e}")

    class LocalResWidget(QWidget):
        def __init__(self, tool_instance, session, parent=None):
            super().__init__(parent)
            self.tool_instance = tool_instance
            self.session = session
            self._last_model_ids = set()
            self.initUI()
            
            self.auto_refresh_timer = QTimer(self)
            self.auto_refresh_timer.timeout.connect(self.populate_models)
            self.auto_refresh_timer.start(1000)

        def initUI(self):
            main_layout = QHBoxLayout(self)
            splitter = QSplitter(Qt.Orientation.Horizontal)
            
            left_widget = QWidget()
            left_layout = QVBoxLayout(left_widget)
            
            grp_models = QGroupBox("Input")
            form_models = QFormLayout()
            
            self.combo_map = QComboBox()
            self.combo_locres = QComboBox()
            self.combo_pdb = QComboBox()
            
            form_models.addRow("Map:", self.combo_map)
            form_models.addRow("LocRes Map:", self.combo_locres)
            form_models.addRow("Model:", self.combo_pdb)
            
            grp_models.setLayout(form_models)
            left_layout.addWidget(grp_models)
            
            grp_params = QGroupBox("Color Threshold")
            form_params = QFormLayout()
            
            self.entry_res = QLineEdit("3.0")
            self.entry_range = QLineEdit("0.2")
            
            form_params.addRow("Resolution (Å):", self.entry_res)
            form_params.addRow("Range:", self.entry_range)
            
            grp_params.setLayout(form_params)
            left_layout.addWidget(grp_params)
            
            grp_preview = QGroupBox("Preview Control")
            v_preview = QVBoxLayout()
            
            self.chk_overlay = QCheckBox("Overlay LocRes on Map")
            self.chk_overlay.setChecked(True)
            v_preview.addWidget(self.chk_overlay)
            
            level_layout = QHBoxLayout()
            self.lbl_level = QLabel("Level: ")
            self.lbl_level.setFixedWidth(120)
            self.slider_level = QSlider(Qt.Orientation.Horizontal)
            self.slider_level.setRange(0, 1000)
            self.slider_level.setEnabled(False)
            level_layout.addWidget(self.lbl_level)
            level_layout.addWidget(self.slider_level)
            v_preview.addLayout(level_layout)
            
            rms_layout = QHBoxLayout()
            self.lbl_rms = QLabel("RMS Level: 5.00")
            self.lbl_rms.setFixedWidth(120)
            self.slider_rms = QSlider(Qt.Orientation.Horizontal)
            self.slider_rms.setRange(0, 400)
            self.slider_rms.setEnabled(False)
            rms_layout.addWidget(self.lbl_rms)
            rms_layout.addWidget(self.slider_rms)
            v_preview.addLayout(rms_layout)
            
            grp_preview.setLayout(v_preview)
            left_layout.addWidget(grp_preview)
            
            left_layout.addStretch()
            
            self.btn_run = QPushButton("Run")
            self.btn_run.setStyleSheet("font-weight: bold; padding: 10px;")
            self.btn_run.clicked.connect(self.execute_commands)
            left_layout.addWidget(self.btn_run)
            
            right_widget = QWidget()
            right_layout = QVBoxLayout(right_widget)
            
            self.fig = Figure(figsize=(3, 3), dpi=100)
            self.fig.patch.set_facecolor('#1e1e1e')
            self.ax = self.fig.add_subplot(111)
            self.ax.set_facecolor('#1e1e1e')
            self.ax.axis('off')
            self.ax.text(0.5, 0.5, "Select Maps to Preview", color='white', 
                         fontsize=12, ha='center', va='center', transform=self.ax.transAxes)
            
            self.canvas = FigureCanvas(self.fig)
            right_layout.addWidget(self.canvas)
            
            splitter.addWidget(left_widget)
            splitter.addWidget(right_widget)
            splitter.setStretchFactor(0, 6)
            splitter.setStretchFactor(1, 4)
            
            main_layout.addWidget(splitter)
            
            self.slider_level.valueChanged.connect(self.on_level_changed)
            self.slider_rms.valueChanged.connect(self.on_rms_changed)
            self.chk_overlay.toggled.connect(self.draw_plot)
            self.entry_res.textChanged.connect(self.draw_plot)
            self.entry_range.textChanged.connect(self.draw_plot)
            self.combo_map.currentIndexChanged.connect(self.draw_plot)
            self.combo_locres.currentIndexChanged.connect(self.draw_plot)
            self.combo_pdb.currentIndexChanged.connect(self.draw_plot)
            
            self.data_map = None
            self.data_locres = None
            self.map_min = 0
            self.map_max = 1
            self.map_rms = 1.0
            self._updating_sliders = False

        def on_level_changed(self, value):
            if self._updating_sliders: return
            self._updating_sliders = True
            
            if hasattr(self, 'map_rms') and self.map_rms > 0:
                slider_frac = value / 1000.0
                current_threshold = self.map_min + slider_frac * (self.map_max - self.map_min)
                current_rms_level = current_threshold / self.map_rms
                
                rms_val_int = int(current_rms_level * 40)
                self.slider_rms.blockSignals(True)
                self.slider_rms.setValue(max(0, min(400, rms_val_int)))
                self.slider_rms.blockSignals(False)
                self.lbl_rms.setText(f"RMS Level: {current_rms_level:.2f}")

            self.draw_plot()
            self._updating_sliders = False

        def on_rms_changed(self, value):
            if self._updating_sliders: return
            self._updating_sliders = True
            
            if hasattr(self, 'map_rms') and self.map_rms > 0:
                snapped_value = round(value / 10) * 10
                if snapped_value != value:
                    self.slider_rms.blockSignals(True)
                    self.slider_rms.setValue(snapped_value)
                    self.slider_rms.blockSignals(False)
                    
                rms_level = snapped_value / 40.0
                current_threshold = rms_level * self.map_rms
                
                if self.map_max > self.map_min:
                    slider_frac = (current_threshold - self.map_min) / (self.map_max - self.map_min)
                    level_val_int = int(slider_frac * 1000)
                    self.slider_level.blockSignals(True)
                    self.slider_level.setValue(max(0, min(1000, level_val_int)))
                    self.slider_level.blockSignals(False)
                    
                self.lbl_rms.setText(f"RMS Level: {rms_level:.2f}")

            self.draw_plot()
            self._updating_sliders = False

        def populate_models(self):
            current_models = self.session.models.list()
            current_ids = {(m.id_string, id(m), m.name) for m in current_models}
            
            if current_ids != self._last_model_ids:
                vol_models = [m for m in current_models if hasattr(m, 'data')]
                atomic_models = [m for m in current_models if hasattr(m, 'atoms')]

                def update_combo(combo, allowed_models):
                    curr_data = combo.currentData()
                    curr_id = curr_data.id_string if curr_data else None
                    combo.blockSignals(True)
                    combo.clear()
                    combo.addItem("--- Select ---", userData=None)
                    for m in allowed_models:
                        combo.addItem(f"#{m.id_string} {m.name}", userData=m)
                    
                    if curr_id:
                        for i in range(combo.count()):
                            m_data = combo.itemData(i)
                            if m_data and m_data.id_string == curr_id:
                                combo.setCurrentIndex(i)
                                break
                    combo.blockSignals(False)

                update_combo(self.combo_map, vol_models)
                update_combo(self.combo_locres, vol_models)
                update_combo(self.combo_pdb, atomic_models)
                
                for i in range(1, self.combo_map.count()):
                    m = self.combo_map.itemData(i)
                    if not m: continue
                    
                    name_lower = getattr(m, 'name', '').lower()
                    locres_keywords = ["locres", "localres", "local_res", "loc_res"]
                    is_locres_map = any(kw in name_lower for kw in locres_keywords)

                    if is_locres_map:
                        curr = self.combo_locres.currentData()
                        if not curr or not any(kw in getattr(curr, 'name', '').lower() for kw in locres_keywords):
                            self.combo_locres.setCurrentIndex(i)
                    else:
                        curr = self.combo_map.currentData()
                        if not curr or any(kw in getattr(curr, 'name', '').lower() for kw in locres_keywords):
                            self.combo_map.setCurrentIndex(i)
                            
                for i in range(1, self.combo_pdb.count()):
                    m = self.combo_pdb.itemData(i)
                    if not m: continue
                    curr = self.combo_pdb.currentData()
                    if not curr:
                        self.combo_pdb.setCurrentIndex(i)

                self._last_model_ids = current_ids
                self.draw_plot()

        def load_data(self):
            map_model = self.combo_map.currentData()
            locres_model = self.combo_locres.currentData()
            
            if not map_model or not locres_model: return False
            
            if getattr(self, '_loaded_map', None) == map_model and getattr(self, '_loaded_locres', None) == locres_model:
                return True
            
            try:
                self.data_map = map_model.data.matrix().astype(np.float32)
                self.data_locres = locres_model.data.matrix().astype(np.float32)
                
                self.map_min = float(np.min(self.data_map))
                self.map_max = float(np.max(self.data_map))
                if self.map_max == self.map_min:
                    self.map_max = self.map_min + 1.0
                    
                self.map_rms = float(np.sqrt(np.mean(np.square(self.data_map))))
                default_rms_level = 5.0
                default_thresh = float(default_rms_level * self.map_rms)
                
                default_slider_val = int(1000 * (default_thresh - self.map_min) / (self.map_max - self.map_min))
                default_slider_val = max(0, min(1000, default_slider_val))
                
                self.slider_level.blockSignals(True)
                self.slider_level.setValue(default_slider_val)
                self.slider_level.setEnabled(True)
                self.slider_level.blockSignals(False)
                
                self.slider_rms.blockSignals(True)
                self.slider_rms.setValue(int(default_rms_level * 40))
                self.slider_rms.setEnabled(True)
                self.slider_rms.blockSignals(False)
                self.lbl_rms.setText(f"RMS Level: {default_rms_level:.2f}")
                
                if self.data_map.shape == self.data_locres.shape:
                    mask = self.data_map > default_thresh
                    valid_locres = self.data_locres[mask]
                    if len(valid_locres) > 0:
                        optimal_res = round(float(np.percentile(valid_locres, 30)), 1)
                        self.entry_res.blockSignals(True)
                        self.entry_res.setText(f"{optimal_res:.1f}")
                        self.entry_res.blockSignals(False)
                
                self._loaded_map = map_model
                self._loaded_locres = locres_model
                return True
            except Exception:
                return False

        def draw_plot(self):
            if not self.load_data() or self.data_map is None or self.data_locres is None:
                self.ax.clear()
                self.ax.axis('off')
                self.ax.text(0.5, 0.5, "Select Maps to Preview", color='white', 
                             fontsize=12, ha='center', va='center', transform=self.ax.transAxes)
                self.canvas.draw_idle()
                
                self.slider_level.blockSignals(True)
                self.slider_level.setValue(0)
                self.slider_level.setEnabled(False)
                self.slider_level.blockSignals(False)
                self.lbl_level.setText("Level: ")
                
                self.slider_rms.blockSignals(True)
                self.slider_rms.setValue(200)
                self.slider_rms.setEnabled(False)
                self.slider_rms.blockSignals(False)
                self.lbl_rms.setText("RMS Level: 5.00")
                
                return
                
            import matplotlib.colors as mcolors
            
            try:
                r = float(self.entry_res.text().strip())
                rg = float(self.entry_range.text().strip())
            except ValueError:
                return

            max_z = self.data_map.shape[0]
            current_z = max_z // 2
            slice_thickness = 5
            
            slider_frac = self.slider_level.value() / 1000.0
            current_threshold = self.map_min + slider_frac * (self.map_max - self.map_min)
            self.lbl_level.setText(f"Level: {current_threshold:.4f}")

            z_start = max(0, current_z - slice_thickness // 2)
            z_end = min(max_z, current_z + slice_thickness // 2 + 1)

            slice_map = np.mean(self.data_map[z_start:z_end, :, :], axis=0)
            
            if self.data_map.shape == self.data_locres.shape:
                slice_locres = np.mean(self.data_locres[z_start:z_end, :, :], axis=0)
            else:
                self.ax.clear()
                self.ax.axis('off')
                self.ax.text(0.5, 0.5, "Maps do not match!", color='red', 
                             fontsize=12, ha='center', va='center', transform=self.ax.transAxes)
                self.canvas.draw_idle()
                return

            self.ax.clear()
            self.ax.axis('off')

            import matplotlib.patches as patches
            
            center_y, center_x = slice_map.shape[0] / 2, slice_map.shape[1] / 2
            radius = min(slice_map.shape) * 0.3 
            
            clip_circle = patches.Circle((center_x, center_y), radius, transform=self.ax.transData)
            
            img_base = self.ax.imshow(slice_map, cmap='gray', origin='lower')
            img_base.set_clip_path(clip_circle)

            if self.chk_overlay.isChecked():
                cmap = mcolors.LinearSegmentedColormap.from_list("rwb", ["red", "white", "blue"])
                vmin = r - rg
                vmax = r + rg
                norm = mcolors.Normalize(vmin=vmin, vmax=vmax)

                import scipy.ndimage
                base_mask = np.where(slice_map > current_threshold, 1.0, 0.0)
                smooth_mask = scipy.ndimage.gaussian_filter(base_mask, sigma=0.8)
                alpha_mask = np.clip(smooth_mask * 1.5, 0.0, 0.9) 
                
                rgba_img = cmap(norm(slice_locres))
                rgba_img[..., 3] = alpha_mask 

                img_overlay = self.ax.imshow(rgba_img, origin='lower', interpolation='bilinear')
                img_overlay.set_clip_path(clip_circle)

            self.ax.set_xlim(center_x - radius * 1.05, center_x + radius * 1.05)
            self.ax.set_ylim(center_y - radius * 1.05, center_y + radius * 1.05)

            self.canvas.draw_idle()

        def execute_commands(self):
            map_model = self.combo_map.currentData()
            locres_model = self.combo_locres.currentData()
            pdb_model = self.combo_pdb.currentData()

            if not map_model or not locres_model or not pdb_model:
                QMessageBox.warning(self, "Missing Selection", "Please ensure all inputs are provided")
                return

            map_id = map_model.id_string
            locres_id = locres_model.id_string
            pdb_id = pdb_model.id_string

            try:
                res_val = float(self.entry_res.text().strip())
                range_val = float(self.entry_range.text().strip())
            except ValueError:
                QMessageBox.critical(self, "Error", "Resolution and Range must be numbers")
                return

            low = res_val - range_val
            mid = res_val
            high = res_val + range_val

            commands = f"""
color sample #{map_id} map #{locres_id} palette {low:.1f},#ff0000:{mid:.1f},#ffffff:{high:.1f},#0000ff 
surface dust #{map_id}
set bgColor #ffffff00
hide #{locres_id} models 
graphics silhouettes true
graphics silhouette width 20
ui tool show "Surface Color"
key red-white-blue :{low:.1f}Å :{mid:.1f}Å :{high:.1f}Å
key fontSize 14
key pos 0.8000,0.06000
key size 0.15000,0.03000
key borderWidth 3.0
surface zone #{map_id} near #{pdb_id}
hide #{pdb_id} models 
hide #{locres_id} models
view orient
zoom pixelSize 0.18
scalebar 10 xpos 0.01 ypos 0.01
2dlabels create scalebar_legend text "10 Å" xpos 0.01 ypos 0.02 size 20
volume #{map_id} step 1
volume #{map_id} rmsLevel {self.slider_rms.value() / 40.0:.2f}
lighting soft
"""
            from chimerax.core.commands import run
            for cmd in commands.strip().split('\n'):
                if cmd.strip():
                    try:
                        run(self.session, cmd.strip())
                    except Exception as e:
                        print(f"Command failed: {cmd}\n{e}")

    class ModifierToolTabs(QTabWidget):
        def __init__(self, tool_instance, session):
            super().__init__()
            self.pdb_widget = PDBLayerIdentifier(tool_instance, session)
            self.cif_widget = CIFLayerIdentifier(tool_instance, session)
            self.localres_widget = LocalResWidget(tool_instance, session)
            
            self.addTab(self.pdb_widget, "Modifier")
            self.addTab(self.cif_widget, "Layer Viewer")
            self.addTab(self.localres_widget, "LocalRes")
            
            max_min_height = max(
                self.pdb_widget.minimumSizeHint().height(),
                self.cif_widget.minimumSizeHint().height(),
                self.localres_widget.minimumSizeHint().height()
            )
            self.setMinimumHeight(max_min_height)
            

    return ModifierToolTabs


class PDBModifierTool(ToolInstance):
    SESSION_ENDURING = False

    @property
    def tool_info(self):
        # Standalone launches have no bundle, but ChimeraX still queries tool_info.
        bundle = self.bundle_info
        return next((t for t in bundle.tools if t.name == self.tool_name), None) if bundle is not None else None
    
    def __init__(self, session, tool_name):
        super().__init__(session, tool_name)
        self.display_name = "Amyloid Modifier"
        
        from chimerax.ui import MainToolWindow
        self.tool_window = MainToolWindow(self)
        
        WidgetClass = open_chain_modifier()
        self.widget = WidgetClass(self, session)
        
        from PyQt6.QtWidgets import QVBoxLayout
        layout = QVBoxLayout()
        layout.setContentsMargins(0, 0, 0, 0)
        layout.addWidget(self.widget)
        self.tool_window.ui_area.setLayout(layout)

        self.tool_window.manage("right")
        
        self.widget.setStyleSheet("""
            QMainWindow, QDialog, QWidget { background-color: #1e1e1e; color: #d4d4d4; font-family: Arial; font-size: 8pt; }
            
            QPushButton { background-color: #3e3e42; color: #d4d4d4; border: 1px solid #3e3e42; padding: 5px 15px; border-radius: 4px; }
            QPushButton:hover { background-color: #4e4e52; border: 1px solid #98c379; }
            QPushButton:pressed, QPushButton:checked { background-color: #98c379; color: #1e1e1e; border: 1px solid #98c379; }
            
            QLineEdit, QSpinBox, QDoubleSpinBox, QTextEdit, QPlainTextEdit, QTextBrowser, QTableWidget, QComboBox { 
                background-color: #1e1e1e; border: 1px solid #3e3e42; color: #cccccc; padding: 2px; border-radius: 2px;
                selection-background-color: #98c379; selection-color: #1e1e1e;
            }
            
            QSpinBox::up-button, QDoubleSpinBox::up-button,
            QSpinBox::down-button, QDoubleSpinBox::down-button {
                width: 20px; 
            }

            QComboBox::drop-down { border: none; }
            QComboBox QAbstractItemView { 
                background-color: #252526; color: #d4d4d4; border: 1px solid #3e3e42; outline: none;
                selection-background-color: #98c379; selection-color: #1e1e1e;
            }
            QComboBox QAbstractItemView::item { padding: 2px 5px; min-height: 16px; border: none !important; }
            QComboBox QAbstractItemView::item:hover,
            QComboBox QAbstractItemView::item:selected { 
                background-color: #98c379; color: #1e1e1e; border: none !important; outline: none !important;
            }

            QTreeWidget, QListWidget, QTreeView { 
                background-color: #252526; border: 1px solid #3e3e42; color: #cccccc; outline: none; 
            }
            QTreeWidget::item, QListWidget::item { padding: 5px; }
            QTreeWidget::item:selected, QListWidget::item:selected { 
                background-color: #3e3e42; color: #ffffff; border-left: 3px solid #98c379; 
            }

            QHeaderView::section { background-color: #252526; color: #98c379; border: 1px solid #3e3e42; padding: 4px; font-weight: bold; }
            QTableCornerButton::section { background-color: #252526; border: 1px solid #3e3e42; }

            QMenu { background-color: #252526; color: #d4d4d4; border: 1px solid #3e3e42; }
            QMenu::item:selected { background-color: #3e3e42; color: #98c379; }

            QProgressBar { 
                border: 1px solid #3e3e42; border-radius: 2px; 
                background-color: #1e1e1e; text-align: center; color: #d4d4d4; 
            }
            QProgressBar::chunk { background-color: #98c379; border-radius: 2px; color: #1e1e1e; }

            QGroupBox { border: 1px solid #3e3e42; border-radius: 4px; margin-top: 1.0em; font-weight: bold; color: #98c379; }
            QGroupBox::title { subcontrol-origin: margin; subcontrol-position: top left; padding: 0 5px; }
            
            QCheckBox::indicator { width: 14px; height: 14px; border: 1px solid #555; border-radius: 2px; background-color: #1e1e1e; }
            QCheckBox::indicator:checked { background-color: #98c379; border: 1px solid #98c379; }
            QRadioButton::indicator { width: 12px; height: 12px; border: 1px solid #555; border-radius: 7px; background-color: #1e1e1e; }
            QRadioButton::indicator:checked { background-color: #98c379; border: 1px solid #98c379; }
            
            QSlider:vertical { min-width: 20px; }
            QSlider::groove:vertical { background: #3c3c3c; width: 6px; border-radius: 3px; }
            QSlider::handle:vertical { background: #98c379; height: 14px; width: 14px; margin: 0 -4px; border-radius: 7px; }
            QSlider::handle:vertical:hover { background: #b5e890; }
            
            QSlider:horizontal { min-height: 20px; }
            QSlider::groove:horizontal { background: #3c3c3c; height: 6px; border-radius: 3px; }
            QSlider::handle:horizontal { background: #98c379; height: 14px; width: 14px; margin: -4px 0; border-radius: 7px; }
            QSlider::handle:horizontal:hover { background: #b5e890; }
            
            QTabWidget::tab-bar { alignment: left; }
            QTabWidget::pane { border: 1px solid #3e3e42; top: -1px; }
            QTabBar::tab { background: #252526; border: 1px solid #3e3e42; padding: 6px 12px; color: #d4d4d4; }
            QTabBar::tab:selected { background: #3e3e42; color: #98c379; font-weight: bold; border-bottom: 1px solid #3e3e42; }
            QTabBar::tab:hover { background: #4e4e52; }
            
            QScrollBar:vertical { background: transparent; width: 12px; margin: 0px; }
            QScrollBar::handle:vertical { background-color: #4e4e52; min-height: 20px; border-radius: 5px; margin: 2px; }
            QScrollBar::handle:vertical:hover, QScrollBar::handle:vertical:pressed { background-color: #98c379; }
            QScrollBar::add-line:vertical, QScrollBar::sub-line:vertical { height: 0px; background: transparent; }
            QScrollBar::add-page:vertical, QScrollBar::sub-page:vertical { background: transparent; }
            
            QScrollBar:horizontal { background: transparent; height: 12px; margin: 0px; }
            QScrollBar::handle:horizontal { background-color: #4e4e52; min-width: 20px; border-radius: 5px; margin: 2px; }
            QScrollBar::handle:horizontal:hover, QScrollBar::handle:horizontal:pressed { background-color: #98c379; }
            QScrollBar::add-line:horizontal, QScrollBar::sub-line:horizontal { width: 0px; background: transparent; }
            QScrollBar::add-page:horizontal, QScrollBar::sub-page:horizontal { background: transparent; }
        """)
        
    def delete(self):
        for sub_widget in (self.widget.pdb_widget, self.widget.cif_widget):
            sub_widget.auto_refresh_timer.stop()
            sub_widget.working_structure = None
        super().delete()
