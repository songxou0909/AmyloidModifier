import os
import math
import shutil
import atexit
import uuid

from PyQt6.QtWidgets import (QApplication, QWidget, QLabel, QLineEdit, 
                             QPushButton, QTextEdit, QMessageBox, QFileDialog, 
                             QVBoxLayout, QHBoxLayout, QSizePolicy, QSlider,
                             QSplitter, QSpinBox, QGroupBox, QFormLayout, 
                             QCheckBox, QTableWidget, QHeaderView, QTableWidgetItem, 
                             QTextBrowser, QDialog, QComboBox, QTabWidget)
from PyQt6.QtCore import Qt, QTimer
from chimerax.core.tools import ToolInstance
from chimerax.core.commands import run

toplevel_windows = []

def show_error_message(message):
    msg = QMessageBox()
    msg.setIcon(QMessageBox.Icon.Critical)
    msg.setWindowTitle("Error")
    msg.setText(message)
    msg.exec()

# ====================== Multi-layer Chain Modifier Class ======================
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

    # ====================== CIF Logic Modules ======================
    def generate_chain_id(index):
        chars = "ABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789"
        if index < len(chars): return chars[index]
        n = index - len(chars); base = len(chars)
        return chars[n // base] + chars[n % base]

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
            layout.addWidget(QLabel("Select the layers you want to KEEP\nChains will be automatically renamed to maximum 2 letters upon saving."))
            self.table = QTableWidget()
            self.table.setColumnCount(4)
            self.table.setHorizontalHeaderLabels(["Keep", "Rank", "Original Chains", "Residues"])
            self.table.horizontalHeader().setSectionResizeMode(QHeaderView.ResizeMode.Stretch)
            self.table.horizontalHeader().setSectionResizeMode(0, QHeaderView.ResizeMode.ResizeToContents)
            self.table.horizontalHeader().setSectionResizeMode(1, QHeaderView.ResizeMode.ResizeToContents)
            self.populate_table()
            layout.addWidget(self.table)
            btn_layout = QHBoxLayout(); btn_layout.addStretch()
            self.btn_save = QPushButton("Save New CIF")
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

                rank_item = QTableWidgetItem(str(i + 1))
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
            import uuid
            import os
            from chimerax.core.commands import run

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
            name_part, ext = os.path.splitext(self.original_filename)
            
            folder = os.path.join(os.path.expanduser("~"), ".chimerax_modifier_temp")
            os.makedirs(folder, exist_ok=True)
            unique_suffix = uuid.uuid4().hex[:8]
            path = os.path.join(folder, f"{name_part}_trimmed_{unique_suffix}.cif")
            
            try:
                self.write_cif(path, mapping)
                
                parent_widget = self.parent()
                if parent_widget:
                    parent_widget.created_temp_files.append(path)
                    run(parent_widget.session, f'open "{path}"')
                    
                QMessageBox.information(self, "Success", f"Trimmed model loaded into ChimeraX!\n(Renamed {len(mapping)} chains)")
                self.accept()
            except Exception as e: 
                QMessageBox.critical(self, "Error", str(e))

        def write_cif(self, path, mapping):
            data = self.full_cif_data
            cm = data['col_map']
            idx_auth = cm.get('auth_asym_id')
            idx_label = cm.get('label_asym_id')
            target_idx = idx_auth if idx_auth is not None else idx_label
            if target_idx is None: raise ValueError("No chain ID columns found.")
            with open(path, 'w') as f:
                f.writelines(data['pre_loop'])
                f.writelines(data['loop_headers'])
                for old_id, lines in data['chain_lines'].items():
                    if old_id in mapping:
                        new_id = mapping[old_id]
                        for line in lines:
                            parts = line.strip().split()
                            if len(parts) > target_idx:
                                if idx_auth is not None: parts[idx_auth] = new_id
                                if idx_label is not None: parts[idx_label] = new_id
                                f.write(" ".join(parts) + "\n")
                            else: f.write(line)

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
            self.created_temp_files = [] 
            self._last_model_ids = set()
            
            def cleanup_files():
                for f in self.created_temp_files:
                    if os.path.exists(f):
                        try: os.remove(f)
                        except: pass
            atexit.register(cleanup_files)
            
            self.working_file_path = None; self.layers_result = []; self.full_cif_data = {}
            self.initUI()
            qt_app = QApplication.instance()
            if qt_app: 
                self.setStyleSheet(qt_app.styleSheet())

        def initUI(self):
            main_l = QHBoxLayout()
            splitter = QSplitter(Qt.Orientation.Horizontal)
            left = QWidget(); l_lay = QVBoxLayout(left); l_lay.setSpacing(10); l_lay.setContentsMargins(10,10,10,10)
            
            hb_load = QHBoxLayout()
            self.cmb_models = QComboBox()
            self.populate_models()
            hb_load.addWidget(self.cmb_models)
            
            self.auto_refresh_timer = QTimer(self)
            self.auto_refresh_timer.timeout.connect(self.populate_models)
            self.auto_refresh_timer.start(1000)
            
            self.btn_load_model = QPushButton('Load and Copy from ChimeraX')
            self.btn_load_model.clicked.connect(self.load_from_chimerax)
            hb_load.addWidget(self.btn_load_model)
            l_lay.addLayout(hb_load)
            
            self.lbl_status = QLabel("No file loaded"); self.lbl_status.setStyleSheet("color: gray; font-style: italic;")
            l_lay.addWidget(self.lbl_status)
            
            self.grp_ops = QGroupBox("Operations"); self.grp_ops.setEnabled(False)
            ops_l = QVBoxLayout()
            self.btn_trim = QPushButton("Trim Best Layers")
            self.btn_trim.setToolTip("Select the layers wanted and save as a new CIF file")
            self.btn_trim.clicked.connect(self.open_trim_dialog)
            ops_l.addWidget(self.btn_trim)
            self.grp_ops.setLayout(ops_l); l_lay.addWidget(self.grp_ops)
            
            grp_p = QGroupBox("Detection Parameters"); form = QFormLayout()
            self.ed_chn = QLineEdit("4"); self.ed_chn.setValidator(QIntValidator())
            form.addRow("Chains per Layer:", self.ed_chn)
            self.ed_zmn = QLineEdit("0.0"); self.ed_zmn.setValidator(QDoubleValidator())
            form.addRow("Z Shift Min (Å):", self.ed_zmn)
            self.ed_zmx = QLineEdit("4.0"); self.ed_zmx.setValidator(QDoubleValidator())
            form.addRow("Z Shift Max (Å):", self.ed_zmx)
            self.ed_chn.editingFinished.connect(self.recalc); self.ed_zmn.editingFinished.connect(self.recalc); self.ed_zmx.editingFinished.connect(self.recalc)
            grp_p.setLayout(form); l_lay.addWidget(grp_p)
            
            self.txt = QTextBrowser(); self.txt.setOpenLinks(False); self.txt.anchorClicked.connect(self.link_clk)
            self.txt.setStyleSheet("font-family: Consolas, monospace;")
            l_lay.addWidget(self.txt)
            
            right = QWidget(); r_lay = QHBoxLayout(right); r_lay.setContentsMargins(0,0,0,0)
            self.cvs = CIFMol3DCanvas(self, dpi=100); r_lay.addWidget(self.cvs, 1)
            sl = QSlider(Qt.Orientation.Vertical); sl.setRange(1, 100); sl.setValue(50); sl.valueChanged.connect(lambda v: self.cvs.set_zoom_from_slider(v))
            self.cvs.slider_callback = lambda v: sl.setValue(v)
            r_lay.addWidget(sl)
            
            splitter.addWidget(left); splitter.addWidget(right)
            splitter.setSizes([250, 300]) 
            splitter.setStretchFactor(0, 4); splitter.setStretchFactor(1, 6)
            main_l.addWidget(splitter); self.setLayout(main_l)

        def populate_models(self):
            from chimerax.atomic import AtomicStructure
            current_models = self.session.models.list(type=AtomicStructure)
            current_ids = {m.id_string for m in current_models}
            if current_ids != self._last_model_ids:
                current_sel = self.cmb_models.currentData()
                current_sel_id = current_sel.id_string if current_sel else None
                self.cmb_models.blockSignals(True)
                self.cmb_models.clear()
                for model in current_models:
                    self.cmb_models.addItem(f"#{model.id_string} {model.name}", userData=model)
                if current_sel_id:
                    for i in range(self.cmb_models.count()):
                        m = self.cmb_models.itemData(i)
                        if m and m.id_string == current_sel_id:
                            self.cmb_models.setCurrentIndex(i)
                            break
                self.cmb_models.blockSignals(False)
                self._last_model_ids = current_ids

        def load_from_chimerax(self):
            model = self.cmb_models.currentData()
            if not model: return
            try:
                self.loaded_filename = model.name
                self.lbl_status.setText(f"Loaded: {self.loaded_filename}")
                
                folder = os.path.join(os.path.expanduser("~"), ".chimerax_modifier_temp")
                os.makedirs(folder, exist_ok=True)
                unique_suffix = uuid.uuid4().hex[:8]
                temp_path = os.path.join(folder, f"working_copy_{unique_suffix}.cif")
                
                run(self.session, f'save "{temp_path}" models #{model.id_string} format mmcif')
                self.created_temp_files.append(temp_path)
                
                self.working_file_path = temp_path
                self.grp_ops.setEnabled(True)
                self.txt.clear(); self.txt.append(f"Loaded {model.name}")
                self.process_cif(temp_path)
            except Exception as e:
                self.txt.append(f"Load Error: {e}")

        def recalc(self):
            if self.working_file_path: QTimer.singleShot(10, lambda: self.process_cif(self.working_file_path))

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
            chains_c = {}; waters = []; ions = []; full_d = {'pre_loop':[], 'loop_headers':[], 'chain_lines':{}, 'col_map':{}}
            with open(path,'r') as f: lines=f.readlines()
            in_h = False; buf = []; d_idx = -1; headers = []
            for i, l in enumerate(lines):
                s = l.strip()
                if d_idx == -1 and not in_h:
                    if s=="loop_": in_h=True; buf=[l]
                    else: full_d['pre_loop'].append(l)
                    continue
                if in_h:
                    if s.startswith("_"): buf.append(l)
                    else:
                        if any("_atom_site." in x for x in buf):
                            headers = [x.strip() for x in buf if x.strip().startswith("_")]
                            full_d['loop_headers']=buf; d_idx=i; break
                        else:
                            full_d['pre_loop'].extend(buf); in_h=False
                            if s=="loop_": in_h=True; buf=[l]
                            else: full_d['pre_loop'].append(l)
            if not headers: return {}, [], [], {}
            cm = {h.split('.')[1] if '.' in h else h: i for i, h in enumerate(headers)}
            for k, v in [("auth_asym_id", "auth_asym_id"), ("label_asym_id", "label_asym_id"), ("label_atom_id", "label_atom_id")]:
                for h in headers: 
                    if k in h: cm[k] = headers.index(h)
            
            full_d['col_map'] = cm
            ic = cm.get("auth_asym_id", cm.get("label_asym_id"))
            ia = cm.get("label_atom_id")
            ix, iy, iz = cm.get("Cartn_x"), cm.get("Cartn_y"), cm.get("Cartn_z")
            ig = cm.get("group_PDB")
            ir = cm.get("label_comp_id", cm.get("auth_comp_id"))
            
            for i in range(d_idx, len(lines)):
                ln = lines[i]; s = ln.strip()
                if not s or s.startswith(("#", "_", "loop_")): continue
                p = s.split()
                if len(p) < len(headers): continue
                cid = p[ic]
                if cid not in full_d['chain_lines']: full_d['chain_lines'][cid] = []
                full_d['chain_lines'][cid].append(ln)
                
                try: coord = np.array([float(p[ix]), float(p[iy]), float(p[iz])])
                except: continue
                
                if ia is not None and p[ia] == "CA":
                    if cid not in chains_c: chains_c[cid]=[]
                    chains_c[cid].append(coord)
                elif ig is not None and p[ig] == "HETATM":
                    res_name = p[ir] if ir is not None else ""
                    atom_name = p[ia] if ia is not None else ""
                    if res_name in ["HOH", "WAT"]:
                        if atom_name.startswith("O"): waters.append(tuple(coord))
                    else:
                        ions.append(tuple(coord))
            return chains_c, waters, ions, full_d

        def process_cif(self, path):
            chn_max = int(self.ed_chn.text()); zmn = float(self.ed_zmn.text()); zmx = float(self.ed_zmx.text())
            try:
                self.txt.append("Parsing..."); coords, waters, ions, self.full_cif_data = self.parse_cif(path)
                if not coords: self.txt.append("No CA atoms found."); return
                self.chains_data_plot = {c: [tuple(p) for p in coords[c]] for c in coords}
                self.current_waters = waters
                self.current_ions = ions
                cents = {c: np.mean(coords[c], axis=0) for c in coords}
                pool = sorted(list(coords.keys()), key=lambda c: cents[c][2])
                self.layers_result = []
                while pool:
                    ref = pool.pop(0); cur = [ref]; rz = cents[ref][2]
                    cands = sorted([(abs(cents[c][2]-rz), c) for c in pool if zmn <= abs(cents[c][2]-rz) <= zmx], key=lambda x:x[0])
                    sel = [c[1] for c in cands[:chn_max-1]]
                    cur.extend(sel); 
                    for c in sel: pool.remove(c)
                    res_cnt = sum(len(coords[c]) for c in cur)
                    self.layers_result.append({'chains': cur, 'count': len(cur), 'residues': res_cnt})
                self.layers_result.sort(key=lambda x: (x['count'], x['residues']), reverse=True)
                
                h_lines = []
                high_ids = set()
                for i, lay in enumerate(self.layers_result):
                    lnk = f"<a href='{i}' style='color: #97c379; font-weight: bold; text-decoration: none;'>Layer {i+1}</a>"
                    h_lines.append(f"{lnk}: {', '.join(lay['chains'])} [Chains, Residues]:{lay['count']}, {lay['residues']}<br>")
                    if i==0: high_ids.update(lay['chains'])
                self.txt.setHtml("<br>".join(h_lines))
                self.cvs.plot_chains(self.chains_data_plot, label_ids=high_ids, waters=waters, ions=ions)
                
                if hasattr(self, 'trim_dialog') and self.trim_dialog and self.trim_dialog.isVisible():
                    self.trim_dialog.update_data(self.layers_result, self.full_cif_data, self.loaded_filename)

            except Exception as e: self.txt.append(str(e))

    # ====================== Multi-layer Chain Modifier Class (PDB) ======================
    class PDBLayerIdentifier(QWidget):
        
        CHAIN_RECORD_SPECS = {
            "DBREF":  [(11, 13)], "DBREF1": [(11, 13)], "DBREF2": [(11, 13)],
            "SEQADV": [(15, 17)], "MODRES": [(15, 17)], "SEQRES": [(10, 12)],
            "HET":    [(11, 13)], "SSBOND": [(14, 16), (28, 30)], "CISPEP": [(14, 16), (28, 30)],
            "LINK":   [(20, 22), (50, 52)], "SITE":   [(21, 23), (32, 34), (43, 45), (53, 56)],
            "ATOM":   [(20, 22)], "ANISOU": [(20, 22)], "TER":    [(20, 22)],
            "HETATM": [(20, 22)], "HELIX":  [(18, 20), (30, 32)], "SHEET":  [(20, 22), (31, 33)]
        }

        def __init__(self, tool_instance, session):
            super().__init__()
            self.tool_instance = tool_instance
            self.session = session
            self.working_model_id = None
            self.final_sandwiches = [] 
            self.detected_layers = 0 
            self.created_temp_files = [] 
            self._last_model_ids = set() 
            
            def cleanup_files():
                for f in self.created_temp_files:
                    if os.path.exists(f):
                        try: os.remove(f)
                        except: pass
            atexit.register(cleanup_files)

            self.working_file_path = None 
            self.original_filename_display = ""
            self.initUI()

        def initUI(self):
            main_layout = QHBoxLayout()
            left_widget = QWidget()
            layout = QVBoxLayout(left_widget)
            layout.setSpacing(10)
            layout.setContentsMargins(10, 10, 10, 10)

            hb_load = QHBoxLayout()
            self.cmb_models = QComboBox()
            self.populate_models()
            hb_load.addWidget(self.cmb_models)
            
            self.auto_refresh_timer = QTimer(self)
            self.auto_refresh_timer.timeout.connect(self.populate_models)
            self.auto_refresh_timer.start(1000)
            
            self.btn_load_model = QPushButton('Load and Copy from ChimeraX')
            self.btn_load_model.clicked.connect(self.load_from_chimerax)
            hb_load.addWidget(self.btn_load_model)
            hb_load.addStretch()

            self.chk_anti = QCheckBox("Anti")
            self.chk_anti.setToolTip("Check this if loading an anti-parallel amyloid")
            
            def toggle_anti_mode(state):
                self.edit_z_min.blockSignals(True)
                self.edit_z_max.blockSignals(True)
                self.edit_xy_limit.blockSignals(True)
                
                if state:
                    self.edit_z_min.setText("9.0") 
                    self.edit_z_max.setText("11.5")
                    self.edit_xy_limit.setText("6.0")
                else:
                    self.edit_z_min.setText("4.6")
                    self.edit_z_max.setText("5.0")
                    self.edit_xy_limit.setText("3.0")
                    
                self.edit_z_min.blockSignals(False)
                self.edit_z_max.blockSignals(False)
                self.edit_xy_limit.blockSignals(False)
                self.on_param_change()
                
            self.chk_anti.stateChanged.connect(toggle_anti_mode)
            hb_load.addWidget(self.chk_anti)
            layout.addLayout(hb_load)

            self.grp_ops = QGroupBox("Operations (Applied to Temporary Copy)")
            self.grp_ops.setEnabled(False) 
            v_ops = QVBoxLayout()

            hb_trim = QHBoxLayout()
            self.btn_trim = QPushButton("Trim/Expand Layers")
            self.btn_trim.setToolTip("Trim or expand the structure to the target number of layers")
            self.btn_trim.clicked.connect(self.action_trim)
            hb_trim.addWidget(self.btn_trim)
            hb_trim.addWidget(QLabel("Target Layers:"))
            self.spin_trim = QSpinBox()
            self.spin_trim.setRange(1, 500)
            self.spin_trim.setValue(5)
            hb_trim.addWidget(self.spin_trim)
            v_ops.addLayout(hb_trim)

            hb_ren = QHBoxLayout()
            self.btn_rename = QPushButton("Rename Chains")
            self.btn_rename.clicked.connect(self.action_rename)
            hb_ren.addWidget(self.btn_rename)
            hb_ren.addWidget(QLabel("(Rename Chain ID from the Middle Out)"))
            v_ops.addLayout(hb_ren)

            self.btn_beta = QPushButton("Re-evaluate β-sheet structure")
            self.btn_beta.clicked.connect(self.action_evaluate_beta)
            self.btn_beta.setEnabled(False) 
            v_ops.addWidget(self.btn_beta)
            
            self.chk_numbers_first = QCheckBox("Use Numbers (0-9) after Z")
            v_ops.addWidget(self.chk_numbers_first)

            self.chk_only_numbers = QCheckBox("Use only numbers")
            v_ops.addWidget(self.chk_only_numbers)
            self.chk_numbers_first.toggled.connect(lambda state: state and self.chk_only_numbers.setChecked(False))
            self.chk_only_numbers.toggled.connect(lambda state: state and self.chk_numbers_first.setChecked(False))

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
            params_layout.setContentsMargins(5, 5, 5, 5)

            def create_param_input(default_val, tooltip_msg, is_int=False):
                line_edit = QLineEdit()
                line_edit.setText(str(default_val))
                if is_int: line_edit.setValidator(QIntValidator())
                else: line_edit.setValidator(QDoubleValidator())
                line_edit.editingFinished.connect(self.on_param_change) 
                line_edit.setToolTip(tooltip_msg) 
                return line_edit

            self.edit_z_min = create_param_input(4.6, "Minimum vertical distance to be considered a layer stack")
            params_layout.addRow("Z Shift Min (Å):", self.edit_z_min)
            self.edit_z_max = create_param_input(5.0, "Maximum vertical distance to be considered a layer stack")
            params_layout.addRow("Z Shift Max (Å):", self.edit_z_max)
            self.edit_xy_limit = create_param_input(3.0, "Maximum horizontal drift allowed between stacked subunits")
            params_layout.addRow("XY Shift Limit (Å):", self.edit_xy_limit)
            self.edit_neighbor_count = create_param_input(6, "How many closest chains to check for stacking", is_int=True)
            params_layout.addRow("Neighbor Search Count:", self.edit_neighbor_count)
            self.edit_water_z = create_param_input(2.5, "Maximum vertical distance from the protein layer centroids for a water/ion to be included")
            params_layout.addRow("Water/Ion Z Shift:", self.edit_water_z)

            self.grp_params.setLayout(params_layout)
            layout.addWidget(self.grp_params)

            self.text_area = QTextEdit()
            self.text_area.setReadOnly(True)
            font = self.text_area.font()
            font.setFamily("Courier New")
            font.setPointSize(9)
            self.text_area.setFont(font)
            layout.addWidget(self.text_area)

            # --- 1. SETUP THE RIGHT PANEL (3D VIEWER) ---
            right_widget = QWidget()
            r_layout = QHBoxLayout(right_widget)
            r_layout.setContentsMargins(0, 0, 0, 0)
            r_layout.setSpacing(2)

            self.mol_canvas = Mol3DCanvas(self, width=3, height=4, dpi=100)
            r_layout.addWidget(self.mol_canvas, 1)

            self.zoom_slider = QSlider(Qt.Orientation.Vertical)
            self.zoom_slider.setRange(1, 100)
            self.zoom_slider.setValue(50)
            self.zoom_slider.setTickPosition(QSlider.TickPosition.TicksRight)
            self.zoom_slider.setTickInterval(10)
            self.zoom_slider.valueChanged.connect(lambda v: self.mol_canvas.set_zoom_from_slider(v))
            self.mol_canvas.slider_callback = lambda v: self.zoom_slider.setValue(v)
            r_layout.addWidget(self.zoom_slider)

            # --- 2. WRAP BOTH SIDES IN A SPLITTER ---
            splitter = QSplitter(Qt.Orientation.Horizontal)
            splitter.addWidget(left_widget)
            splitter.addWidget(right_widget)
            
            splitter.setSizes([250, 300]) 
            splitter.setStretchFactor(0, 4)
            splitter.setStretchFactor(1, 6)

            main_layout.addWidget(splitter)
            self.setLayout(main_layout)
        
        def on_param_change(self):
            if self.working_file_path:
                self.text_area.append("\n--- Parameters changed: Re-calculating ---")
                QTimer.singleShot(10, lambda: self.process_pdb(self.working_file_path))

        def get_param(self, widget, default):
            try: return float(widget.text())
            except ValueError: return default
        
        def get_kept_heteroatoms(self, keep_layer_indices, xy_limit, water_z_limit):
            if not self.working_file_path or not self.final_sandwiches: return set()
            target_z_planes = []
            chain_centroids = {cid: self.results[cid]['centroid'] for cid in self.results}
            max_depth = self.detected_layers
            
            for layer_idx in keep_layer_indices:
                if layer_idx < 0 or layer_idx >= max_depth: continue
                layer_z_values = []
                for pf in self.final_sandwiches:
                    if layer_idx < len(pf):
                        cid = pf[layer_idx]
                        if cid in chain_centroids:
                            layer_z_values.append(chain_centroids[cid][2])
                if layer_z_values:
                    target_z_planes.append(sum(layer_z_values) / len(layer_z_values))

            if not target_z_planes: return set()

            het_atoms = []
            molecule_map = {}
            with open(self.working_file_path, 'r') as f:
                for i, line in enumerate(f):
                    if line.startswith("HETATM") or line.startswith("ATOM"):
                        res_name = line[17:20].strip()
                        chain_id = line[20:22].strip()
                        res_seq = line[22:26].strip()
                        if chain_id in self.results: continue 
                        try:
                            x, y, z = float(line[30:38]), float(line[38:46]), float(line[46:54])
                            atom_data = {'idx': i, 'x': x, 'y': y, 'z': z, 'res': res_name, 'chain': chain_id, 'res_seq': res_seq}
                            het_atoms.append(atom_data)
                            
                            mol_key = (chain_id, res_seq)
                            if mol_key not in molecule_map:
                                molecule_map[mol_key] = []
                            molecule_map[mol_key].append(i)
                        except: continue

            water_columns = [] 
            processed_indices = set()
            het_atoms.sort(key=lambda a: a['z'])
            
            for atom in het_atoms:
                if atom['idx'] in processed_indices: continue
                column = [atom]
                processed_indices.add(atom['idx'])
                matched_col = False
                for col in water_columns:
                    ref = col[0] 
                    if abs(atom['x'] - ref['x']) < 2.0 and abs(atom['y'] - ref['y']) < 2.0:
                        col.append(atom)
                        processed_indices.add(atom['idx'])
                        matched_col = True
                        break
                if not matched_col: water_columns.append(column)

            kept_line_indices = set()
            for col in water_columns:
                for target_z in target_z_planes:
                    best_atom = None
                    min_dist = float('inf')
                    for atom in col:
                        dist = abs(atom['z'] - target_z)
                        if dist < min_dist:
                            min_dist = dist
                            best_atom = atom
                    if best_atom and min_dist <= water_z_limit: 
                        mol_key = (best_atom['chain'], best_atom['res_seq'])
                        for mol_idx in molecule_map[mol_key]:
                            kept_line_indices.add(mol_idx)

            return kept_line_indices

        def populate_models(self):
            from chimerax.atomic import AtomicStructure
            current_models = self.session.models.list(type=AtomicStructure)
            current_ids = {m.id_string for m in current_models}
            
            if current_ids != self._last_model_ids:
                current_sel = self.cmb_models.currentData()
                current_sel_id = current_sel.id_string if current_sel else None
                
                self.cmb_models.blockSignals(True)
                self.cmb_models.clear()
                
                for model in current_models:
                    self.cmb_models.addItem(f"#{model.id_string} {model.name}", userData=model)
                
                if current_sel_id:
                    for i in range(self.cmb_models.count()):
                        m = self.cmb_models.itemData(i)
                        if m and m.id_string == current_sel_id:
                            self.cmb_models.setCurrentIndex(i)
                            break
                            
                self.cmb_models.blockSignals(False)
                self._last_model_ids = current_ids

        def load_from_chimerax(self):
            model = self.cmb_models.currentData()
            if not model: return
            
            try:
                self.original_filename_display = model.name
                folder = os.path.join(os.path.expanduser("~"), ".chimerax_modifier_temp")
                os.makedirs(folder, exist_ok=True)
                
                name_part, _ = os.path.splitext(self.original_filename_display)
                temp_path = os.path.join(folder, f"{name_part}_modified.pdb")
                
                run(self.session, f'save "{temp_path}" models #{model.id_string} format pdb')
                self.created_temp_files.append(temp_path)
                self.working_file_path = temp_path
                self.grp_ops.setEnabled(True)
                self.reload_working_model()
                self.text_area.clear()
                self.text_area.append(f"Created PDB working copy from {model.name}")
                self.check_spatial_and_id_duplicates(self.working_file_path)
                self.process_pdb(self.working_file_path)
            except Exception as e:
                self.text_area.append(f"Load Error: {e}")

        def reload_working_model(self):
            if self.working_model_id:
                try: run(self.session, f"close #{self.working_model_id}")
                except: pass
                
            try:
                models = run(self.session, f'open "{self.working_file_path}"')
                if models:
                    first_item = models[0]
                    first_model = first_item[0] if isinstance(first_item, (list, tuple)) else first_item
                        
                    if hasattr(first_model, 'id_string'):
                        self.working_model_id = first_model.id_string
                        run(self.session, f"color #{self.working_model_id} bypolymer")
                        run(self.session, f"view #{self.working_model_id}")
            except Exception as e:
                self.text_area.append(f"\nModel display warning: {e}")

        def action_rename(self):
            if not self.working_file_path: return
            try:
                old_sandwiches = self.final_sandwiches
                mapping = self.create_renaming_mapping(old_sandwiches)
                
                new_temp = self.working_file_path + ".tmp"
                self.write_renamed_pdb(self.working_file_path, new_temp, self.final_sandwiches)
                shutil.move(new_temp, self.working_file_path)
                
                self.text_area.append("\n>>>> Applied Renaming")
                self.reload_working_model()
                self.process_pdb(self.working_file_path, old_sandwiches=old_sandwiches, rename_map=mapping)
            except Exception as e: QMessageBox.critical(self, "Error", f"Rename failed: {e}")

        class ExpandDialog(QDialog):
            def __init__(self, current_layers, parent=None):
                super().__init__(parent)
                self.setWindowTitle("Expand Layers")
                self.current_layers = current_layers
                self.initUI()
                
            def initUI(self):
                layout = QVBoxLayout()
                self.chk_auto = QCheckBox("Auto Computing")
                self.chk_auto.setChecked(True)
                layout.addWidget(self.chk_auto)
                
                self.chk_alt = QCheckBox("Alternating Layers (A-B-A-B)")
                self.chk_alt.setToolTip("Expand using a 2-layer step to preserve alternating conformation differences.")
                layout.addWidget(self.chk_alt)
                
                form_layout = QFormLayout()
                self.edit_twist = QLineEdit("0.0")
                self.edit_rise = QLineEdit("4.8")
                self.edit_twist.setEnabled(False)
                self.edit_rise.setEnabled(False)
                self.edit_twist.setStyleSheet("color: gray;")
                self.edit_rise.setStyleSheet("color: gray;")
                
                form_layout.addRow("Twist (°):", self.edit_twist)
                form_layout.addRow("Rise (Å):", self.edit_rise)
                layout.addLayout(form_layout)
                
                self.btn_expand = QPushButton("Expand")
                self.btn_expand.clicked.connect(self.accept)
                layout.addWidget(self.btn_expand)
                self.setLayout(layout)
                
                self.chk_auto.stateChanged.connect(self.toggle_manual)
                if self.current_layers <= 1:
                    self.chk_auto.setChecked(False)
                    self.chk_auto.setEnabled(False)
                    self.chk_alt.setChecked(False)
                    self.chk_alt.setEnabled(False)
                    self.toggle_manual()
                    
            def toggle_manual(self):
                is_auto = self.chk_auto.isChecked()
                self.edit_twist.setEnabled(not is_auto)
                self.edit_rise.setEnabled(not is_auto)
                if is_auto:
                    self.edit_twist.setStyleSheet("color: gray;")
                    self.edit_rise.setStyleSheet("color: gray;")
                else:
                    self.edit_twist.setStyleSheet("")
                    self.edit_rise.setStyleSheet("")

        def action_trim(self):
            if not self.working_file_path: return
            target_layers = self.spin_trim.value()
            total_layers = self.detected_layers
            
            if target_layers == total_layers:
                self.text_area.append(f"\nYour model is already having {total_layers} layer(s)")
                return
            
            try:
                new_temp = self.working_file_path + ".tmp"
                
                if target_layers < total_layers:
                    start_idx = (total_layers - target_layers) // 2
                    end_idx = start_idx + target_layers
                    valid_layer_indices = range(start_idx, end_idx)
                    
                    keep_chains = set()
                    for sandwich in self.final_sandwiches:
                        for i, chain_id in enumerate(sandwich):
                            if i in valid_layer_indices: keep_chains.add(chain_id)
                    
                    self.text_area.append("Calculating water/ion retention based on layer centroids...")
                    xy_lim = self.get_param(self.edit_xy_limit, 3.0)
                    water_z_limit = self.get_param(self.edit_water_z, 4.0)
                    keep_het_lines = self.get_kept_heteroatoms(valid_layer_indices, xy_lim, water_z_limit)
                    
                    self.write_trimmed_pdb(self.working_file_path, new_temp, keep_chains, keep_het_lines)
                    shutil.move(new_temp, self.working_file_path)
                    self.text_area.append(f"\n>>>> Applied Trimming (Kept middle {target_layers} layers + associated waters)")
                
                else:
                    is_anti = getattr(self, 'chk_anti', None) and self.chk_anti.isChecked()
                    if is_anti:
                        QMessageBox.warning(self, "Action Restricted", "Expansion is currently not supported in Anti-parallel mode. You can only trim layers.")
                        return

                    dialog = self.ExpandDialog(total_layers, self)
                    if dialog.exec() != QDialog.DialogCode.Accepted:
                        return
                    
                    use_auto = dialog.chk_auto.isChecked()
                    use_alt = dialog.chk_alt.isChecked()
                    manual_twist = 0.0
                    manual_rise = 0.0
                    if not use_auto:
                        try:
                            manual_twist = float(dialog.edit_twist.text())
                            manual_rise = float(dialog.edit_rise.text())
                        except ValueError:
                            QMessageBox.warning(self, "Invalid Input", "Please enter valid numbers for Twist and Rise.")
                            return

                    self.text_area.append("Calculating expansion transformations...")
                    layers_to_add = target_layers - total_layers
                    water_z_limit = self.get_param(self.edit_water_z, 2.5)
                    self.write_expanded_pdb(self.working_file_path, new_temp, layers_to_add, use_auto, manual_twist, manual_rise, water_z_limit, use_alt)
                    shutil.move(new_temp, self.working_file_path)
                    self.text_area.append(f"\n>>>> Applied Expansion (Added {layers_to_add} layers, Alternating: {use_alt})")

                self.reload_working_model()
                self.process_pdb(self.working_file_path)
                
            except Exception as e: 
                QMessageBox.critical(self, "Error", f"Action failed: {e}")

        def get_transform(self, coords_from, coords_to):
            import numpy as np
            min_len = min(len(coords_from), len(coords_to))
            if min_len < 3:
                return np.eye(3), np.zeros(3)
            P = np.array(coords_from[:min_len])
            Q = np.array(coords_to[:min_len])
            centroid_P = np.mean(P, axis=0)
            centroid_Q = np.mean(Q, axis=0)
            P_centered = P - centroid_P
            Q_centered = Q - centroid_Q
            H = P_centered.T @ Q_centered
            U, S, Vt = np.linalg.svd(H)
            R = Vt.T @ U.T
            if np.linalg.det(R) < 0:
                Vt[2, :] *= -1
                R = Vt.T @ U.T
            t = centroid_Q - R @ centroid_P
            return R, t

        def write_expanded_pdb(self, input_path, output_path, layers_to_add, use_auto, manual_twist, manual_rise, water_z_limit=2.5, use_alt=False):
            import numpy as np
            import math
            top_layers_to_add = layers_to_add // 2
            bottom_layers_to_add = layers_to_add - top_layers_to_add

            used_chain_ids = set()
            max_res_seq = {} 
            
            orig_atoms_ters = []
            orig_hetatms = []
            headers = []
            orig_conects = []

            with open(input_path, 'r') as fin:
                for line in fin:
                    if line.startswith("MASTER") or line.startswith("END") or line.startswith("SEQRES") or line.startswith("SHEET") or line.startswith("HELIX"):
                        continue
                    elif line.startswith("ATOM") or line.startswith("TER") or line.startswith("ANISOU"):
                        orig_atoms_ters.append(line)
                        if len(line) >= 22:
                            chain_id = line[20:22].strip()
                            used_chain_ids.add(chain_id)
                    elif line.startswith("HETATM"):
                        orig_hetatms.append(line)
                    elif line.startswith("CONECT"):
                        orig_conects.append(line)
                    else:
                        headers.append(line)

            het_types = set()
            for line in orig_hetatms:
                if len(line) >= 20:
                    het_types.add(line[17:20].strip())
            
            reserved_het_chains = {}
            alphabet_backwards = "ZYXWVUTSRQPONMLKJIHGFEDCBA"
            for i, htype in enumerate(sorted(list(het_types))):
                if i < len(alphabet_backwards):
                    c = alphabet_backwards[i]
                    reserved_het_chains[htype] = c
                    used_chain_ids.add(c)
            
            updated_orig_hetatms = []
            het_counters = {c: 1 for c in reserved_het_chains.values()}
            het_molecule_map = {}
            
            for line in orig_hetatms:
                if len(line) < 26:
                    updated_orig_hetatms.append(line)
                    continue
                res_name = line[17:20].strip()
                if res_name in reserved_het_chains:
                    new_chain = reserved_het_chains[res_name]
                    orig_chain = line[20:22].strip()
                    try: orig_res_seq = int(line[22:26])
                    except: orig_res_seq = line[22:26]
                    
                    mol_key = (orig_chain, orig_res_seq, res_name)
                    if mol_key not in het_molecule_map:
                        het_molecule_map[mol_key] = het_counters[new_chain]
                        het_counters[new_chain] += 1
                        
                    new_seq = het_molecule_map[mol_key]
                    max_res_seq[new_chain] = new_seq
                    
                    l = list(line)
                    l[21] = new_chain; l[20] = ' '
                    l[22:26] = list(f"{new_seq:>4}")
                    updated_orig_hetatms.append("".join(l))
                else:
                    updated_orig_hetatms.append(line)
                    
            orig_hetatms = updated_orig_hetatms
            
            def get_new_chain_id():
                idx = 0
                while True:
                    label = self.generate_label(idx)
                    if label not in used_chain_ids:
                        used_chain_ids.add(label)
                        return label
                    idx += 1

            core_chains = set()
            for s in self.final_sandwiches:
                core_chains.update(s)

            chain_atoms = {}
            chain_ca_coords = {}
            standalone_atoms = orig_hetatms 
            
            for line in orig_atoms_ters:
                if line.startswith("ATOM") and len(line) >= 54:
                    chain_id = line[20:22].strip()
                    if chain_id in core_chains:
                        if chain_id not in chain_atoms:
                            chain_atoms[chain_id] = []
                            chain_ca_coords[chain_id] = []
                        chain_atoms[chain_id].append(line)
                        if line[12:16].strip() == "CA":
                            try:
                                x, y, z = float(line[30:38]), float(line[38:46]), float(line[46:54])
                                chain_ca_coords[chain_id].append([x, y, z])
                            except ValueError: pass

            chain_centroids = {}
            for cid, coords in chain_ca_coords.items():
                if coords:
                    chain_centroids[cid] = np.mean(coords, axis=0)
                    
            core_ca_coords = []
            for cid in core_chains:
                if cid in chain_ca_coords:
                    core_ca_coords.extend(chain_ca_coords[cid])
            global_centroid = np.mean(core_ca_coords, axis=0) if core_ca_coords else np.zeros(3)

            associated_atoms = {cid: [] for cid in core_chains}
            for line in standalone_atoms:
                if len(line) < 54: continue
                try:
                    x, y, z = float(line[30:38]), float(line[38:46]), float(line[46:54])
                    best_chain = None
                    min_dist = float('inf')
                    for cid, centroid in chain_centroids.items():
                        dz = abs(z - centroid[2])
                        if dz <= water_z_limit:
                            dist = math.sqrt((x - centroid[0])**2 + (y - centroid[1])**2 + (z - centroid[2])**2)
                            if dist < min_dist:
                                min_dist = dist
                                best_chain = cid
                    if best_chain:
                        associated_atoms[best_chain].append(line)
                except ValueError:
                    pass

            expansions = [] 
            for sandwich in self.final_sandwiches:
                if len(sandwich) == 0: continue
                
                bottom_chain = sandwich[0]
                top_chain = sandwich[-1]
                
                if use_alt:
                    # Alternating expansion (A-B-A-B). Uses a 2-layer step transform.
                    if use_auto:
                        if len(sandwich) >= 3:
                            # Calculate a 2-layer step directly
                            R_top, t_top = self.get_transform(chain_ca_coords[sandwich[-3]], chain_ca_coords[top_chain])
                            R_bottom, t_bottom = self.get_transform(chain_ca_coords[sandwich[2]], chain_ca_coords[bottom_chain])
                        elif len(sandwich) == 2:
                            # Square the 1-layer transform to approximate a 2-layer step if only 2 layers exist
                            R_t1, t_t1 = self.get_transform(chain_ca_coords[sandwich[-2]], chain_ca_coords[top_chain])
                            R_top = R_t1 @ R_t1
                            t_top = R_t1 @ t_t1 + t_t1
                            R_b1, t_b1 = self.get_transform(chain_ca_coords[sandwich[1]], chain_ca_coords[bottom_chain])
                            R_bottom = R_b1 @ R_b1
                            t_bottom = R_b1 @ t_b1 + t_b1
                        else:
                            R_top, t_top = np.eye(3), np.zeros(3)
                            R_bottom, t_bottom = np.eye(3), np.zeros(3)
                    else:
                        rad = math.radians(manual_twist * 2.0)
                        cos_t, sin_t = math.cos(rad), math.sin(rad)
                        rise = manual_rise * 2.0
                        R_top = np.array([[cos_t, -sin_t, 0], [sin_t, cos_t, 0], [0, 0, 1]])
                        base_t_top = np.array([0.0, 0.0, rise])
                        t_top = global_centroid - R_top @ global_centroid + base_t_top
                        
                        R_bottom = np.array([[cos_t, sin_t, 0], [-sin_t, cos_t, 0], [0, 0, 1]])
                        base_t_bottom = np.array([0.0, 0.0, -rise])
                        t_bottom = global_centroid - R_bottom @ global_centroid + base_t_bottom
                    
                    if len(sandwich) >= 2:
                        prev_b1 = sandwich[1]
                        prev_b2 = sandwich[0]
                    else:
                        prev_b1 = prev_b2 = sandwich[0]
                        
                    for _ in range(bottom_layers_to_add):
                        new_chain_id = get_new_chain_id()
                        expansions.append((prev_b1, new_chain_id, R_bottom, t_bottom))
                        prev_b1 = prev_b2
                        prev_b2 = new_chain_id

                    if len(sandwich) >= 2:
                        prev_t1 = sandwich[-2]
                        prev_t2 = sandwich[-1]
                    else:
                        prev_t1 = prev_t2 = sandwich[-1]
                        
                    for _ in range(top_layers_to_add):
                        new_chain_id = get_new_chain_id()
                        expansions.append((prev_t1, new_chain_id, R_top, t_top))
                        prev_t1 = prev_t2
                        prev_t2 = new_chain_id

                else:
                    if use_auto and len(sandwich) >= 2:
                        bottom_next = sandwich[1]
                        R_bottom, t_bottom = self.get_transform(chain_ca_coords[bottom_next], chain_ca_coords[bottom_chain])
                        top_prev = sandwich[-2]
                        R_top, t_top = self.get_transform(chain_ca_coords[top_prev], chain_ca_coords[top_chain])
                    else:
                        rad = math.radians(manual_twist)
                        cos_t = math.cos(rad)
                        sin_t = math.sin(rad)
                        R_top = np.array([[cos_t, -sin_t, 0], [sin_t, cos_t, 0], [0, 0, 1]])
                        base_t_top = np.array([0.0, 0.0, manual_rise])
                        t_top = global_centroid - R_top @ global_centroid + base_t_top
                        
                        R_bottom = np.array([[cos_t, sin_t, 0], [-sin_t, cos_t, 0], [0, 0, 1]])
                        base_t_bottom = np.array([0.0, 0.0, -manual_rise])
                        t_bottom = global_centroid - R_bottom @ global_centroid + base_t_bottom
                    
                    current_ref_chain = bottom_chain
                    for _ in range(bottom_layers_to_add):
                        new_chain_id = get_new_chain_id()
                        expansions.append((current_ref_chain, new_chain_id, R_bottom, t_bottom))
                        current_ref_chain = new_chain_id
                        
                    current_ref_chain = top_chain
                    for _ in range(top_layers_to_add):
                        new_chain_id = get_new_chain_id()
                        expansions.append((current_ref_chain, new_chain_id, R_top, t_top))
                        current_ref_chain = new_chain_id

            new_chain_lines = {}
            all_new_hetatms = []

            for ref_chain, new_chain, R, t in expansions:
                lines = chain_atoms.get(ref_chain, [])
                new_lines = []
                for line in lines:
                    try:
                        orig_x = float(line[30:38])
                        orig_y = float(line[38:46])
                        orig_z = float(line[46:54])
                        orig_vec = np.array([orig_x, orig_y, orig_z])
                        new_vec = R @ orig_vec + t
                        
                        l = list(line)
                        if len(new_chain) == 1:
                            l[21] = new_chain; l[20] = ' '
                        elif len(new_chain) >= 2:
                            l[20] = new_chain[0]; l[21] = new_chain[1]
                        
                        l[30:38] = list(f"{new_vec[0]:8.3f}")
                        l[38:46] = list(f"{new_vec[1]:8.3f}")
                        l[46:54] = list(f"{new_vec[2]:8.3f}")
                        
                        new_lines.append("".join(l))
                    except Exception:
                        new_lines.append(line)
                chain_atoms[new_chain] = new_lines
                new_chain_lines[new_chain] = new_lines

                het_lines = associated_atoms.get(ref_chain, [])
                het_res_map = {} 
                new_het_lines_for_next = []
                
                for line in het_lines:
                    try:
                        orig_x = float(line[30:38])
                        orig_y = float(line[38:46])
                        orig_z = float(line[46:54])
                        orig_vec = np.array([orig_x, orig_y, orig_z])
                        new_vec = R @ orig_vec + t
                        
                        l = list(line)
                        
                        orig_chain = line[20:22].strip()
                        orig_res_seq = int(line[22:26])
                        
                        res_key = (orig_chain, orig_res_seq)
                        if res_key not in het_res_map:
                            max_res_seq[orig_chain] = max_res_seq.get(orig_chain, 0) + 1
                            het_res_map[res_key] = max_res_seq[orig_chain]
                            
                        new_res_seq = het_res_map[res_key]
                        
                        if len(new_chain) == 1:
                            l[21] = new_chain; l[20] = ' '
                        elif len(new_chain) >= 2:
                            l[20] = new_chain[0]; l[21] = new_chain[1]
                            
                        l[22:26] = list(f"{new_res_seq:>4}")
                        l[30:38] = list(f"{new_vec[0]:8.3f}")
                        l[38:46] = list(f"{new_vec[1]:8.3f}")
                        l[46:54] = list(f"{new_vec[2]:8.3f}")
                        
                        new_line_str = "".join(l)
                        all_new_hetatms.append(new_line_str)
                        new_het_lines_for_next.append(new_line_str)
                    except Exception:
                        all_new_hetatms.append(line)
                        new_het_lines_for_next.append(line)
                
                associated_atoms[new_chain] = new_het_lines_for_next

            all_new_hetatms.sort(key=lambda line: (line[20:22], int(line[22:26]) if line[22:26].strip().isdigit() else 0))

            serial_counter = 1
            serial_map = {}

            with open(output_path, 'w') as fout:
                for line in headers:
                    fout.write(line)
                    
                for line in orig_atoms_ters:
                    try:
                        old_serial = int(line[6:11])
                        l = list(line)
                        l[6:11] = list(f"{serial_counter:>5}")
                        fout.write("".join(l))
                        serial_map[old_serial] = serial_counter
                        serial_counter += 1
                    except: fout.write(line)
                
                for ref_chain, new_chain, R, t in expansions:
                    last_line = ""
                    for line in new_chain_lines[new_chain]:
                        try:
                            l = list(line)
                            l[6:11] = list(f"{serial_counter:>5}")
                            fout.write("".join(l))
                            serial_counter += 1
                            last_line = line
                        except: fout.write(line)
                    
                    if last_line:
                        ter_line = list("TER   " + " " * 74)
                        ter_line[6:11] = list(f"{serial_counter:>5}")
                        ter_line[17:20] = list(last_line[17:20]) 
                        ter_line[20:22] = list(last_line[20:22]) 
                        ter_line[22:26] = list(last_line[22:26]) 
                        ter_line[26] = last_line[26]             
                        fout.write("".join(ter_line).rstrip() + "\n")
                        serial_counter += 1
                
                for line in orig_hetatms:
                    try:
                        old_serial = int(line[6:11])
                        l = list(line)
                        l[6:11] = list(f"{serial_counter:>5}")
                        fout.write("".join(l))
                        serial_map[old_serial] = serial_counter
                        serial_counter += 1
                    except: fout.write(line)

                for line in all_new_hetatms:
                    try:
                        l = list(line)
                        l[6:11] = list(f"{serial_counter:>5}")
                        fout.write("".join(l))
                        serial_counter += 1
                    except: fout.write(line)

                for line in orig_conects:
                    try:
                        parts = line.split()
                        if len(parts) < 2: continue
                        old_source = int(parts[1])
                        if old_source not in serial_map: continue
                        new_source = serial_map[old_source]
                        
                        valid_targets = []
                        for p in parts[2:]:
                            try:
                                tgt_old = int(p)
                                if tgt_old in serial_map:
                                    valid_targets.append(serial_map[tgt_old])
                            except: pass
                        
                        if valid_targets:
                            out_line = "CONECT" + f"{new_source:>5}"
                            for tgt in valid_targets:
                                out_line += f"{tgt:>5}"
                            fout.write(out_line + "\n")
                    except: pass
                    
                fout.write("END   \n")

        def action_restrain(self):
            if not self.final_sandwiches or self.detected_layers == 0: return
            file_path, _ = QFileDialog.getSaveFileName(self, "Save Restrain File", "restrain_torsion.cxc", "ChimeraX Command (*.cxc);;All Files (*)")
            if file_path:
                try:
                    mapping = self.create_renaming_mapping(self.final_sandwiches)
                    mid_layer_idx = self.detected_layers // 2 
                    with open(file_path, 'w') as f:
                        for sandwich in self.final_sandwiches:
                            mid_chain_original = sandwich[mid_layer_idx]
                            mid_id_new = mapping[mid_chain_original]
                            for layer_k in range(self.detected_layers):
                                if layer_k == mid_layer_idx: continue 
                                target_chain_original = sandwich[layer_k]
                                target_id_new = mapping[target_chain_original]
                                f.write(f"isolde restrain torsions /{target_id_new} template /{mid_id_new} angleRange 180\n")
                    self.text_area.append(f"\n[Export] Restrain file generated: {os.path.basename(file_path)}")
                except Exception as e: QMessageBox.critical(self, "Error", f"Failed to save restrain file:\n{str(e)}")

        def action_evaluate_beta(self):
            if not self.working_file_path or self.detected_layers <= 1: return
            if not self.working_model_id: return
            
            try:
                self.text_area.append(f"\n>>> Re-evaluating secondary structure using ChimeraX dssp...")
                run(self.session, f"dssp #{self.working_model_id}")
                
                run(self.session, f'save "{self.working_file_path}" models #{self.working_model_id} format pdb')
                self.reload_working_model()
                self.process_pdb(self.working_file_path)
                self.text_area.append(f">>> β-sheet evaluation completed \n>>> Reloaded PDB file")
                    
            except Exception as e:
                self.text_area.append(f"ChimeraX DSSP Error: {e}")

        def check_spatial_and_id_duplicates(self, file_path):
            split_line_indices = set()
            segments_per_chain_id = {} 
            last_chain_id = None
            last_ca_coord = None
            
            with open(file_path, 'r') as f:
                for i, line in enumerate(f):
                    if line.startswith("ATOM"):
                        current_chain_id = line[20:22].strip() 
                        is_ca = line[12:16].strip() == "CA"
                        current_coord = None
                        if is_ca:
                            try: current_coord = (float(line[30:38]), float(line[38:46]), float(line[46:54]))
                            except: pass

                        is_new_segment = False
                        if current_chain_id != last_chain_id:
                            is_new_segment = True
                            last_ca_coord = None 
                        
                        if is_new_segment:
                            split_line_indices.add(i)
                            if current_chain_id not in segments_per_chain_id:
                                segments_per_chain_id[current_chain_id] = 0
                            segments_per_chain_id[current_chain_id] += 1
                            last_chain_id = current_chain_id
                        
                        if is_ca: last_ca_coord = current_coord

            has_duplicates = any(count > 1 for count in segments_per_chain_id.values())
            if has_duplicates:
                reply = QMessageBox.question(
                    self, 
                    "Broken Chains / Duplicates Detected", 
                    "Detected multiple segments using the same Chain ID.\n"
                    "(Based on file order OR spatial gaps > 15Å)\n\n"
                    "Rename all segments to unique IDs to proceed?",
                    QMessageBox.StandardButton.Yes | QMessageBox.StandardButton.No, 
                    QMessageBox.StandardButton.Yes
                )
                if reply == QMessageBox.StandardButton.Yes:
                    self.rewrite_split_chains(file_path, split_line_indices)
                    self.text_area.append(">>> Auto-corrected duplicate/broken chains.")

        def rewrite_split_chains(self, file_path, split_indices):
            temp_out = file_path + ".tmp"
            current_segment_idx = -1 
            current_label = "A"
            with open(file_path, 'r') as fin, open(temp_out, 'w') as fout:
                for i, line in enumerate(fin):
                    if i in split_indices:
                        current_segment_idx += 1
                        current_label = self.generate_label(current_segment_idx)
                    if line.startswith("ATOM") or line.startswith("HETATM"):
                        line_list = list(line)
                        if len(current_label) == 1: 
                            line_list[21] = current_label
                            line_list[20] = ' '
                        elif len(current_label) >= 2: 
                            line_list[20] = current_label[0]
                            line_list[21] = current_label[1]
                        fout.write("".join(line_list))
                    else:
                        fout.write(line)
            shutil.move(temp_out, file_path)

        def process_pdb(self, filename, old_sandwiches=None, rename_map=None):
            Z_MIN = self.get_param(self.edit_z_min, 4.6)
            Z_MAX = self.get_param(self.edit_z_max, 5.0)
            XY_SHIFT_LIMIT = self.get_param(self.edit_xy_limit, 3.0)
            NEIGHBOR_SEARCH_COUNT = int(self.get_param(self.edit_neighbor_count, 6))
            
            try:
                raw_chains, waters, ions = self.parse_pdb_manual_logic(filename)
                if not raw_chains:
                    self.text_area.append("Error: No Alpha Carbons found.")
                    return
                
                chain_ids = list(raw_chains.keys())
                total_chains = len(chain_ids)

                self.chains_data_plot = {}
                for cid, residues in raw_chains.items():
                    self.chains_data_plot[cid] = [tuple(r['coord']) for r in residues]
                self.current_waters = waters
                self.current_ions = ions

                is_anti = getattr(self, 'chk_anti', None) and self.chk_anti.isChecked()
                
                if is_anti:
                    anti_pairing_max_xy = XY_SHIFT_LIMIT * 2.0
                    anti_pairing_max_z = Z_MAX * 0.7
                else:
                    anti_pairing_max_xy = 0
                    anti_pairing_max_z = 0

                base_centroids = {cid: self.get_centroid(raw_chains[cid]) for cid in chain_ids}
                unit_ids = []
                unit_chains_map = {} 
                
                if is_anti:
                    self.text_area.append("Anti mode: Grouping chains into structural pairs...")
                    paired = set()
                    for cid in chain_ids:
                        if cid in paired: continue
                        best_match = None
                        min_dist = float('inf')
                        c1 = base_centroids[cid]
                        for other_cid in chain_ids:
                            if other_cid == cid or other_cid in paired: continue
                            c2 = base_centroids[other_cid]
                            d_xy = math.hypot(c1[0] - c2[0], c1[1] - c2[1])
                            d_z = abs(c1[2] - c2[2])
                            if d_xy > anti_pairing_max_xy or d_z > anti_pairing_max_z: continue
                            d = np.linalg.norm(c1 - c2)
                            if d < min_dist:
                                min_dist = d; best_match = other_cid
                        if best_match:
                            paired.add(cid); paired.add(best_match)
                            u_id = f"{cid}|{best_match}"
                            unit_ids.append(u_id); unit_chains_map[u_id] = [cid, best_match]
                        else:
                            unit_ids.append(cid); unit_chains_map[cid] = [cid]
                else:
                    unit_ids = chain_ids[:]
                    for cid in chain_ids: unit_chains_map[cid] = [cid]

                self.results = {}
                for u_id in unit_ids:
                    chains_in_unit = unit_chains_map[u_id]
                    avg_c = np.mean([base_centroids[c] for c in chains_in_unit], axis=0)
                    self.results[u_id] = {
                        'centroid': avg_c, 'protofilament': None, 'layer': None,
                        'neighbors': {'up': None, 'down': None}, 'original_chains': chains_in_unit
                    }

                self.text_area.append("Identifying Protofilaments...")
                adjacency = {u_id: [] for u_id in unit_ids}
                
                for id_a in unit_ids:
                    center_a = self.results[id_a]['centroid']
                    distances = []
                    for id_b in unit_ids:
                        if id_a == id_b: continue
                        dist = np.linalg.norm(center_a - self.results[id_b]['centroid'])
                        distances.append((dist, id_b))
                    
                    distances.sort(key=lambda x: x[0])
                    closest = [x[1] for x in distances[:NEIGHBOR_SEARCH_COUNT]]
                    
                    for id_b in closest:
                        best_relation = 0
                        for c_a in self.results[id_a]['original_chains']:
                            for c_b in self.results[id_b]['original_chains']:
                                rel = self.check_stacking(raw_chains[c_a], raw_chains[c_b], Z_MIN, Z_MAX, XY_SHIFT_LIMIT)
                                if rel != 0: best_relation = rel; break
                            if best_relation != 0: break
                        
                        relation = best_relation
                        if relation == 1: 
                            self.results[id_a]['neighbors']['up'] = id_b
                            self.results[id_b]['neighbors']['down'] = id_a
                            adjacency[id_a].append(id_b); adjacency[id_b].append(id_a)
                        elif relation == -1: 
                            self.results[id_a]['neighbors']['down'] = id_b
                            self.results[id_b]['neighbors']['up'] = id_a
                            adjacency[id_a].append(id_b); adjacency[id_b].append(id_a)

                protofilaments = []
                visited = set()
                for u_id in unit_ids:
                    if u_id not in visited:
                        component = []
                        stack = [u_id]
                        visited.add(u_id)
                        while stack:
                            curr = stack.pop()
                            component.append(curr)
                            for neighbor in adjacency[curr]:
                                if neighbor not in visited:
                                    visited.add(neighbor)
                                    stack.append(neighbor)
                        protofilaments.append(component)

                protofilaments.sort(key=len, reverse=True)
                for idx, group in enumerate(protofilaments):
                    pf_id = idx + 1
                    for u_id in group: self.results[u_id]['protofilament'] = pf_id

                self.text_area.append("Aligning Layers...")
                layer_map = {} 
                if protofilaments:
                    ref_pf = protofilaments[0]
                    anchor_unit = ref_pf[len(ref_pf)//2]
                    layer_map[anchor_unit] = 0
                    anchor_centroid = self.results[anchor_unit]['centroid']
                    initial_layer_units = [anchor_unit]

                    for pf in protofilaments[1:]:
                        best_match = None
                        min_dist = float('inf')
                        for u_id in pf:
                            dist = np.linalg.norm(self.results[u_id]['centroid'] - anchor_centroid)
                            if dist < min_dist: min_dist = dist; best_match = u_id
                        if best_match:
                            layer_map[best_match] = 0
                            initial_layer_units.append(best_match)

                    queue = initial_layer_units[:] 
                    processed = set(initial_layer_units)
                    
                    while queue:
                        curr = queue.pop(0)
                        curr_layer = layer_map[curr]
                        up = self.results[curr]['neighbors']['up']
                        if up and up not in processed:
                            layer_map[up] = curr_layer + 1; processed.add(up); queue.append(up)
                        down = self.results[curr]['neighbors']['down']
                        if down and down not in processed:
                            layer_map[down] = curr_layer - 1; processed.add(down); queue.append(down)
                    
                    if layer_map:
                        min_L = min(layer_map.values())
                        offset = 1 - min_L
                        for k in layer_map: layer_map[k] += offset

                unpacked_protofilaments = []
                unpacked_layer_map = {}
                
                for pf in protofilaments:
                    if is_anti:
                        pf_chain1, pf_chain2 = [], []
                        for u_id in pf:
                            chains = self.results[u_id]['original_chains']
                            c1 = chains[0]; c2 = chains[1] if len(chains) > 1 else None
                            pf_chain1.append(c1)
                            if c2: pf_chain2.append(c2)
                            if u_id in layer_map:
                                unpacked_layer_map[c1] = layer_map[u_id]
                                if c2: unpacked_layer_map[c2] = layer_map[u_id]
                        if pf_chain1: unpacked_protofilaments.append(pf_chain1)
                        if pf_chain2: unpacked_protofilaments.append(pf_chain2)
                    else:
                        pf_chains = []
                        for u_id in pf:
                            c1 = self.results[u_id]['original_chains'][0]
                            pf_chains.append(c1)
                            if u_id in layer_map: unpacked_layer_map[c1] = layer_map[u_id]
                        if pf_chains: unpacked_protofilaments.append(pf_chains)

                protofilaments = unpacked_protofilaments
                layer_map = unpacked_layer_map

                final_results = {}
                for u_id in unit_ids:
                    for c_id in self.results[u_id]['original_chains']:
                        final_results[c_id] = {'centroid': base_centroids[c_id]}
                self.results = final_results

                self.final_sandwiches = []
                if protofilaments:
                    pf_centroids = []
                    for pf in protofilaments:
                        group_centroids = [self.results[c]['centroid'] for c in pf]
                        pf_centroids.append(np.mean(group_centroids, axis=0))
                    
                    global_center = np.mean(pf_centroids, axis=0)
                    indexed_centroids = list(enumerate(pf_centroids))
                    start_idx, start_centroid = min(indexed_centroids, key=lambda p: p[1][0])
                    sorted_indices = [start_idx]; visited_indices = {start_idx}
                    
                    def get_xy_dist(c1, c2): return math.sqrt((c1[0]-c2[0])**2 + (c1[1]-c2[1])**2)
                    def find_symmetry(ref_idx, tolerance=15.0):
                        ref_c = pf_centroids[ref_idx]
                        vec_x = ref_c[0] - global_center[0]; vec_y = ref_c[1] - global_center[1]
                        target_x = global_center[0] - vec_x; target_y = global_center[1] - vec_y
                        target = (target_x, target_y)
                        best_sym_idx = -1; min_sym_dist = float('inf')
                        for idx, c in indexed_centroids:
                            if idx in visited_indices: continue
                            d = math.sqrt((c[0]-target[0])**2 + (c[1]-target[1])**2)
                            if d < min_sym_dist: min_sym_dist = d; best_sym_idx = idx
                        if best_sym_idx != -1 and min_sym_dist <= tolerance: return best_sym_idx
                        return -1

                    sym_a = find_symmetry(start_idx)
                    if sym_a != -1: sorted_indices.append(sym_a); visited_indices.add(sym_a)

                    remaining = []
                    for idx, c in indexed_centroids:
                        if idx not in visited_indices:
                            d = get_xy_dist(start_centroid, c)
                            remaining.append((d, idx))
                    
                    remaining.sort(key=lambda x: x[0]) 
                    
                    for dist, idx in remaining:
                        if idx in visited_indices: continue
                        sorted_indices.append(idx); visited_indices.add(idx)
                        sym_partner = find_symmetry(idx)
                        if sym_partner != -1: sorted_indices.append(sym_partner); visited_indices.add(sym_partner)

                    protofilaments = [protofilaments[i] for i in sorted_indices]

                for group in protofilaments:
                    group_with_layers = [c for c in group if c in layer_map]
                    if not group_with_layers: continue
                    sorted_group = sorted(group_with_layers, key=lambda c: layer_map[c])
                    self.final_sandwiches.append(sorted_group)

                if self.final_sandwiches: self.detected_layers = max(len(s) for s in self.final_sandwiches)
                else: self.detected_layers = 0
                
                self.btn_beta.setEnabled(self.detected_layers > 1)

                output_lines = []
                output_lines.append("--------------- Analysis Result ---------------")
                output_lines.append(f"Total Chains: {total_chains}")
                output_lines.append(f"Max Layers: {self.detected_layers}")
                
                is_anti = getattr(self, 'chk_anti', None) and self.chk_anti.isChecked()
                
                if rename_map is not None and old_sandwiches is not None: display_sandwiches = old_sandwiches; is_renaming = True
                else: display_sandwiches = self.final_sandwiches; is_renaming = False

                if is_anti:
                    merged_pfs = []
                    i = 0
                    while i < len(display_sandwiches):
                        s1 = display_sandwiches[i]
                        if i + 1 < len(display_sandwiches) and len(s1) == len(display_sandwiches[i+1]):
                            s2 = display_sandwiches[i+1]
                            merged_old = [f"({c1},{c2})" for c1, c2 in zip(s1, s2)]
                            if is_renaming: merged_new = [f"({rename_map.get(c1, '?')},{rename_map.get(c2, '?')})" for c1, c2 in zip(s1, s2)]
                            else: merged_new = []
                            merged_pfs.append((s1[len(s1)//2], merged_old, merged_new))
                            i += 2 
                        else:
                            merged_pfs.append((s1[len(s1)//2], list(s1), []))
                            i += 1
                    
                    output_lines.append(f"Protofilaments: {len(merged_pfs)}")
                    for i, (mid_chain, m_old, m_new) in enumerate(merged_pfs):
                        if is_renaming:
                            new_mid = rename_map.get(mid_chain, "?")
                            output_lines.append(f"\nAnti-Protofilament #{i+1} (Contains Chain {mid_chain} -> {new_mid})")
                            output_lines.append(f"Old Units: {', '.join(m_old)}")
                            output_lines.append(f"New Units: {', '.join(m_new)}")
                        else:
                            output_lines.append(f"\nAnti-Protofilament #{i+1} (Contains Chain {mid_chain})")
                            output_lines.append(", ".join(m_old))
                else:
                    output_lines.append(f"Protofilaments: {len(display_sandwiches)}")
                    for i, sandwich in enumerate(display_sandwiches):
                        mid_chain = sandwich[len(sandwich)//2]
                        if is_renaming:
                            new_mid = rename_map.get(mid_chain, "?")
                            output_lines.append(f"\nProtofilament #{i+1} (Contains Chain {mid_chain} -> {new_mid})")
                            output_lines.append(f"Old: {', '.join(sandwich)}")
                            new_seq = [rename_map.get(c, "?") for c in sandwich]
                            output_lines.append(f"New: {', '.join(new_seq)}")
                        else:
                            output_lines.append(f"\nProtofilament #{i+1} (Contains Chain {mid_chain})")
                            output_lines.append(", ".join(sandwich))

                self.text_area.append("\n".join(output_lines))
                
                highlight_ids = set()
                if self.final_sandwiches:
                    mid_idx = self.detected_layers // 2
                    for s in self.final_sandwiches:
                        if mid_idx < len(s): 
                            highlight_ids.add(s[mid_idx])
                
                self.mol_canvas.plot_chains(
                    self.chains_data_plot, 
                    label_ids=highlight_ids, 
                    waters=self.current_waters, 
                    ions=self.current_ions
                )

            except Exception as e:
                self.text_area.append(f"Analysis Error: {str(e)}")

        def parse_pdb_manual_logic(self, filepath):
            chains = {}; waters = []; ions = []
            with open(filepath, 'r') as f:
                for line in f:
                    if line.startswith("ATOM") or line.startswith("HETATM"):
                        atom_name = line[12:16].strip()
                        res_name = line[17:20].strip()
                        try: x, y, z = float(line[30:38]), float(line[38:46]), float(line[46:54])
                        except ValueError: continue

                        if "CA" == atom_name or "CA " in line[12:16]:
                            chain_id = line[20:22].strip()
                            res_id = line[22:27].strip()
                            if chain_id not in chains: chains[chain_id] = []
                            chains[chain_id].append({'id': res_id, 'coord': np.array([x, y, z])})
                        elif line.startswith("HETATM"):
                            if res_name in ["HOH", "WAT"]:
                                if atom_name.startswith("O"): waters.append((x, y, z))
                            else: ions.append((x, y, z))
            return chains, waters, ions

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
            n_layers = self.detected_layers
            layer_order = self.get_expansion_order(n_layers)
            
            skip_chars = set()
            if self.working_file_path:
                het_types = set()
                try:
                    with open(self.working_file_path, 'r') as f:
                        for line in f:
                            if line.startswith("HETATM") and len(line) >= 20:
                                het_types.add(line[17:20].strip())
                    alphabet_backwards = "ZYXWVUTSRQPONMLKJIHGFEDCBA"
                    for i, htype in enumerate(sorted(list(het_types))):
                        if i < len(alphabet_backwards):
                            skip_chars.add(alphabet_backwards[i])
                except: pass

            mapping = {}; global_index = 0
            for layer_idx in layer_order:
                for sandwich in sandwiches:
                    if layer_idx < len(sandwich):
                        chain_id = sandwich[layer_idx]
                        while True:
                            new_label = self.generate_label(global_index)
                            global_index += 1
                            if new_label not in skip_chars:
                                break
                        mapping[chain_id] = new_label
            return mapping

        def write_renamed_pdb(self, input_filename, output_filename, sandwiches):
            prot_mapping = self.create_renaming_mapping(sandwiches)
            het_types = set()
            with open(input_filename, 'r') as f:
                for line in f:
                    if line.startswith("HETATM") and len(line) >= 20:
                        het_types.add(line[17:20].strip())
            
            het_type_map = {}
            alphabet_backwards = "ZYXWVUTSRQPONMLKJIHGFEDCBA"
            for i, htype in enumerate(sorted(list(het_types))):
                if i < len(alphabet_backwards):
                    het_type_map[htype] = alphabet_backwards[i]

            het_counters = {cid: 1 for cid in het_type_map.values()}
            last_seen_residue = {cid: (None, None) for cid in het_type_map.values()}

            with open(input_filename, 'r') as fin, open(output_filename, 'w') as fout:
                for line in fin:
                    record = line[0:6].strip()
                    is_het = (record == "HETATM")
                    is_anisou = (record == "ANISOU")
                    res_name = line[17:20].strip() if len(line) >= 20 else ""

                    if (is_het) or (is_anisou and res_name in het_type_map):
                        if res_name in het_type_map:
                            new_chain = het_type_map[res_name]
                            current_input_key = line[20:27] 
                            last_key, last_seq = last_seen_residue[new_chain]
                            if current_input_key == last_key: new_seq = last_seq
                            else:
                                new_seq = het_counters[new_chain]; het_counters[new_chain] += 1
                                last_seen_residue[new_chain] = (current_input_key, new_seq)
                            
                            l = list(line)
                            if len(new_chain) == 1: l[21] = new_chain; l[20] = ' ' 
                            elif len(new_chain) >= 2: l[20] = new_chain[0]; l[21] = new_chain[1]
                            new_seq_str = f"{new_seq:>4}"[-4:] 
                            l[22:26] = list(new_seq_str)
                            fout.write("".join(l))
                            continue 
                    
                    if record in self.CHAIN_RECORD_SPECS:
                        line_chars = list(line)
                        slice_list = self.CHAIN_RECORD_SPECS[record]
                        for start, end in slice_list:
                            if end <= len(line):
                                old_id = line[start:end].strip()
                                if old_id in prot_mapping:
                                    new_id = prot_mapping[old_id]
                                    width = end - start
                                    if len(new_id) <= width:
                                        formatted = f"{new_id:>{width}}"
                                        line_chars[start:end] = list(formatted)
                        fout.write("".join(line_chars))
                    else:
                        fout.write(line)

        def write_trimmed_pdb(self, input_path, output_path, keep_chains, keep_atom_indices=None):
            if keep_atom_indices is None: keep_atom_indices = set()
            serial_map = {}; new_serial_counter = 1
            master_counts = {
                "numRemark": 0, "numHet": 0, "numHelix": 0, "numSheet": 0,
                "numTurn": 0, "numSite": 0, "numXform": 0, "numCoord": 0,
                "numTer": 0, "numConect": 0, "numSeq": 0
            }
            sheet_buffer = []; helix_counter = 1

            def update_master_count(record_name):
                if record_name == "REMARK": master_counts["numRemark"] += 1
                elif record_name == "HET":    master_counts["numHet"] += 1
                elif record_name == "HELIX":  master_counts["numHelix"] += 1
                elif record_name == "SHEET":  master_counts["numSheet"] += 1
                elif record_name == "TURN":   master_counts["numTurn"] += 1
                elif record_name == "SITE":   master_counts["numSite"] += 1
                elif record_name in ["ATOM", "HETATM"]: master_counts["numCoord"] += 1
                elif record_name == "TER":    master_counts["numTer"] += 1
                elif record_name == "CONECT": master_counts["numConect"] += 1
                elif record_name == "SEQRES": master_counts["numSeq"] += 1
                elif record_name in ["ORIGX1", "ORIGX2", "ORIGX3", "SCALE1", "SCALE2", "SCALE3", "MTRIX1", "MTRIX2", "MTRIX3"]:
                    master_counts["numXform"] += 1

            with open(input_path, 'r') as fin, open(output_path, 'w') as fout:
                for i, line in enumerate(fin):
                    record = line[0:6].strip()
                    if record in ["MASTER", "END"]: continue

                    if record == "SHEET":
                        sheet_buffer.append(line)
                        continue
                    if sheet_buffer and record != "SHEET":
                        self.process_and_write_sheets(sheet_buffer, keep_chains, fout, master_counts)
                        sheet_buffer = []

                    if record == "HELIX":
                        if len(line) > 31:
                            c1, c2 = line[19:21].strip(), line[31:33].strip()
                            if c1 in keep_chains and c2 in keep_chains:
                                new_id_str = f"{helix_counter:>3}"
                                new_line = line[:7] + new_id_str + line[10:]
                                fout.write(new_line)
                                update_master_count("HELIX")
                                helix_counter += 1
                        continue

                    if record in ["ATOM", "HETATM", "TER", "ANISOU"]:
                        if len(line) < 22: continue
                        chain_id = line[20:22].strip() 
                        is_kept_protein = chain_id in keep_chains
                        is_kept_water = i in keep_atom_indices
                        
                        if is_kept_protein or is_kept_water:
                            try:
                                old_serial = int(line[6:11])
                                if record == "ANISOU" and old_serial in serial_map: current_new_id = serial_map[old_serial]
                                else:
                                    current_new_id = new_serial_counter
                                    serial_map[old_serial] = current_new_id
                                    if record != "ANISOU": new_serial_counter += 1
                                new_line = line[:6] + f"{current_new_id:>5}" + line[11:]
                                fout.write(new_line)
                                update_master_count(record)
                            except ValueError:
                                fout.write(line)
                                update_master_count(record)
                        continue

                    if record == "CONECT":
                        try:
                            parts = line.split() 
                            if len(parts) < 2: continue
                            old_source = int(parts[1])
                            if old_source not in serial_map: continue
                            new_source = serial_map[old_source]
                            
                            valid_targets = []
                            for p in parts[2:]:
                                try:
                                    if int(p) in serial_map: valid_targets.append(serial_map[int(p)])
                                except: pass
                            
                            if not valid_targets: continue
                                
                            out_line = "CONECT" + f"{new_source:>5}"
                            for tgt in valid_targets: out_line += f"{tgt:>5}"
                            fout.write(out_line + "\n")
                            update_master_count("CONECT")
                        except: pass
                        continue

                    if record == "LINK":
                        try:
                            if len(line) > 51:
                                if line[21] in keep_chains and line[51] in keep_chains: fout.write(line)
                        except: pass
                        continue

                    should_write = True
                    if record in self.CHAIN_RECORD_SPECS:
                        slice_list = self.CHAIN_RECORD_SPECS[record]
                        chains_in_line = []
                        for start, end in slice_list:
                            if end <= len(line):
                                c = line[start:end].strip()
                                if c: chains_in_line.append(c)
                        if chains_in_line:
                            if not all(c in keep_chains for c in chains_in_line): should_write = False
                    
                    if should_write:
                        fout.write(line)
                        update_master_count(record)

                if sheet_buffer: self.process_and_write_sheets(sheet_buffer, keep_chains, fout, master_counts)

                master_line = "MASTER    {:>5}{:>5}{:>5}{:>5}{:>5}{:>5}{:>5}{:>5}{:>5}{:>5}{:>5}{:>5}".format(
                    master_counts["numRemark"], "0", master_counts["numHet"], master_counts["numHelix"],
                    master_counts["numSheet"], master_counts["numTurn"], master_counts["numSite"],
                    master_counts["numXform"], master_counts["numCoord"], master_counts["numTer"],
                    master_counts["numConect"], master_counts["numSeq"]
                )
                fout.write(master_line + "\n")
                fout.write("END   \n")

        def process_and_write_sheets(self, sheet_lines, keep_chains, fout, master_counts):
            from collections import OrderedDict
            grouped_sheets = OrderedDict()
            for line in sheet_lines:
                if len(line) < 22: continue
                chain_id = line[20:22].strip() 
                if chain_id not in keep_chains: continue
                try: sheet_id = line.split()[2]
                except: sheet_id = line[11:14].strip()
                
                if sheet_id not in grouped_sheets: grouped_sheets[sheet_id] = []
                grouped_sheets[sheet_id].append(line)
            
            new_sheet_id_counter = 1
            for old_id, lines in grouped_sheets.items():
                total_strands = len(lines)
                if total_strands == 0: continue
                new_sheet_id_str = f"{new_sheet_id_counter:>3}"
                new_num_strands_str = f"{total_strands:>2}"
                current_strand_id = 1 
                for line in lines:
                    new_strand_id_str = f"{current_strand_id:>3}"
                    chars = list(line)
                    chars[7:10] = list(new_strand_id_str)
                    chars[11:14] = list(new_sheet_id_str)
                    chars[14:16] = list(new_num_strands_str)
                    
                    if current_strand_id == 1:
                        if len(chars) > 40: chars[38:40] = list(" 0")
                        if len(chars) > 41:
                            limit = min(len(chars), 70)
                            for i in range(41, limit): chars[i] = ' '
                    
                    fout.write("".join(chars))
                    master_counts["numSheet"] += 1
                    current_strand_id += 1
                new_sheet_id_counter += 1

    # ================== Master Tab Wrapper Class ==================
    class ModifierToolTabs(QTabWidget):
        def __init__(self, tool_instance, session):
            super().__init__()
            self.pdb_widget = PDBLayerIdentifier(tool_instance, session)
            self.cif_widget = CIFLayerIdentifier(tool_instance, session)
            self.addTab(self.pdb_widget, "PDB Mode")
            self.addTab(self.cif_widget, "CIF Mode")

    return ModifierToolTabs


class PDBModifierTool(ToolInstance):
    SESSION_ENDURING = False
    
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
        
        self.widget.setMinimumWidth(350)
        
        self.tool_window.manage("right")
        
        self.widget.setStyleSheet("""
            QMainWindow, QDialog, QWidget { background-color: #1e1e1e; color: #d4d4d4; font-family: Arial; font-size: 10pt; }
            QPushButton { background-color: #3e3e42; color: #d4d4d4; border: 1px solid #3e3e42; padding: 5px 15px; border-radius: 4px; }
            QPushButton:hover { background-color: #4e4e52; border: 1px solid #98c379; }
            QPushButton:pressed, QPushButton:checked { background-color: #98c379; color: #1e1e1e; border: 1px solid #98c379; }
            QLineEdit, QSpinBox, QTextEdit, QTextBrowser, QTableWidget { background-color: #3c3c3c; border: 1px solid #3c3c3c; color: #cccccc; padding: 4px;}
            QGroupBox { border: 1px solid #3e3e42; border-radius: 4px; margin-top: 1.5em; font-weight: bold; color: #98c379; }
            QGroupBox::title { subcontrol-origin: margin; subcontrol-position: top left; padding: 0 5px; }
            QCheckBox::indicator { width: 14px; height: 14px; border: 1px solid #555; border-radius: 2px; background-color: #1e1e1e; }
            QCheckBox::indicator:checked { background-color: #98c379; border: 1px solid #98c379; }
            QSlider:vertical { min-width: 20px; }
            QSlider::groove:vertical { background: #3c3c3c; width: 6px; border-radius: 3px; }
            QSlider::handle:vertical { background: #98c379; height: 14px; width: 14px; margin: 0 -4px; border-radius: 7px; }
            QSlider::handle:vertical:hover { background: #b5e890; }
            QTabWidget::pane { border: 1px solid #3e3e42; top: -1px; }
            QTabBar::tab { background: #252526; border: 1px solid #3e3e42; padding: 6px 12px; color: #d4d4d4; }
            QTabBar::tab:selected { background: #3e3e42; color: #98c379; font-weight: bold; border-bottom: 1px solid #3e3e42; }
            QTabBar::tab:hover { background: #4e4e52; }
        """)
        
    def delete(self):
        import os
        for sub_widget in [self.widget.pdb_widget, self.widget.cif_widget]:
            if hasattr(sub_widget, 'created_temp_files'):
                for f in sub_widget.created_temp_files:
                    if os.path.exists(f):
                        try: os.remove(f)
                        except: pass
                    
        super().delete()
