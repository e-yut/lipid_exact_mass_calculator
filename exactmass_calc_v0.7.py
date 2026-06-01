import tkinter as tk
from tkinter import messagebox, ttk
import re
import webbrowser
from datetime import datetime

# ── Monoisotopic atomic masses ────────────────────────────────────────────────
ELEMENTS = {
    'C': 12.00000000000,
    'H':  1.00782503223,
    'O': 15.99491461957,
    'N': 14.00307400443,
    'P': 30.97376199842,
}

# ── Adduct mass offsets ───────────────────────────────────────────────────────
ADDUCTS = {
    '[M+H]+'      :  ELEMENTS['H'],
    '[M+Na]+'     :  22.98976928,   # NIST AME2020: Na-23 monoisotopic mass
    '[M+K]+'      :  38.96370649,   # NIST AME2020: K-39 monoisotopic mass
    '[M+NH4]+'    :  ELEMENTS['N'] + 4 * ELEMENTS['H'],
    '[M+2H]2+'    :  2 * ELEMENTS['H'],
    '[M-H]-'      : -ELEMENTS['H'],
    '[M+HCOO]-'   :  ELEMENTS['H'] + ELEMENTS['C'] + 2 * ELEMENTS['O'],
    '[M+CH3COO]-' :  2 * ELEMENTS['H'] + 2 * ELEMENTS['C'] + 2 * ELEMENTS['O'],
    '[M+Cl]-'     :  34.96885271,   # NIST AME2020: Cl-35 monoisotopic mass
    '[M-2H]2-'    : -2 * ELEMENTS['H'],
}

# ── Lipid class definitions ───────────────────────────────────────────────────
HEADGROUPS = {
    'PC': {'name': 'Phosphatidylcholine'},
    'PE': {'name': 'Phosphatidylethanolamine'},
    'PS': {'name': 'Phosphatidylserine'},
    'PG': {'name': 'Phosphatidylglycerol'},
    'PI': {'name': 'Phosphatidylinositol'},
    'PA': {'name': 'Phosphatidic Acid'},
    'SM': {'name': 'Sphingomyelin'},
    'CL': {'name': 'Cardiolipin'},
    'TG': {'name': 'Triacylglycerol'},
    'DG': {'name': 'Diacylglycerol'},
    'MG': {'name': 'Monoacylglycerol'},
    'FA': {'name': 'Fatty Acid'},
}

GP_CLASSES = {'PC', 'PE', 'PS', 'PG', 'PI', 'PA'}

# Glycerophospholipid cores (glycerol + phosphate + headgroup, two free sn-OH)
GP_CORES = {
    'PC': {'C': 8, 'H': 20, 'N': 1, 'O': 6, 'P': 1},
    'PE': {'C': 5, 'H': 14, 'N': 1, 'O': 6, 'P': 1},
    'PS': {'C': 6, 'H': 14, 'N': 1, 'O': 8, 'P': 1},
    'PG': {'C': 6, 'H': 15,           'O': 8, 'P': 1},
    'PI': {'C': 9, 'H': 19,           'O': 11,'P': 1},
    'PA': {'C': 3, 'H':  9,           'O': 6, 'P': 1},
}

GLYCEROL            = {'C': 3, 'H':  8, 'O':  3}
CL_CORE             = {'C': 9, 'H': 22, 'O': 13, 'P': 2}
PHOSPHOCHOLINE_ACID = {'C': 5, 'H': 14, 'N':  1, 'O': 4, 'P': 1}

LCB_PREFIXES = {
    'D': {'H_OFFSET': 3, 'O_COUNT': 2},   # dihydroxy (e.g. d18:1)
    'T': {'H_OFFSET': 3, 'O_COUNT': 3},   # trihydroxy
    'M': {'H_OFFSET': 5, 'O_COUNT': 1},   # monohydroxy
}

# ── Chemistry helpers ─────────────────────────────────────────────────────────
def combine_formulas(*parts):
    out = {}
    for part in parts:
        for elem, count in part.items():
            out[elem] = out.get(elem, 0) + count
    return {k: v for k, v in out.items() if v > 0}

def subtract_water(formula, n=1):
    out = dict(formula)
    out['H'] = out.get('H', 0) - 2 * n
    out['O'] = out.get('O', 0) - n
    return {k: v for k, v in out.items() if v > 0}

def calculate_formula_mass(formula_dict):
    return sum(ELEMENTS[e] * c for e, c in formula_dict.items() if e in ELEMENTS)

def fatty_acid_formula(carbons, double_bonds):
    return {'C': carbons, 'H': 2 * carbons - 2 * double_bonds, 'O': 2}

def fatty_alcohol_formula(carbons, double_bonds):
    return {'C': carbons, 'H': 2 * carbons + 2 - 2 * double_bonds, 'O': 1}

def alkenyl_alcohol_formula(carbons, double_bonds):
    return {'C': carbons, 'H': 2 * carbons - 2 * double_bonds, 'O': 1}

def lcb_formula(prefix, carbons, double_bonds):
    prefix = (prefix or 'D').upper()
    if prefix not in LCB_PREFIXES:
        raise ValueError(
            f"Unsupported sphingoid base prefix: '{prefix}'. Use d-, t-, or m- notation.")
    info = LCB_PREFIXES[prefix]
    return {
        'C': carbons,
        'H': 2 * carbons - 2 * double_bonds + info['H_OFFSET'],
        'N': 1,
        'O': info['O_COUNT'],
    }

# ── Class-specific builders ───────────────────────────────────────────────────
def build_gp_formula(base_head, linkage, chains, lyso_prefix):
    core = GP_CORES[base_head]

    # Lyso glycerophospholipids (one radyl chain)
    if lyso_prefix:
        if len(chains) != 1:
            raise ValueError(f"L{base_head} expects exactly 1 chain, got {len(chains)}.")
        _, c1, d1 = chains[0]
        if linkage == 'ester':
            radyl = fatty_acid_formula(c1, d1)
        elif linkage == 'O':
            radyl = fatty_alcohol_formula(c1, d1)
        elif linkage == 'P':
            radyl = alkenyl_alcohol_formula(c1, d1)
        else:
            raise ValueError(f"Unsupported linkage: {linkage}")
        return subtract_water(combine_formulas(core, radyl), 1)

    # Summed composition: PC(38:6), PC(O-38:6), PE(P-38:6)
    if len(chains) == 1:
        _, total_c, total_db = chains[0]
        if linkage == 'ester':
            radyls = {'C': total_c, 'H': 2 * total_c - 2 * total_db, 'O': 4}
        elif linkage == 'O':
            radyls = {'C': total_c, 'H': 2 * total_c + 2 - 2 * total_db, 'O': 3}
        elif linkage == 'P':
            radyls = {'C': total_c, 'H': 2 * total_c - 2 * total_db, 'O': 3}
        else:
            raise ValueError(f"Unsupported linkage: {linkage}")
        return subtract_water(combine_formulas(core, radyls), 2)

    # Explicit 2-chain GP lipids
    if len(chains) == 2:
        (_, c1, d1), (_, c2, d2) = chains
        if linkage == 'ester':
            parts = [core, fatty_acid_formula(c1, d1), fatty_acid_formula(c2, d2)]
        elif linkage == 'O':
            parts = [core, fatty_alcohol_formula(c1, d1), fatty_acid_formula(c2, d2)]
        elif linkage == 'P':
            parts = [core, alkenyl_alcohol_formula(c1, d1), fatty_acid_formula(c2, d2)]
        else:
            raise ValueError(f"Unsupported linkage: {linkage}")
        return subtract_water(combine_formulas(*parts), 2)

    raise ValueError(f"{base_head} expects 2 chains or 1 summed composition, got {len(chains)}.")

def build_glycerolipid_formula(base_head, chains):
    expected = {'MG': 1, 'DG': 2, 'TG': 3}[base_head]
    if len(chains) not in (1, expected):
        raise ValueError(
            f"{base_head} expects {expected} chain(s) or 1 summed composition, got {len(chains)}.")
    if len(chains) == 1:
        _, total_c, total_db = chains[0]
        acyl_pool = {'C': total_c, 'H': 2 * total_c - 2 * total_db, 'O': 2 * expected}
        return subtract_water(combine_formulas(GLYCEROL, acyl_pool), expected)
    parts = [GLYCEROL] + [fatty_acid_formula(c, db) for _, c, db in chains]
    return subtract_water(combine_formulas(*parts), expected)

def build_cl_formula(chains):
    if len(chains) not in (1, 4):
        raise ValueError(f"CL expects 4 chains or 1 summed composition, got {len(chains)}.")
    if len(chains) == 1:
        _, total_c, total_db = chains[0]
        acyl_pool = {'C': total_c, 'H': 2 * total_c - 2 * total_db, 'O': 8}
        return subtract_water(combine_formulas(CL_CORE, acyl_pool), 4)
    parts = [CL_CORE] + [fatty_acid_formula(c, db) for _, c, db in chains]
    return subtract_water(combine_formulas(*parts), 4)

def build_sm_formula(chains):
    if len(chains) != 2:
        raise ValueError(f"SM expects 2 chains (LCB/acyl), got {len(chains)}.")
    prefix, lcb_c, lcb_db = chains[0]
    _,      acyl_c, acyl_db = chains[1]
    lcb  = lcb_formula(prefix or 'D', lcb_c, lcb_db)
    acyl = fatty_acid_formula(acyl_c, acyl_db)
    return subtract_water(combine_formulas(lcb, acyl, PHOSPHOCHOLINE_ACID), 2)

# ── Main parser ───────────────────────────────────────────────────────────────
def parse_lipid_name(raw_name):
    """
    Supported notations:
      GP lipids     : PE(16:0/22:6), PC(O-16:0/18:1), PE(P-16:0/20:4)
      Lyso GP       : LPC(16:0), LPC(O-16:0), LPE(P-18:0)
      Sphingolipid  : SM(d18:1/16:0)
      Glycerolipids : TG(16:0/18:1/18:1), DG(16:0/18:1), MG(16:0)
      Cardiolipin   : CL(18:1/18:1/18:1/18:1)
      Fatty acid    : FA(16:0), FA(20:4)
      Sum comp.     : PC(38:4), PC(O-38:6), PE(P-38:6), TG(52:3), CL(72:8)
    """
    name = re.sub(r'\s+', '', raw_name).upper()
    lyso_prefix = bool(re.match(r'^L[A-Z]', name))
    if lyso_prefix:
        name = name[1:]

    match = re.match(r'^([A-Z]{1,3})\((.*)\)$', name)
    if not match:
        raise ValueError("Could not identify lipid class/headgroup.")

    raw_head, inner = match.groups()
    linkage = 'ester'

    if raw_head.endswith('O') and raw_head[:-1] in HEADGROUPS:
        base_head = raw_head[:-1]
        linkage   = 'O'
    elif raw_head in HEADGROUPS:
        base_head = raw_head
    else:
        raise ValueError(
            f"Unsupported lipid class: '{raw_head}'.\n"
            f"Supported: {', '.join(sorted(HEADGROUPS.keys()))}")

    if inner.startswith('O-'):
        linkage = 'O';  inner = inner[2:]
    elif inner.startswith('P-'):
        linkage = 'P';  inner = inner[2:]

    chain_matches = re.findall(r'([DMT]?)(\d+):(\d+)', inner)
    if not chain_matches:
        raise ValueError("No chain info found. Use format like 16:0/22:6 or 38:4.")
    chains = [(prefix.upper(), int(c), int(d)) for prefix, c, d in chain_matches]

    if base_head in GP_CLASSES:
        if lyso_prefix and len(chains) != 1:
            raise ValueError(f"L{base_head} expects 1 chain, got {len(chains)}.")
        if (not lyso_prefix) and len(chains) not in (1, 2):
            raise ValueError(f"{base_head} expects 2 chains or 1 summed composition.")
        formula = build_gp_formula(base_head, linkage, chains, lyso_prefix)
    elif base_head == 'SM':
        formula = build_sm_formula(chains);  linkage = 'amide'
    elif base_head in {'TG', 'DG', 'MG'}:
        formula = build_glycerolipid_formula(base_head, chains)
    elif base_head == 'CL':
        formula = build_cl_formula(chains)
    elif base_head == 'FA':
        if len(chains) != 1:
            raise ValueError(f"FA expects 1 chain, got {len(chains)}.")
        _, c1, d1 = chains[0]
        formula = fatty_acid_formula(c1, d1)
    else:
        raise ValueError(f"No formula builder implemented for class {base_head}.")

    return base_head, linkage, formula

def formula_to_hill(formula):
    hill = ['C', 'H'] + sorted(k for k in formula if k not in ('C', 'H'))
    return "".join(f"{e}{formula[e]}" for e in hill if e in formula)

# ── GUI constants ─────────────────────────────────────────────────────────────
LINK_LABELS = {
    'ester': 'Diacyl Ester',
    'O'    : 'Alkyl Ether (O-)',
    'P'    : 'Alkenyl Ether / Plasmalogen (P-)',
    'amide': 'Amide-linked Sphingolipid',
}

ADDUCT_ROWS = [
    ('[M+H]+'      , 'Positive', 1, 'pos' ),
    ('[M+Na]+'     , 'Positive', 1, 'pos' ),
    ('[M+K]+'      , 'Positive', 1, 'pos' ),
    ('[M+NH4]+'    , 'Positive', 1, 'pos' ),
    ('[M+2H]2+'    , 'Positive', 2, 'pos2'),
    ('[M-H]-'      , 'Negative', 1, 'neg' ),
    ('[M+HCOO]-'   , 'Negative', 1, 'neg' ),
    ('[M+CH3COO]-' , 'Negative', 1, 'neg' ),
    ('[M+Cl]-'     , 'Negative', 1, 'neg' ),
    ('[M-2H]2-'    , 'Negative', 2, 'neg2'),
]

# ── Main application ──────────────────────────────────────────────────────────
class LipidMassCalculatorGUI:
    def __init__(self, root):
        self.root = root
        self.root.title("Lipid Exact Mass Calculator v0.7 by Eylan Yutuc")
        self.root.geometry("920x620")
        self.root.minsize(760, 500)
        self.root.configure(bg='#f0f0f0')
        self.history = []
        self._current_formula_str = ""
        self._current_mass = None
        self._build_ui()

    def _build_ui(self):
        root = self.root

        # Title bar
        top = tk.Frame(root, bg='#2c3e50', pady=10)
        top.pack(fill='x')
        tk.Label(top, text="Lipid Exact Mass Calculator v0.7",
                 font=('Arial', 15, 'bold'), fg='white', bg='#2c3e50').pack()

        # Input row
        inp = tk.Frame(root, bg='#f0f0f0', pady=6)
        inp.pack(fill='x', padx=16)
        tk.Label(inp, text="Lipid Name:", font=('Arial', 10, 'bold'),
                 bg='#f0f0f0').pack(side='left', padx=(0, 6))
        self.entry_name = tk.Entry(inp, font=('Arial', 12), width=28,
                                   justify='center', relief='solid', bd=1)
        self.entry_name.pack(side='left')
        self.entry_name.insert(0, "PE(16:0/22:6)")
        self.entry_name.bind('<Return>', lambda e: self.on_calculate())

        tk.Button(inp, text="Calculate", font=('Arial', 10, 'bold'),
                  bg='#27ae60', fg='white', relief='flat', padx=10, pady=3,
                  cursor='hand2', command=self.on_calculate).pack(side='left', padx=8)

        self.btn_pubchem = tk.Button(inp, text="\U0001f50d Check PubChem",
                                     font=('Arial', 10, 'bold'),
                                     bg='#2980b9', fg='white', relief='flat',
                                     padx=10, pady=3, cursor='hand2',
                                     state='disabled',
                                     command=self.on_pubchem_lookup)
        self.btn_pubchem.pack(side='left', padx=(0, 8))

        tk.Button(inp, text="Clear History", font=('Arial', 9),
                  bg='#c0392b', fg='white', relief='flat', padx=8, pady=3,
                  cursor='hand2', command=self.clear_history).pack(side='left')

        tk.Label(inp,
                 text="e.g.  PE(16:0/22:6)  PC(O-16:0/18:1)  SM(d18:1/16:0)"
                      "  TG(16:0/18:1/18:1)  LPC(16:0)  CL(18:1/18:1/18:1/18:1)  FA(20:4)",
                 font=('Arial', 8), fg='#777', bg='#f0f0f0').pack(side='left', padx=10)

        # Body: left panel + separator + right history panel
        body = tk.Frame(root, bg='#f0f0f0')
        body.pack(fill='both', expand=True, padx=12, pady=(0, 8))

        # LEFT — results
        left = tk.Frame(body, bg='#f0f0f0')
        left.pack(side='left', fill='both', expand=True)

        res_frame = tk.LabelFrame(left, text=" Results ",
                                  font=('Arial', 9, 'bold'),
                                  bg='#f0f0f0', padx=10, pady=6)
        res_frame.pack(fill='x', pady=(4, 4))

        self.lbl_class = tk.Label(res_frame, text="Lipid Class:  --",
                                  font=('Arial', 10), bg='#f0f0f0', anchor='w')
        self.lbl_class.grid(row=0, column=0, sticky='w', pady=1)

        self.lbl_formula = tk.Label(res_frame, text="Molecular Formula:  --",
                                    font=('Arial', 10), bg='#f0f0f0', anchor='w')
        self.lbl_formula.grid(row=1, column=0, sticky='w', pady=1)

        self.lbl_mass = tk.Label(res_frame, text="Exact Mass (M):  -- Da",
                                 font=('Arial', 11, 'bold'), fg='#1a5276',
                                 bg='#f0f0f0', anchor='w')
        self.lbl_mass.grid(row=2, column=0, sticky='w', pady=3)

        adduct_frame = tk.LabelFrame(left, text=" Adduct m/z Values ",
                                     font=('Arial', 9, 'bold'),
                                     bg='#f0f0f0', padx=6, pady=6)
        adduct_frame.pack(fill='both', expand=True, pady=(0, 4))

        style = ttk.Style()
        style.configure('Lipid.Treeview',         font=('Courier', 9), rowheight=22)
        style.configure('Lipid.Treeview.Heading', font=('Arial', 9, 'bold'))
        style.configure('Hist.Treeview',          font=('Arial', 9),   rowheight=20)
        style.configure('Hist.Treeview.Heading',  font=('Arial', 8,    'bold'))

        cols = ('Adduct', 'm/z (Da)', 'Mode', 'z')
        self.adduct_tree = ttk.Treeview(adduct_frame, columns=cols,
                                        show='headings', height=11,
                                        style='Lipid.Treeview')
        for c, w in zip(cols, [168, 148, 80, 36]):
            self.adduct_tree.heading(c, text=c)
            self.adduct_tree.column(c, width=w, anchor='center', stretch=False)

        self.adduct_tree.tag_configure('pos',  background='#d5f5e3')
        self.adduct_tree.tag_configure('pos2', background='#a9dfbf')
        self.adduct_tree.tag_configure('neg',  background='#d6eaf8')
        self.adduct_tree.tag_configure('neg2', background='#a9cce3')

        vsb_a = ttk.Scrollbar(adduct_frame, orient='vertical',
                               command=self.adduct_tree.yview)
        self.adduct_tree.configure(yscrollcommand=vsb_a.set)
        self.adduct_tree.pack(side='left', fill='both', expand=True)
        vsb_a.pack(side='right', fill='y')

        tk.Frame(body, bg='#cccccc', width=1).pack(side='left', fill='y', padx=6)

        # RIGHT — history
        right = tk.Frame(body, bg='#f0f0f0', width=230)
        right.pack(side='left', fill='y')
        right.pack_propagate(False)

        tk.Label(right, text="History", font=('Arial', 10, 'bold'),
                 bg='#2c3e50', fg='white', pady=5).pack(fill='x')

        hcols = ('Input', 'M (Da)')
        self.hist_tree = ttk.Treeview(right, columns=hcols,
                                      show='headings', style='Hist.Treeview')
        self.hist_tree.heading('Input',  text='Input')
        self.hist_tree.heading('M (Da)', text='M (Da)')
        self.hist_tree.column('Input',   width=138, anchor='w',      stretch=True)
        self.hist_tree.column('M (Da)',  width=88,  anchor='center', stretch=False)

        self.hist_tree.tag_configure('even', background='#ffffff')
        self.hist_tree.tag_configure('odd',  background='#eaf4fb')

        vsb_h = ttk.Scrollbar(right, orient='vertical',
                               command=self.hist_tree.yview)
        self.hist_tree.configure(yscrollcommand=vsb_h.set)
        self.hist_tree.pack(side='left', fill='both', expand=True)
        vsb_h.pack(side='right', fill='y')

        self.hist_tree.bind('<<TreeviewSelect>>', self.on_history_select)

        self.status = tk.Label(
            root,
            text="Ready — enter a lipid name and press Calculate or Enter.",
            font=('Arial', 8), fg='#555', bg='#dde1e6',
            anchor='w', padx=8, pady=3)
        self.status.pack(fill='x', side='bottom')

    # ── Event handlers ────────────────────────────────────────────────────────
    def on_calculate(self):
        raw = self.entry_name.get().strip()
        if not raw:
            messagebox.showwarning("Input Error", "Please enter a lipid name.")
            return
        try:
            base_head, linkage, formula = parse_lipid_name(raw)
            M = calculate_formula_mass(formula)
            self._display_result(raw, base_head, linkage, formula, M)
            self._add_to_history(raw, base_head, linkage, formula, M)
            self.status.config(
                text=f"OK  |  '{raw}'  →  M = {M:.6f} Da  |  "
                     f"{datetime.now().strftime('%H:%M:%S')}")
        except Exception as e:
            messagebox.showerror("Calculation Error", str(e))
            self.status.config(text=f"Error: {e}")

    def on_pubchem_lookup(self):
        """Search PubChem by molecular formula in the default web browser."""
        if not self._current_formula_str:
            return
        url = (
            "https://pubchem.ncbi.nlm.nih.gov/#query="
            + self._current_formula_str
            + "&input_type=formula"
        )
        try:
            webbrowser.open(url)
            self.status.config(
                text=f"Opened PubChem search for {self._current_formula_str}")
        except Exception as e:
            messagebox.showerror("Browser Error", f"Could not open browser:\n{e}")

    def _display_result(self, raw, base_head, linkage, formula, M):
        formula_str = formula_to_hill(formula)
        class_label = HEADGROUPS[base_head]['name']
        if base_head in GP_CLASSES or base_head == 'SM':
            class_label += f"  [{LINK_LABELS.get(linkage, linkage)}]"

        self._current_formula_str = formula_str
        self._current_mass = M
        self.btn_pubchem.config(state='normal')

        self.lbl_class.config(text=f"Lipid Class:  {class_label}")
        self.lbl_formula.config(text=f"Molecular Formula:  {formula_str}")
        self.lbl_mass.config(text=f"Exact Mass (M):  {M:.6f} Da")

        for row in self.adduct_tree.get_children():
            self.adduct_tree.delete(row)
        for adduct, mode, z, tag in ADDUCT_ROWS:
            mz = (M + ADDUCTS[adduct]) / z
            z_label = f"{z}+" if mode == 'Positive' else f"{z}-"
            self.adduct_tree.insert('', 'end',
                                    values=(adduct, f"{mz:.6f}", mode, z_label),
                                    tags=(tag,))

    def _add_to_history(self, raw, base_head, linkage, formula, M):
        if self.history and self.history[-1]['input'].upper() == raw.upper():
            return
        record = {'input': raw, 'base_head': base_head,
                  'linkage': linkage, 'formula': formula, 'mass': M}
        self.history.append(record)
        tag = 'even' if len(self.history) % 2 == 0 else 'odd'
        self.hist_tree.insert('', 0,
                              iid=str(len(self.history) - 1),
                              values=(raw, f"{M:.4f}"),
                              tags=(tag,))

    def on_history_select(self, event):
        sel = self.hist_tree.selection()
        if not sel:
            return
        idx = int(sel[0])
        r = self.history[idx]
        self.entry_name.delete(0, 'end')
        self.entry_name.insert(0, r['input'])
        self._display_result(r['input'], r['base_head'],
                             r['linkage'], r['formula'], r['mass'])
        self.status.config(
            text=f"Loaded from history  |  '{r['input']}'  →  M = {r['mass']:.6f} Da")

    def clear_history(self):
        if not self.history:
            return
        if messagebox.askyesno("Clear History", "Clear all calculation history?"):
            self.history.clear()
            for row in self.hist_tree.get_children():
                self.hist_tree.delete(row)
            self.status.config(text="History cleared.")


if __name__ == "__main__":
    root = tk.Tk()
    app = LipidMassCalculatorGUI(root)
    root.mainloop()
