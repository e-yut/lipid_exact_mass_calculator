"""Lipid Exact Mass Calculator - browser-friendly core (runs in Pyodide).

No GUI dependencies. The web page calls `calculate_json(name)`.
"""
import json
import re

__version__ = "0.8"

ELEMENTS = {
    'C': 12.00000000000,
    'H':  1.00782503223,
    'O': 15.99491461957,
    'N': 14.00307400443,
    'P': 30.97376199842,
}

ADDUCTS = {
    '[M+H]+'      :  ELEMENTS['H'],
    '[M+Na]+'     :  22.98976928,
    '[M+K]+'      :  38.96370649,
    '[M+NH4]+'    :  ELEMENTS['N'] + 4 * ELEMENTS['H'],
    '[M+2H]2+'    :  2 * ELEMENTS['H'],
    '[M-H]-'      : -ELEMENTS['H'],
    '[M+HCOO]-'   :  ELEMENTS['H'] + ELEMENTS['C'] + 2 * ELEMENTS['O'],
    '[M+CH3COO]-' :  2 * ELEMENTS['H'] + 2 * ELEMENTS['C'] + 2 * ELEMENTS['O'],
    '[M+Cl]-'     :  34.96885271,
    '[M-2H]2-'    : -2 * ELEMENTS['H'],
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

HEADGROUPS = {
    'PC': 'Phosphatidylcholine',
    'PE': 'Phosphatidylethanolamine',
    'PS': 'Phosphatidylserine',
    'PG': 'Phosphatidylglycerol',
    'PI': 'Phosphatidylinositol',
    'PA': 'Phosphatidic Acid',
    'SM': 'Sphingomyelin',
    'CL': 'Cardiolipin',
    'TG': 'Triacylglycerol',
    'DG': 'Diacylglycerol',
    'MG': 'Monoacylglycerol',
    'FA': 'Fatty Acid',
}

GP_CLASSES = {'PC', 'PE', 'PS', 'PG', 'PI', 'PA'}

GP_CORES = {
    'PC': {'C': 8, 'H': 20, 'N': 1, 'O': 6, 'P': 1},
    'PE': {'C': 5, 'H': 14, 'N': 1, 'O': 6, 'P': 1},
    'PS': {'C': 6, 'H': 14, 'N': 1, 'O': 8, 'P': 1},
    'PG': {'C': 6, 'H': 15,           'O': 8, 'P': 1},
    'PI': {'C': 9, 'H': 19,           'O': 11, 'P': 1},
    'PA': {'C': 3, 'H':  9,           'O': 6, 'P': 1},
}

GLYCEROL            = {'C': 3, 'H':  8, 'O':  3}
CL_CORE             = {'C': 9, 'H': 22, 'O': 13, 'P': 2}
PHOSPHOCHOLINE_ACID = {'C': 5, 'H': 14, 'N':  1, 'O': 4, 'P': 1}

LCB_PREFIXES = {
    'D': {'H_OFFSET': 3, 'O_COUNT': 2},
    'T': {'H_OFFSET': 3, 'O_COUNT': 3},
    'M': {'H_OFFSET': 3, 'O_COUNT': 1},
}

LINK_LABELS = {
    'ester': 'Diacyl Ester',
    'O'    : 'Alkyl Ether (O-)',
    'P'    : 'Alkenyl Ether / Plasmalogen (P-)',
    'amide': 'Amide-linked Sphingolipid',
}


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


def calculate_formula_mass(formula):
    return sum(ELEMENTS[e] * c for e, c in formula.items() if e in ELEMENTS)


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


def build_gp_formula(base_head, linkage, chains, lyso_prefix):
    core = GP_CORES[base_head]

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
    _, acyl_c, acyl_db = chains[1]
    lcb = lcb_formula(prefix or 'D', lcb_c, lcb_db)
    acyl = fatty_acid_formula(acyl_c, acyl_db)
    return subtract_water(combine_formulas(lcb, acyl, PHOSPHOCHOLINE_ACID), 2)


def parse_lipid_name(raw_name):
    """Parse a lipid shorthand name -> (class, linkage, formula dict)."""
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
        linkage = 'O'
    elif raw_head in HEADGROUPS:
        base_head = raw_head
    else:
        raise ValueError(
            f"Unsupported lipid class: '{raw_head}'. "
            f"Supported: {', '.join(sorted(HEADGROUPS))}")

    if inner.startswith('O-'):
        linkage = 'O'
        inner = inner[2:]
    elif inner.startswith('P-'):
        linkage = 'P'
        inner = inner[2:]

    chain_matches = re.findall(r'([DMT]?)(\d+):(\d+)', inner)
    if not chain_matches:
        raise ValueError("No chain info found. Use format like 16:0/22:6 or 38:4.")
    chains = [(p.upper(), int(c), int(d)) for p, c, d in chain_matches]

    if base_head in GP_CLASSES:
        if lyso_prefix and len(chains) != 1:
            raise ValueError(f"L{base_head} expects 1 chain, got {len(chains)}.")
        if (not lyso_prefix) and len(chains) not in (1, 2):
            raise ValueError(f"{base_head} expects 2 chains or 1 summed composition.")
        formula = build_gp_formula(base_head, linkage, chains, lyso_prefix)
    elif lyso_prefix:
        raise ValueError(f"Lyso (L) prefix is only supported for glycerophospholipids, not {base_head}.")
    elif base_head == 'SM':
        formula = build_sm_formula(chains)
        linkage = 'amide'
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

    for _, c, d in chains:
        if c <= 0 or d < 0 or d > c:
            raise ValueError(f"Invalid chain {c}:{d}.")

    return base_head, linkage, formula, chains, lyso_prefix


# ── Structure (SMILES) generation ─────────────────────────────────────────────
# Shorthand names do not give double-bond positions, so common positions are
# assumed (e.g. 20:4 -> Δ5,8,11,14); otherwise methylene-interrupted placement.
COMMON_DB_POSITIONS = {
    (14, 1): [9], (16, 1): [9], (18, 1): [9], (18, 2): [9, 12],
    (18, 3): [9, 12, 15], (20, 1): [11], (20, 2): [11, 14],
    (20, 3): [8, 11, 14], (20, 4): [5, 8, 11, 14], (20, 5): [5, 8, 11, 14, 17],
    (22, 1): [13], (22, 4): [7, 10, 13, 16], (22, 5): [7, 10, 13, 16, 19],
    (22, 6): [4, 7, 10, 13, 16, 19], (24, 1): [15],
}

PHOSPHOCHOLINE = 'P(=O)([O-])OCC[N+](C)(C)C'
GP_HEAD_SMILES = {
    'PC': PHOSPHOCHOLINE,
    'PE': 'P(=O)(O)OCCN',
    'PS': 'P(=O)(O)OCC(N)C(=O)O',
    'PG': 'P(=O)(O)OCC(O)CO',
    'PI': 'P(=O)(O)OC1C(O)C(O)C(O)C(O)C1O',
    'PA': 'P(=O)(O)O',
}


def db_positions(carbons, double_bonds, min_start=2, common=True):
    if double_bonds == 0:
        return []
    if common and (carbons, double_bonds) in COMMON_DB_POSITIONS:
        pos = COMMON_DB_POSITIONS[(carbons, double_bonds)]
        if pos[0] >= min_start and pos[-1] < carbons:
            return list(pos)
    for last in (carbons - 6, carbons - 3, carbons - 1):
        start = last - 3 * (double_bonds - 1)
        if start >= min_start:
            return list(range(start, last + 1, 3))
    pos = list(range(min_start, min_start + 2 * double_bonds, 2))
    return pos if pos[-1] < carbons else None


def chain_smiles(carbons, dbs, subs=None):
    """Linear chain C1..Cn; `dbs` are Δ positions; `subs` maps index->branches."""
    subs = subs or {}
    out = []
    for i in range(1, carbons + 1):
        out.append('C' + ''.join(f'({s})' for s in subs.get(i, [])))
        if i < carbons and i in dbs:
            out.append('=')
    return ''.join(out)


def acyl_smiles(c, d):
    dbs = db_positions(c, d)
    if dbs is None or c < 1:
        return None
    return chain_smiles(c, dbs, {1: ['=O']})


def radyl_smiles(c, d, linkage):
    if linkage == 'ester':
        return acyl_smiles(c, d)
    if linkage == 'O':
        dbs = db_positions(c, d, min_start=2)
        return None if dbs is None else chain_smiles(c, dbs)
    if linkage == 'P':
        if c < 2:
            return None
        dbs = db_positions(c, d, min_start=3)
        return None if dbs is None else chain_smiles(c, [1] + dbs)
    return None


def lcb_smiles(prefix, c, d, acyl):
    prefix = (prefix or 'D').upper()
    if c < 3:
        return None
    if prefix == 'D' and d >= 1:
        rest = db_positions(c, d - 1, min_start=6, common=False)
        dbs = None if rest is None else [4] + rest
    else:
        dbs = db_positions(c, d, min_start=4 if prefix == 'D' else 5 if prefix == 'T' else 3,
                           common=False)
    if dbs is None:
        return None
    subs = {1: ['O' + PHOSPHOCHOLINE], 2: ['N' + acyl]}
    if prefix in ('D', 'T'):
        subs[3] = ['O']
    if prefix == 'T':
        subs[4] = ['O']
    return chain_smiles(c, dbs, subs)


def lipid_smiles(base_head, linkage, chains, lyso):
    """Representative SMILES, or None if not drawable (e.g. summed composition)."""
    if base_head == 'FA':
        a = acyl_smiles(chains[0][1], chains[0][2])
        return None if a is None else 'O' + a

    if base_head == 'SM':
        (p, c1, d1), (_, c2, d2) = chains
        acyl = acyl_smiles(c2, d2)
        return None if acyl is None else lcb_smiles(p, c1, d1, acyl)

    if base_head in GP_CLASSES:
        if lyso:
            r1 = radyl_smiles(chains[0][1], chains[0][2], linkage)
            return None if r1 is None else f'C(O{r1})C(O)CO{GP_HEAD_SMILES[base_head]}'
        if len(chains) != 2:
            return None
        r1 = radyl_smiles(chains[0][1], chains[0][2], linkage)
        r2 = acyl_smiles(chains[1][1], chains[1][2])
        if r1 is None or r2 is None:
            return None
        return f'C(O{r1})C(O{r2})CO{GP_HEAD_SMILES[base_head]}'

    expected = {'MG': 1, 'DG': 2, 'TG': 3, 'CL': 4}.get(base_head)
    if expected is None or len(chains) != expected:
        return None
    acyls = [acyl_smiles(c, d) for _, c, d in chains]
    if any(a is None for a in acyls):
        return None
    if base_head == 'MG':
        return f'C(O{acyls[0]})C(O)CO'
    if base_head == 'DG':
        return f'C(O{acyls[0]})C(O{acyls[1]})CO'
    if base_head == 'TG':
        return f'C(O{acyls[0]})C(O{acyls[1]})CO{acyls[2]}'
    r1, r2, r3, r4 = acyls
    return (f'OC(COP(=O)(O)OCC(O{r2})CO{r1})'
            f'COP(=O)(O)OCC(O{r3})CO{r4}')


def formula_to_hill(formula):
    hill = ['C', 'H'] + sorted(k for k in formula if k not in ('C', 'H'))
    return "".join(f"{e}{formula[e]}" for e in hill if e in formula)


def calculate(raw_name):
    """Return a JSON-serialisable dict; never raises."""
    try:
        base_head, linkage, formula, chains, lyso = parse_lipid_name(raw_name)
        mass = calculate_formula_mass(formula)
        try:
            smiles = lipid_smiles(base_head, linkage, chains, lyso)
        except Exception:
            smiles = None
        label = HEADGROUPS[base_head]
        if base_head in GP_CLASSES or base_head == 'SM':
            label += f" [{LINK_LABELS.get(linkage, linkage)}]"
        formula_str = formula_to_hill(formula)
        return {
            'ok': True,
            'input': raw_name.strip(),
            'class': label,
            'formula': formula_str,
            'mass': mass,
            'smiles': smiles,
            'summed': smiles is None and len(chains) == 1 and base_head not in ('FA', 'MG')
                      and not lyso,
            'adducts': [
                {'adduct': a, 'mode': m, 'z': z, 'tag': t,
                 'mz': (mass + ADDUCTS[a]) / z}
                for a, m, z, t in ADDUCT_ROWS
            ],
        }
    except Exception as exc:
        return {'ok': False, 'error': str(exc)}


def calculate_json(raw_name):
    return json.dumps(calculate(raw_name))
