"""Lipid Exact Mass Calculator - browser-friendly core (runs in Pyodide).

No GUI dependencies. The web page calls `calculate_json(name)`.
"""
import json
import re

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
    'M': {'H_OFFSET': 5, 'O_COUNT': 1},
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

    return base_head, linkage, formula


def formula_to_hill(formula):
    hill = ['C', 'H'] + sorted(k for k in formula if k not in ('C', 'H'))
    return "".join(f"{e}{formula[e]}" for e in hill if e in formula)


def calculate(raw_name):
    """Return a JSON-serialisable dict; never raises."""
    try:
        base_head, linkage, formula = parse_lipid_name(raw_name)
        mass = calculate_formula_mass(formula)
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
