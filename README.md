# Lipid Exact Mass Calculator

[![Python 3.8+](https://img.shields.io/badge/python-3.8%2B-blue.svg)](https://www.python.org/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)
[![Platform](https://img.shields.io/badge/platform-Windows%20%7C%20macOS%20%7C%20Linux-lightgrey.svg)]()
[![No dependencies](https://img.shields.io/badge/dependencies-none-brightgreen.svg)]()

A desktop GUI tool for calculating **monoisotopic exact masses** and **adduct m/z values** for lipid species from their shorthand nomenclature. Built with Python and Tkinter — no third-party packages required.

Developed by **Eylan Yutuc** for use in lipidomics and LC-MS/MS research workflows.

---

## Features

- Computes **exact neutral mass (M)** from first-principles molecular formula construction — no lookup tables, no hardcoded masses
- Calculates **10 adduct m/z values** simultaneously (positive and negative ion modes)
- Supports **12 lipid classes** across glycerophospholipids, sphingolipids, glycerolipids, and fatty acids
- Accepts both **summed composition** (e.g. `PC(38:4)`) and **sn-resolved** notation (e.g. `PC(16:0/22:6)`)
- **Lyso variants** supported for all glycerophospholipid classes
- **Ether linkage notation** — diacyl ester, alkyl ether (O-), and alkenyl ether / plasmalogen (P-)
- **🔍 Check PubChem** button — opens a PubChem molecular formula search in your browser after each calculation to cross-validate results
- **Calculation history panel** — click any previous entry to instantly reload it
- All atomic masses and ionic constants sourced from **NIST AME2020** high-precision values

---

## Requirements

- **Python 3.8 or higher**
- **No third-party packages** — only standard library modules are used:
  - `tkinter` (GUI)
  - `re` (regex parsing)
  - `webbrowser` (PubChem lookup)
  - `datetime` (status bar timestamps)

> **Note for Linux users:** `tkinter` is not always bundled with the system Python. If you get a `ModuleNotFoundError`, install it with:
> ```bash
> sudo apt install python3-tk        # Debian / Ubuntu
> sudo dnf install python3-tkinter   # Fedora / RHEL
> ```

---

## Installation & Usage

No installation is required. Download or clone the repository, then run the script directly:

```bash
# Clone
git clone https://github.com/YOUR_USERNAME/lipid-exact-mass-calculator.git
cd lipid-exact-mass-calculator

# Run
python exactmass_calc_v0.7.py
```

On **Windows**, you can double-click the `.py` file if Python is associated with `.py` files in your system settings.

---

## Supported Lipid Classes

### Glycerophospholipids

| Class | Full Name | Example Input |
|-------|-----------|---------------|
| `PC` | Phosphatidylcholine | `PC(16:0/22:6)` |
| `PE` | Phosphatidylethanolamine | `PE(18:0/20:4)` |
| `PS` | Phosphatidylserine | `PS(18:0/18:1)` |
| `PG` | Phosphatidylglycerol | `PG(16:0/18:1)` |
| `PI` | Phosphatidylinositol | `PI(18:0/20:4)` |
| `PA` | Phosphatidic Acid | `PA(16:0/18:1)` |

Lyso variants are supported with the `L` prefix: `LPC(16:0)`, `LPE(18:1)`, `LPS(18:0)`, `LPG(16:0)`, `LPI(18:0)`, `LPA(16:0)`.

### Sphingolipids

| Class | Full Name | Example Input |
|-------|-----------|---------------|
| `SM` | Sphingomyelin | `SM(d18:1/16:0)` |

### Glycerolipids

| Class | Full Name | Example Input |
|-------|-----------|---------------|
| `TG` | Triacylglycerol | `TG(16:0/18:1/18:1)` |
| `DG` | Diacylglycerol | `DG(16:0/18:1)` |
| `MG` | Monoacylglycerol | `MG(16:0)` |

### Other

| Class | Full Name | Example Input |
|-------|-----------|---------------|
| `CL` | Cardiolipin | `CL(18:1/18:1/18:1/18:1)` |
| `FA` | Fatty Acid | `FA(20:4)` |

---

## Notation Guide

### Chain Notation

Chains are specified as `C:DB` (carbons:double bonds), separated by `/` for sn-resolved species.

| Notation type | Example | Meaning |
|---|---|---|
| Summed composition | `PC(38:4)` | Total carbons and double bonds across both chains |
| sn-resolved | `PC(16:0/22:6)` | Explicit sn-1 and sn-2 chains |
| Lyso | `LPC(16:0)` | Single radyl chain |

### Glycerophospholipid Linkage Prefixes

| Prefix | Linkage | Example |
|--------|---------|---------|
| *(none)* | Diacyl ester | `PC(16:0/18:1)` |
| `O-` | Alkyl ether (plasmanyl) | `PC(O-16:0/18:1)` |
| `P-` | Alkenyl ether / plasmalogen (plasmenyl) | `PE(P-16:0/20:4)` |

Ether linkage prefixes work for lyso variants too: `LPC(O-16:0)`, `LPE(P-18:0)`.

### Sphingoid Base Prefixes (SM only)

| Prefix | Hydroxylation | Example |
|--------|---------------|---------|
| `d` | Dihydroxy (default) | `SM(d18:1/16:0)` |
| `t` | Trihydroxy | `SM(t18:0/24:1)` |
| `m` | Monohydroxy | `SM(m18:1/16:0)` |

---

## Adduct m/z Values Calculated

| Adduct | Ion mode | Charge (z) |
|--------|----------|------------|
| `[M+H]+` | Positive | 1+ |
| `[M+Na]+` | Positive | 1+ |
| `[M+K]+` | Positive | 1+ |
| `[M+NH4]+` | Positive | 1+ |
| `[M+2H]2+` | Positive | 2+ |
| `[M-H]-` | Negative | 1− |
| `[M+HCOO]-` | Negative | 1− |
| `[M+CH3COO]-` | Negative | 1− |
| `[M+Cl]-` | Negative | 1− |
| `[M-2H]2-` | Negative | 2− |

Results are colour-coded in the GUI: green shades for positive mode, blue shades for negative mode.

---

## Monoisotopic Constants (NIST AME2020)

All atomic masses and ionic adduct offsets use values from the NIST Atomic Mass Evaluation 2020.

| Isotope | Mass (Da) | Role |
|---------|-----------|------|
| ¹H | 1.00782503223 | Elemental mass + `[M+H]+`, `[M-H]-`, `[M+NH4]+`, formate, acetate |
| ¹²C | 12.00000000000 | Elemental mass |
| ¹⁴N | 14.00307400443 | Elemental mass |
| ¹⁶O | 15.99491461957 | Elemental mass |
| ³¹P | 30.97376199842 | Elemental mass |
| ²³Na | 22.98976928 | `[M+Na]+` adduct |
| ³⁹K | 38.96370649 | `[M+K]+` adduct |
| ³⁵Cl | 34.96885271 | `[M+Cl]-` adduct |

---

## Example Calculations

| Input | Formula | Exact Mass M (Da) |
|-------|---------|-------------------|
| `PE(16:0/22:6)` | C43H76NO8P | 765.530... |
| `PC(O-16:0/18:1)` | C42H82NO7P | 743.577... |
| `PE(P-16:0/20:4)` | C39H70NO7P | 699.489... |
| `SM(d18:1/16:0)` | C39H79N2O6P | 702.566... |
| `TG(16:0/18:1/18:1)` | C55H102O6 | 874.760... |
| `LPC(16:0)` | C24H50NO7P | 495.332... |
| `CL(18:1/18:1/18:1/18:1)` | C81H151O17P2 | 1466.01... |
| `FA(20:4)` | C20H32O2 | 304.240... |
| `PI(18:0/20:4)` | C47H83O13P | 886.556... |
| `PC(38:4)` | C46H82NO8P | 811.577... |

---

## How It Works

Formula construction is performed from first principles rather than a database lookup. For each lipid class, the calculator:

1. **Parses** the shorthand name using regex to extract the lipid class, linkage type, and chain composition
2. **Builds** the molecular formula by combining pre-defined backbone/headgroup core atoms with radyl chain atoms and removing condensation water molecules
3. **Calculates** the exact neutral mass by summing NIST monoisotopic atomic masses
4. **Derives** adduct m/z values as `(M + adduct_offset) / z`

This approach means the tool correctly handles structural variation in backbone chemistry (e.g. glycerophospholipid cores differ between PC, PE, PS, PG, PI, and PA) rather than applying a single universal formula offset.

---

## PubChem Cross-Validation

After any calculation, the **🔍 Check PubChem** button opens a PubChem formula search in your default browser:

```
https://pubchem.ncbi.nlm.nih.gov/#query=<FORMULA>&input_type=formula
```

This returns all PubChem compounds matching the molecular formula, ranked by exact mass similarity, allowing direct cross-validation of the computed result.

---

## Changelog

### v0.7 — 2026-06-01
- Updated Na-23, K-39, and Cl-35 adduct constants to full NIST AME2020 precision (8 significant figures)
- Added **🔍 Check PubChem** button — opens molecular formula search in default browser after each calculation
- Button is disabled until a successful calculation has been run in the current session

### v0.6 — 2026-06-01
- Complete rewrite of formula construction logic using explicit class-specific chemistry builders for all 12 lipid classes
- Fixed hydrogen undercount in glycerophospholipid diacyl species
- Added lyso GP support with ester, O-, and P-linkage handling
- Added summed composition parsing for all classes
- Added sphingoid base prefix support (d/t/m) for SM

### v0.5 — 2026-05-31
- Updated all atomic masses to NIST high-precision values
- Initial Tkinter GUI with colour-coded adduct table and history panel

### v0.4
- Initial release with basic GUI and formula calculation

---

## License

MIT License — see [LICENSE](LICENSE) for full terms.  
Copyright © 2026 Eylan Yutuc.
