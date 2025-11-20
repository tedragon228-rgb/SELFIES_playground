"""
Utility functions for SELFIES-based molecule exploration.
"""
from __future__ import annotations

import random
from typing import Dict, List, Optional, Set, Tuple

import numpy as np
import pandas as pd
import selfies as sf
from rdkit import Chem, DataStructs
from rdkit.Chem import AllChem, Descriptors, Draw
from PIL import Image

# --------------------------
# Utility functions
# --------------------------

def smiles_to_mol(smiles: str) -> Optional[Chem.Mol]:
    """Convert SMILES string to RDKit Mol, return None if invalid."""
    if smiles is None:
        return None
    smiles = smiles.strip()
    if not smiles:
        return None
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    try:
        Chem.SanitizeMol(mol)
    except Exception:
        return None
    return mol


def mol_to_smiles(mol: Chem.Mol) -> str:
    """Canonical SMILES from Mol."""
    return Chem.MolToSmiles(mol, canonical=True)


def smiles_list_to_grid_image(
    smiles_list: List[str],
    legends: Optional[List[str]] = None,
    mols_per_row: int = 4,
    sub_img_size: Tuple[int, int] = (250, 250),
) -> Image.Image:
    """Render multiple SMILES strings into a grid image."""
    mols = [smiles_to_mol(smi) for smi in smiles_list]
    mols = [m for m in mols if m is not None]
    if not mols:
        return Image.new("RGB", (sub_img_size[0], sub_img_size[1]), "white")
    if legends is None:
        legends = [mol_to_smiles(m) for m in mols]
    img = Draw.MolsToGridImage(
        mols,
        molsPerRow=mols_per_row,
        subImgSize=sub_img_size,
        legends=legends,
    )
    return img


def mol_fingerprint(mol: Chem.Mol, radius: int = 2, n_bits: int = 2048):
    return AllChem.GetMorganFingerprintAsBitVect(mol, radius, nBits=n_bits)


def tanimoto_sim(m1: Chem.Mol, m2: Chem.Mol) -> float:
    if m1 is None or m2 is None:
        return 0.0
    fp1 = mol_fingerprint(m1)
    fp2 = mol_fingerprint(m2)
    return float(DataStructs.TanimotoSimilarity(fp1, fp2))


def is_charged_smiles(smiles: str) -> bool:
    """Check whether a SMILES string contains ionic annotations."""
    if smiles is None:
        return False
    return ("+" in smiles) or ("-" in smiles)


# --------------------------
# SELFIES helpers
# --------------------------

def smiles_to_selfies(smiles: str) -> Optional[str]:
    """SMILES -> SELFIES (canonical SMILES first)."""
    try:
        mol = smiles_to_mol(smiles)
        if mol is None:
            return None
        can_smi = mol_to_smiles(mol)
        return sf.encoder(can_smi)
    except Exception:
        return None


def selfies_to_smiles(selfies: str) -> Optional[str]:
    """SELFIES -> SMILES."""
    try:
        smi = sf.decoder(selfies)
        mol = smiles_to_mol(smi)
        if mol is None:
            return None
        return mol_to_smiles(mol)
    except Exception:
        return None


def get_selfies_alphabet(seed_smiles: List[str]) -> Set[str]:
    """Build SELFIES alphabet from seed SMILES and robust defaults."""
    alphabet: Set[str] = set()
    for s in seed_smiles:
        sf_str = smiles_to_selfies(s)
        if sf_str is None:
            continue
        for sym in sf.split_selfies(sf_str):
            alphabet.add(sym)
    default_alphabet = sf.get_semantic_robust_alphabet()
    alphabet.update(default_alphabet)
    return alphabet


def get_default_selfies_alphabet() -> Set[str]:
    """Return only the default SELFIES alphabet."""
    return set(sf.get_semantic_robust_alphabet())


def random_selfies_mutation(
    selfies: str,
    alphabet: List[str],
    num_edits: int = 1,
    p_replace: float = 0.5,
    p_insert: float = 0.25,
    p_delete: float = 0.25,
) -> str:
    """Apply token-level replace/insert/delete mutations to a SELFIES string."""
    tokens = list(sf.split_selfies(selfies))

    for _ in range(num_edits):
        if len(tokens) == 0:
            op = "insert"
        else:
            op = random.choices(
                ["replace", "insert", "delete"],
                weights=[p_replace, p_insert, p_delete],
                k=1,
            )[0]

        if op == "replace":
            idx = random.randrange(len(tokens))
            tokens[idx] = random.choice(alphabet)
        elif op == "insert":
            idx = random.randrange(len(tokens) + 1)
            tokens.insert(idx, random.choice(alphabet))
        elif op == "delete":
            idx = random.randrange(len(tokens))
            del tokens[idx]

    return "".join(tokens)


# --------------------------
# 1. Local chemical space
# --------------------------

def generate_local_chemical_space(
    seed_smiles: List[str],
    alphabet: Set[str],
    n_mutations_per_seed: int = 256,
    min_edits: int = 1,
    max_edits: int = 3,
    max_total: int = 2000,
    random_seed: int = 0,
    allow_charged: bool = True,
) -> pd.DataFrame:
    """Generate molecules around seeds via SELFIES mutations."""
    random.seed(random_seed)
    np.random.seed(random_seed)

    alphabet_list = list(alphabet)

    records = []
    seen: Set[str] = set()

    for seed_idx, smi in enumerate(seed_smiles):
        seed_mol = smiles_to_mol(smi)
        if seed_mol is None:
            continue
        seed_sf = smiles_to_selfies(smi)
        if seed_sf is None:
            continue

        for _ in range(n_mutations_per_seed):
            n_edits = random.randint(min_edits, max_edits)
            mutated_sf = random_selfies_mutation(seed_sf, alphabet_list, num_edits=n_edits)
            mutated_smi = selfies_to_smiles(mutated_sf)
            if mutated_smi is None:
                continue

            if (not allow_charged) and is_charged_smiles(mutated_smi):
                continue

            if mutated_smi in seen:
                continue
            seen.add(mutated_smi)

            mol = smiles_to_mol(mutated_smi)
            if mol is None:
                continue

            sim = tanimoto_sim(seed_mol, mol)
            try:
                logp = Descriptors.MolLogP(mol)
                qed = Descriptors.qed(mol)
            except Exception:
                logp, qed = np.nan, np.nan

            records.append(
                {
                    "seed_index": seed_idx,
                    "seed_smiles": mol_to_smiles(seed_mol),
                    "mutated_smiles": mutated_smi,
                    "num_edits": n_edits,
                    "similarity_to_seed": sim,
                    "logP": logp,
                    "QED": qed,
                }
            )
            if len(records) >= max_total:
                break
        if len(records) >= max_total:
            break

    df = pd.DataFrame(records)
    return df.sort_values(["seed_index", "similarity_to_seed"], ascending=[True, False])


# --------------------------
# 2. Greedy chemical path between two molecules
# --------------------------

def greedy_path_between_two(
    start_smiles: str,
    target_smiles: str,
    alphabet: Set[str],
    max_steps: int = 20,
    neighbors_per_step: int = 128,
    random_seed: int = 0,
    allow_charged: bool = True,
) -> pd.DataFrame:
    """Create a greedy SELFIES mutation path from start to target."""
    random.seed(random_seed)
    np.random.seed(random_seed)

    alphabet_list = list(alphabet)

    start_mol = smiles_to_mol(start_smiles)
    target_mol = smiles_to_mol(target_smiles)
    if start_mol is None or target_mol is None:
        return pd.DataFrame()

    target_fp = mol_fingerprint(target_mol)
    current_sf = smiles_to_selfies(start_smiles)
    if current_sf is None:
        return pd.DataFrame()

    path_records = []

    def sim_to_target(mol: Chem.Mol) -> float:
        if mol is None:
            return 0.0
        fp = mol_fingerprint(mol)
        return float(DataStructs.TanimotoSimilarity(fp, target_fp))

    current_smi = selfies_to_smiles(current_sf)
    current_mol = smiles_to_mol(current_smi)
    current_sim = sim_to_target(current_mol)

    path_records.append(
        {
            "step": 0,
            "smiles": current_smi,
            "similarity_to_target": current_sim,
        }
    )

    for step in range(1, max_steps + 1):
        candidates: List[Tuple[float, str, str]] = []
        seen_smis: Set[str] = set()

        for _ in range(neighbors_per_step):
            mutated_sf = random_selfies_mutation(current_sf, alphabet_list, num_edits=1)
            mutated_smi = selfies_to_smiles(mutated_sf)
            if mutated_smi is None:
                continue

            if (not allow_charged) and is_charged_smiles(mutated_smi):
                continue

            if mutated_smi in seen_smis:
                continue
            seen_smis.add(mutated_smi)
            mol = smiles_to_mol(mutated_smi)
            if mol is None:
                continue
            sim = sim_to_target(mol)
            candidates.append((sim, mutated_sf, mutated_smi))

        if not candidates:
            break

        best_sim, best_sf, best_smi = max(candidates, key=lambda x: x[0])

        current_sf = best_sf
        current_smi = best_smi
        current_sim = best_sim

        path_records.append(
            {
                "step": step,
                "smiles": current_smi,
                "similarity_to_target": current_sim,
            }
        )

        if mol_to_smiles(target_mol) == current_smi:
            break

    df = pd.DataFrame(path_records)
    return df


# --------------------------
# 3. Median molecule search (any number of references)
# --------------------------

def median_objective(mol: Chem.Mol, references: List[Chem.Mol]) -> float:
    """Sum of Tanimoto similarities to reference molecules."""
    return float(sum(tanimoto_sim(mol, ref) for ref in references))


def median_molecule_search(
    smiles_list: List[str],
    alphabet: Set[str],
    n_starts: int = 3,
    steps_per_start: int = 20,
    neighbors_per_step: int = 128,
    use_hill_climb: bool = False,
    random_seed: int = 0,
    allow_charged: bool = True,
) -> pd.DataFrame:
    """Search for molecules simultaneously similar to multiple references."""
    random.seed(random_seed)
    np.random.seed(random_seed)

    alphabet_list = list(alphabet)

    ref_mols = [smiles_to_mol(s) for s in smiles_list]
    ref_mols = [m for m in ref_mols if m is not None]
    if len(ref_mols) < 2:
        return pd.DataFrame()

    all_candidates: Dict[str, float] = {}

    for m in ref_mols:
        smi = mol_to_smiles(m)
        all_candidates[smi] = median_objective(m, ref_mols)

    for start_idx in range(min(n_starts, len(ref_mols))):
        start_mol = ref_mols[start_idx]
        start_smi = mol_to_smiles(start_mol)
        current_sf = smiles_to_selfies(start_smi)
        if current_sf is None:
            continue

        current_smi = selfies_to_smiles(current_sf)
        current_mol = smiles_to_mol(current_smi)
        current_score = median_objective(current_mol, ref_mols)

        all_candidates[current_smi] = max(
            all_candidates.get(current_smi, -1e9),
            current_score,
        )

        for _ in range(steps_per_start):
            step_candidates: List[Tuple[float, str, str]] = []
            seen_smis_step: Set[str] = set()

            for _ in range(neighbors_per_step):
                mutated_sf = random_selfies_mutation(current_sf, alphabet_list, num_edits=1)
                mutated_smi = selfies_to_smiles(mutated_sf)
                if mutated_smi is None:
                    continue

                if (not allow_charged) and is_charged_smiles(mutated_smi):
                    continue

                if mutated_smi in seen_smis_step:
                    continue
                seen_smis_step.add(mutated_smi)

                mol = smiles_to_mol(mutated_smi)
                if mol is None:
                    continue

                score = median_objective(mol, ref_mols)
                step_candidates.append((score, mutated_sf, mutated_smi))

                prev = all_candidates.get(mutated_smi, -1e9)
                if score > prev:
                    all_candidates[mutated_smi] = score

            if not step_candidates:
                break

            best_score, best_sf, best_smi = max(step_candidates, key=lambda x: x[0])

            if use_hill_climb and best_score <= current_score:
                break

            current_sf = best_sf
            current_smi = best_smi
            current_mol = smiles_to_mol(current_smi)
            current_score = best_score

            all_candidates[current_smi] = max(
                all_candidates.get(current_smi, -1e9),
                current_score,
            )

    rows = []
    for smi, score in all_candidates.items():
        mol = smiles_to_mol(smi)
        if mol is None:
            continue
        sims = [tanimoto_sim(mol, r) for r in ref_mols]
        row = {
            "smiles": smi,
            "median_score_sum": score,
        }
        for i, s in enumerate(sims):
            row[f"sim_to_ref_{i+1}"] = s
        rows.append(row)

    df = pd.DataFrame(rows)
    if not df.empty:
        df = df.sort_values("median_score_sum", ascending=False)
    return df


__all__ = [
    "smiles_to_mol",
    "mol_to_smiles",
    "smiles_list_to_grid_image",
    "mol_fingerprint",
    "tanimoto_sim",
    "is_charged_smiles",
    "smiles_to_selfies",
    "selfies_to_smiles",
    "get_selfies_alphabet",
    "get_default_selfies_alphabet",
    "random_selfies_mutation",
    "generate_local_chemical_space",
    "greedy_path_between_two",
    "median_molecule_search",
]
