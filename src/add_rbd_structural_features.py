from pathlib import Path
import urllib.request
import math

import numpy as np
import pandas as pd

RAW_DIR = Path("data/raw")
PROCESSED_DIR = Path("data/processed")
OUT_DIR = Path("artifacts/tables")

RAW_DIR.mkdir(parents=True, exist_ok=True)
PROCESSED_DIR.mkdir(parents=True, exist_ok=True)
OUT_DIR.mkdir(parents=True, exist_ok=True)

PDB_ID = "6M0J"
PDB_URL = f"https://files.rcsb.org/download/{PDB_ID}.pdb"
PDB_PATH = RAW_DIR / f"{PDB_ID}.pdb"

IN_PATH = PROCESSED_DIR / "rbd_dms_substitutions_normalized_v0.csv"
OUT_PATH = PROCESSED_DIR / "rbd_dms_substitutions_structural_v0.csv"
SITE_FEATURE_PATH = OUT_DIR / "rbd_site_structural_features_6m0j.csv"
SUMMARY_PATH = OUT_DIR / "rbd_structural_feature_summary.csv"

# Rough receptor-binding motif range used only as a simple annotation feature.
RBM_START = 438
RBM_END = 506

# Distance threshold for defining RBD residues near the ACE2 interface in the 6M0J structure.
CONTACT_DISTANCE_ANGSTROM = 5.0


def download_pdb():
    if PDB_PATH.exists() and PDB_PATH.stat().st_size > 1000:
        print(f"Already exists: {PDB_PATH}")
        return

    print(f"Downloading PDB structure: {PDB_URL}")
    urllib.request.urlretrieve(PDB_URL, PDB_PATH)


def parse_pdb_atoms(path):
    atoms = []

    with open(path, "r", encoding="utf-8", errors="ignore") as f:
        for line in f:
            if not line.startswith("ATOM"):
                continue

            atom_name = line[12:16].strip()
            resname = line[17:20].strip()
            chain = line[21].strip()
            resseq_raw = line[22:26].strip()
            icode = line[26].strip()

            if not resseq_raw:
                continue

            try:
                resseq = int(resseq_raw)
                x = float(line[30:38])
                y = float(line[38:46])
                z = float(line[46:54])
            except ValueError:
                continue

            atoms.append(
                {
                    "atom_name": atom_name,
                    "resname": resname,
                    "chain": chain,
                    "resseq": resseq,
                    "icode": icode,
                    "x": x,
                    "y": y,
                    "z": z,
                }
            )

    atoms_df = pd.DataFrame(atoms)

    if atoms_df.empty:
        raise RuntimeError("No ATOM records were parsed from the PDB file.")

    return atoms_df


def infer_rbd_chain(atoms_df):
    # In 6M0J, the RBD chain is expected to contain many Spike RBD residue numbers.
    chain_scores = []

    for chain, group in atoms_df.groupby("chain"):
        ca = group[group["atom_name"] == "CA"]
        n_rbd_numbered = ca[ca["resseq"].between(319, 541)]["resseq"].nunique()
        n_total = ca["resseq"].nunique()
        chain_scores.append(
            {
                "chain": chain,
                "n_rbd_numbered_residues": n_rbd_numbered,
                "n_total_residues": n_total,
            }
        )

    scores = pd.DataFrame(chain_scores).sort_values(
        ["n_rbd_numbered_residues", "n_total_residues"],
        ascending=False,
    )

    if scores.empty or scores.iloc[0]["n_rbd_numbered_residues"] == 0:
        raise RuntimeError(
            "Could not infer the RBD chain from PDB residue numbering. "
            "Check the PDB file and chain IDs manually."
        )

    rbd_chain = str(scores.iloc[0]["chain"])
    print("\nInferred chain summary:")
    print(scores.to_string(index=False))
    print(f"\nSelected RBD chain: {rbd_chain}")

    return rbd_chain, scores


def compute_site_features(atoms_df, rbd_chain):
    rbd_atoms = atoms_df[
        (atoms_df["chain"] == rbd_chain)
        & (atoms_df["resseq"].between(319, 541))
    ].copy()

    ace2_atoms = atoms_df[atoms_df["chain"] != rbd_chain].copy()

    if rbd_atoms.empty:
        raise RuntimeError("No RBD atoms found after selecting the inferred RBD chain.")
    if ace2_atoms.empty:
        raise RuntimeError("No non-RBD atoms found for interface distance calculation.")

    ace2_xyz = ace2_atoms[["x", "y", "z"]].to_numpy(dtype=float)
    ace2_residue_ids = (
        ace2_atoms["chain"].astype(str)
        + ":"
        + ace2_atoms["resseq"].astype(str)
        + ace2_atoms["icode"].fillna("").astype(str)
    ).to_numpy()

    rows = []

    for site, group in rbd_atoms.groupby("resseq"):
        rbd_xyz = group[["x", "y", "z"]].to_numpy(dtype=float)

        min_dist = math.inf
        contact_residues = set()

        # The structure is small enough for a simple loop over atoms.
        for atom_xyz in rbd_xyz:
            diff = ace2_xyz - atom_xyz
            dists = np.sqrt((diff * diff).sum(axis=1))
            local_min = float(dists.min())
            if local_min < min_dist:
                min_dist = local_min

            contact_mask = dists <= CONTACT_DISTANCE_ANGSTROM
            for residue_id in ace2_residue_ids[contact_mask]:
                contact_residues.add(str(residue_id))

        rows.append(
            {
                "site": int(site),
                "structure_source": PDB_ID,
                "rbd_chain": rbd_chain,
                "is_rbm_rough": int(RBM_START <= int(site) <= RBM_END),
                "min_dist_to_ace2_angstrom": min_dist,
                "is_ace2_interface_5A": int(min_dist <= CONTACT_DISTANCE_ANGSTROM),
                "n_ace2_contact_residues_5A": len(contact_residues),
                "structure_site_available": 1,
            }
        )

    site_features = pd.DataFrame(rows).sort_values("site")
    return site_features


def main():
    print("Reading normalized substitution dataset...")
    df = pd.read_csv(IN_PATH)
    print(f"Input shape: {df.shape}")

    download_pdb()
    atoms_df = parse_pdb_atoms(PDB_PATH)
    rbd_chain, chain_summary = infer_rbd_chain(atoms_df)
    site_features = compute_site_features(atoms_df, rbd_chain)

    site_features.to_csv(SITE_FEATURE_PATH, index=False, encoding="utf-8-sig")

    merged = df.merge(site_features, on="site", how="left")

    merged["structure_site_available"] = merged["structure_site_available"].fillna(0).astype(int)
    merged["is_rbm_rough"] = merged["is_rbm_rough"].fillna(
        ((merged["site"] >= RBM_START) & (merged["site"] <= RBM_END)).astype(int)
    ).astype(int)
    merged["is_ace2_interface_5A"] = merged["is_ace2_interface_5A"].fillna(0).astype(int)
    merged["n_ace2_contact_residues_5A"] = merged["n_ace2_contact_residues_5A"].fillna(0).astype(int)

    # Missing distance means that the residue was not available in the experimental structure.
    # Use a large distance so models treat it as far from ACE2, while keeping availability as a separate feature.
    merged["min_dist_to_ace2_angstrom"] = merged["min_dist_to_ace2_angstrom"].fillna(999.0)

    merged.to_csv(OUT_PATH, index=False, encoding="utf-8-sig")

    summary = pd.DataFrame(
        [
            {
                "input_rows": len(df),
                "output_rows": len(merged),
                "n_sites_in_dataset": df["site"].nunique(),
                "n_sites_with_structure": site_features["site"].nunique(),
                "n_sites_marked_interface_5A": int(site_features["is_ace2_interface_5A"].sum()),
                "n_rows_with_structure": int(merged["structure_site_available"].sum()),
                "n_rows_marked_interface_5A": int(merged["is_ace2_interface_5A"].sum()),
                "rbd_chain": rbd_chain,
                "pdb_id": PDB_ID,
                "contact_distance_angstrom": CONTACT_DISTANCE_ANGSTROM,
            }
        ]
    )

    summary.to_csv(SUMMARY_PATH, index=False, encoding="utf-8-sig")

    print("\nSaved site-level structural features:")
    print(SITE_FEATURE_PATH)

    print("\nSaved merged structural dataset:")
    print(OUT_PATH)

    print("\nSaved structural feature summary:")
    print(SUMMARY_PATH)

    print("\nSummary:")
    print(summary.to_string(index=False))

    print("\nTop interface-like RBD sites by distance:")
    print(
        site_features.sort_values("min_dist_to_ace2_angstrom")
        .head(30)
        .to_string(index=False)
    )


if __name__ == "__main__":
    main()
