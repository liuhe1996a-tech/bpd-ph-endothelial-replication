"""Create a donor-level technical and biological confounding audit."""

from __future__ import annotations

import os

import json
from pathlib import Path

import pandas as pd


PROJECT_ROOT = Path(__file__).resolve().parents[2]
SOURCE_ROOT = Path(os.environ.get("GSE275938_RAW_ROOT", PROJECT_ROOT / "02_raw/GSE275938"))
PROCESSED_DIR = PROJECT_ROOT / "04_processed" / "GSE275938"
RESULT_DIR = PROJECT_ROOT / "07_results" / "human_endothelial_subtypes"
SAMPLES = ["BPD 1", "BPD 2", "BPD+PH 1", "BPD+PH 2"]
GSM = {
    "BPD 1": "GSM8488396",
    "BPD 2": "GSM8488397",
    "BPD+PH 1": "GSM8488398",
    "BPD+PH 2": "GSM8488399",
}


def main() -> int:
    RESULT_DIR.mkdir(parents=True, exist_ok=True)
    all_cells = pd.read_csv(
        SOURCE_ROOT / "GSE275938_cell_metadata.csv.gz",
        low_memory=False,
    )
    disease_endothelial = pd.read_csv(
        PROCESSED_DIR
        / "GSE275938_disease_endothelial_cell_metadata.tsv.gz",
        sep="\t",
    )
    rows = []
    for sample in SAMPLES:
        condition = "BPD+PH" if sample.startswith("BPD+PH") else "BPD"
        all_subset = all_cells.loc[all_cells["dataset"].eq(sample)]
        endo = disease_endothelial.loc[
            disease_endothelial["sample"].eq(sample)
        ]
        rows.append(
            {
                "sample": sample,
                "GEO_accession": GSM[sample],
                "diagnosis_group": condition,
                "tissue_source": "lung tissue from an infant who died",
                "annotated_cells_all_lineages": len(all_subset),
                "recovered_endothelial_cells": len(endo),
                "endothelial_fraction_of_annotated_cells": (
                    len(endo) / len(all_subset)
                ),
                "endothelial_total_UMI": int(endo["cell_total_umi"].sum()),
                "endothelial_median_UMI_per_cell": float(
                    endo["metadata_nCount_RNA"].median()
                ),
                "endothelial_median_features_per_cell": float(
                    endo["metadata_nFeature_RNA"].median()
                ),
                "endothelial_median_percent_mt": float(
                    endo["metadata_percent_mt"].median()
                ),
                "all_cell_median_UMI": float(
                    all_subset["nCount_RNA"].median()
                ),
                "all_cell_median_features": float(
                    all_subset["nFeature_RNA"].median()
                ),
                "all_cell_median_percent_mt": float(
                    all_subset["percent.mt"].median()
                ),
                "platform": "GPL24676",
                "instrument": "Illumina NovaSeq 6000",
                "library_strategy": "10x Chromium 3-prime or 5-prime RNA-seq",
                "sample_specific_3p_or_5p_kit_deposited": False,
                "alignment": "Cell Ranger Count 7.1.0 to GRCh38",
                "ambient_RNA_removal": "CellBender 0.2.1",
                "deposited_QC": (
                    "quality<30 reads removed; cells filtered at >10% mt, "
                    "<500 genes, or >7000 genes"
                ),
            }
        )
    audit = pd.DataFrame(rows)
    audit.to_csv(
        RESULT_DIR / "GSE275938_donor_technical_audit.tsv",
        sep="\t",
        index=False,
    )

    condition_summary = (
        audit.groupby("diagnosis_group", as_index=False)
        .agg(
            donors=("sample", "nunique"),
            median_endothelial_cells=(
                "recovered_endothelial_cells",
                "median",
            ),
            median_endothelial_UMI_per_cell=(
                "endothelial_median_UMI_per_cell",
                "median",
            ),
            median_endothelial_features_per_cell=(
                "endothelial_median_features_per_cell",
                "median",
            ),
            median_endothelial_percent_mt=(
                "endothelial_median_percent_mt",
                "median",
            ),
        )
    )
    condition_summary.to_csv(
        RESULT_DIR / "GSE275938_donor_technical_audit_by_condition.tsv",
        sep="\t",
        index=False,
    )
    bpd = condition_summary.set_index("diagnosis_group").loc["BPD"]
    ph = condition_summary.set_index("diagnosis_group").loc["BPD+PH"]
    summary = {
        "biological_unit": "donor",
        "disease_groups": {"BPD": 2, "BPD+PH": 2},
        "same_deposited_platform_and_instrument": True,
        "sample_specific_3p_or_5p_library_assignment_available": False,
        "median_endothelial_cell_recovery_ratio_BPD_PH_over_BPD": float(
            ph["median_endothelial_cells"]
            / bpd["median_endothelial_cells"]
        ),
        "median_endothelial_UMI_ratio_BPD_over_BPD_PH": float(
            bpd["median_endothelial_UMI_per_cell"]
            / ph["median_endothelial_UMI_per_cell"]
        ),
        "median_endothelial_feature_ratio_BPD_over_BPD_PH": float(
            bpd["median_endothelial_features_per_cell"]
            / ph["median_endothelial_features_per_cell"]
        ),
        "interpretation": (
            "Diagnosis is strongly confounded with endothelial recovery and "
            "per-cell UMI/feature depth. The deposited record does not "
            "identify each sample's 3-prime versus 5-prime kit, so protocol "
            "confounding cannot be excluded. Human results are exploratory "
            "directional support, not an independent disease association."
        ),
    }
    (
        RESULT_DIR / "GSE275938_donor_technical_audit_summary.json"
    ).write_text(
        json.dumps(summary, ensure_ascii=False, indent=2),
        encoding="utf-8",
    )
    print(audit.to_string(index=False))
    print(json.dumps(summary, ensure_ascii=False, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
