"""A1 接收器能力图谱：物理信号接收器在人体细胞类型中的表达分布.

数据源：Human Protein Atlas v25.1, rna_single_cell_type.tsv.zip
  https://www.proteinatlas.org/about/download
  列：Gene / Gene name / Cell type / nCPM（每百万归一化计数）

核心思想（与 vec-trust-kit 的可证伪原则一致）：
  - 不设单一"阳性/阴性"阈值就下结论，而是显式报告阈值敏感性；
  - 特异性是定量的（tau 指数），不是"某细胞会响应"的口头断言；
  - nCPM 只在基因内部跨细胞类型可比，不跨基因比绝对大小。

用法：
  from interface_atlas.receptors import build_matrix, load_matrix, audit, top_cells
  mat = build_matrix("data/rna_single_cell_type.tsv.zip")
  rows = load_matrix()
  report = audit()
"""
from __future__ import annotations

import csv
from pathlib import Path

import pandas as pd

DATA_DIR = Path(__file__).resolve().parents[1] / "data"
HPA_ZIP = DATA_DIR / "rna_single_cell_type.tsv.zip"
MATRIX_CSV = DATA_DIR / "receptor_expression.csv"

# gene -> (modality, receptor label, biological rationale)
RECEPTOR_PANEL: dict[str, tuple[str, str, str]] = {
    "COX4I1": ("photobiomodulation", "CCO subunit IV", "mitochondrial cytochrome c oxidase, red/NIR primary chromophore"),
    "CRY1": ("photobiomodulation / circadian", "Cryptochrome-1", "blue-light flavoprotein, circadian and redox signaling"),
    "CRY2": ("photobiomodulation / circadian", "Cryptochrome-2", "blue-light flavoprotein, circadian and redox signaling"),
    "CHRNA7": ("vagus / acetylcholine", "alpha7 nAChR", "cholinergic anti-inflammatory pathway on macrophages"),
    "PIEZO1": ("mechanical", "PIEZO1", "mechanosensitive Ca2+ channel, bone and cartilage loading"),
    "TRPV4": ("mechanical", "TRPV4", "mechanosensitive channel, osmotic and shear sensing"),
    "TRPV1": ("heat / capsaicin", "TRPV1", "thermal and chemical nocisensor"),
    "TRPM8": ("cold", "TRPM8", "cold-sensing channel"),
    "FFAR2": ("metabolic / microbiota", "GPR43 / FFAR2", "short-chain fatty acid receptor on intestinal L cells"),
    "FFAR3": ("metabolic / microbiota", "GPR41 / FFAR3", "short-chain fatty acid receptor on intestinal L cells"),
    "GLP1R": ("metabolic / incretin", "GLP-1 receptor", "target of endogenous GLP-1 and GLP-1 receptor agonists"),
    "GPBAR1": ("metabolic / bile acid", "TGR5 / GPBAR1", "bile acid receptor driving GLP-1 secretion"),
    "NR1H4": ("metabolic / bile acid", "FXR / NR1H4", "farnesoid X receptor, FGF15/19 and bile acid homeostasis"),
}

THRESHOLDS = (0.5, 1.0, 5.0, 10.0, 20.0)
DEFAULT_THRESHOLD = 1.0
LONG_FIELDS = (
    "gene", "modality", "receptor", "cell_type", "n_cpm", "pct_rank",
    "enrichment", "specificity_tau", "detected_1", "detected_5",
)


def _read_hpa(zip_path: Path | str) -> pd.DataFrame:
    path = Path(zip_path)
    if not path.exists():
        raise FileNotFoundError(f"HPA 数据文件不存在：{path}")
    df = pd.read_csv(
        path, sep="\t", compression="zip",
        usecols=["Gene name", "Cell type", "nCPM"],
        dtype={"Gene name": str, "Cell type": str, "nCPM": float},
    )
    return df[df["Gene name"].isin(RECEPTOR_PANEL)]


def build_matrix(zip_path: Path | str = HPA_ZIP) -> pd.DataFrame:
    """返回 index=细胞类型、columns=接收器基因、values=nCPM 的矩阵."""
    filtered = _read_hpa(zip_path)
    mat = filtered.pivot_table(
        index="Cell type", columns="Gene name", values="nCPM", aggfunc="mean", fill_value=0.0
    )
    return mat.reindex(columns=list(RECEPTOR_PANEL), fill_value=0.0).sort_index()


def to_long(mat: pd.DataFrame, threshold: float = DEFAULT_THRESHOLD) -> pd.DataFrame:
    """长表输出，每行 = 一个基因 × 一个细胞类型，含定量特异性与阈值标记."""
    out = []
    for gene in mat.columns:
        series = mat[gene].astype(float)
        mean = series.mean()
        pct = series.rank(pct=True) * 100
        enrichment = series / (mean + 1e-9)
        tau = float((series - mean).clip(lower=0).sum() / (series.sum() + 1e-9))
        modality, receptor, _ = RECEPTOR_PANEL[gene]
        for cell_type, value in series.items():
            out.append({
                "gene": gene,
                "modality": modality,
                "receptor": receptor,
                "cell_type": cell_type,
                "n_cpm": round(float(value), 4),
                "pct_rank": round(float(pct[cell_type]), 2),
                "enrichment": round(float(enrichment[cell_type]), 3),
                "specificity_tau": round(tau, 4),
                "detected_1": int(value >= 1.0),
                "detected_5": int(value >= 5.0),
            })
    _ = threshold  # 阈值敏感性由 audit() 统一报告，避免硬编码单一阈值
    return pd.DataFrame(out, columns=list(LONG_FIELDS)).sort_values(
        ["gene", "n_cpm"], ascending=[True, False]
    ).reset_index(drop=True)


def save_matrix(mat: pd.DataFrame, path: Path | str = MATRIX_CSV) -> Path:
    p = Path(path)
    p.parent.mkdir(parents=True, exist_ok=True)
    to_long(mat).to_csv(p, index=False, encoding="utf-8-sig")
    return p


def load_matrix(path: Path | str = MATRIX_CSV) -> pd.DataFrame:
    """读取已生成的长表；未生成则明确报错，不静默伪造数据."""
    p = Path(path)
    if not p.exists():
        raise FileNotFoundError(
            f"未找到 {p}。请先运行 python -m interface_atlas.cli receptors --rebuild"
        )
    return pd.read_csv(p, encoding="utf-8-sig")


def top_cells(df: pd.DataFrame, *, gene: str | None = None, modality: str | None = None,
              n: int = 10) -> pd.DataFrame:
    sub = df
    if gene:
        sub = sub[sub["gene"] == gene]
    if modality:
        sub = sub[sub["modality"].str.contains(modality, case=False, na=False)]
    return sub.sort_values("n_cpm", ascending=False).head(n)


def threshold_sensitivity(mat: pd.DataFrame) -> dict[str, dict[float, int]]:
    """每个基因在不同 nCPM 阈值下"阳性"细胞类型数量——显式披露阈值依赖."""
    result: dict[str, dict[float, int]] = {}
    for gene in mat.columns:
        s = mat[gene]
        result[gene] = {t: int((s >= t).sum()) for t in THRESHOLDS}
    return result


def audit(df: pd.DataFrame, mat: pd.DataFrame | None = None) -> dict:
    """数据完整性 + 阈值敏感性 + 选择性预警。不把任何阈值当真理."""
    issues: list[str] = []
    expected = set(RECEPTOR_PANEL)
    found = set(df["gene"].unique())
    if found != expected:
        issues.append(f"基因覆盖不完整：缺 {sorted(expected - found)}，多 {sorted(found - expected)}")
    nulls = int(df.isna().sum().sum())
    if nulls:
        issues.append(f"存在 {nulls} 个空值")
    cell_types = int(df["cell_type"].nunique())
    if cell_types < 100:
        issues.append(f"细胞类型数异常偏少：{cell_types}")

    low_selectivity: list[str] = []
    for gene, group in df.groupby("gene"):
        if group[group["n_cpm"] >= DEFAULT_THRESHOLD].shape[0] / max(cell_types, 1) > 0.8:
            low_selectivity.append(gene)  # 广布接收器，不应宣称细胞特异

    return {
        "status": "PASS" if not issues else "FAIL",
        "issues": issues,
        "total_rows": len(df),
        "genes": len(found),
        "cell_types": cell_types,
        "n_cpm_min": float(df["n_cpm"].min()),
        "n_cpm_max": float(df["n_cpm"].max()),
        "low_selectivity_genes_at_1ncpm": low_selectivity,
        "threshold_sensitivity": threshold_sensitivity(mat) if mat is not None else None,
        "caveat": "nCPM 仅在同一基因内跨细胞类型可比；阈值 1 nCPM 是约定的探索性截断，不是生物学真理。",
    }


def write_audit(report: dict, path: Path | str) -> Path:
    p = Path(path)
    p.parent.mkdir(parents=True, exist_ok=True)
    lines = [
        "# A1 接收器能力图谱审计",
        "",
        f"- 状态：**{report['status']}**",
        f"- 行数：{report['total_rows']}；基因：{report['genes']}；细胞类型：{report['cell_types']}",
        f"- nCPM 范围：{report['n_cpm_min']} – {report['n_cpm_max']}",
        "",
        "## 问题",
        "",
    ]
    lines += [f"- {x}" for x in report["issues"]] or ["- 无"]
    lines += ["", "## 低选择性预警（≥1 nCPM 的细胞类型 > 80%）", ""]
    lines += [f"- {g}" for g in report["low_selectivity_genes_at_1ncpm"]] or ["- 无"]
    lines += ["", "## 阈值敏感性：阳性细胞类型数", "", "| 基因 | " + " | ".join(str(t) for t in THRESHOLDS) + " |", "|---|---|" + "---|" * len(THRESHOLDS)]
    if report["threshold_sensitivity"]:
        for gene, counts in report["threshold_sensitivity"].items():
            lines.append(f"| {gene} | " + " | ".join(str(counts[t]) for t in THRESHOLDS) + " |")
    lines += ["", "> " + report["caveat"], ""]
    p.write_text("\n".join(lines), encoding="utf-8")
    return p


def write_heatmap(df: pd.DataFrame, path: Path | str, top_n: int = 40) -> Path:
    """中文热图：13 接收器 × 综合排名前 top_n 的细胞类型，单元格为基因内百分位."""
    import warnings
    warnings.filterwarnings("ignore", category=UserWarning)
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib import font_manager

    available = {f.name for f in font_manager.fontManager.ttflist}
    for name in ("Microsoft YaHei", "SimHei", "Noto Sans CJK SC"):
        if name in available:
            plt.rcParams["font.sans-serif"] = [name, "DejaVu Sans"]
            break
    plt.rcParams["axes.unicode_minus"] = False

    wide = df.pivot(index="cell_type", columns="gene", values="pct_rank")
    order = wide.mean(axis=1).sort_values(ascending=False).head(top_n).index
    data = wide.loc[order][list(RECEPTOR_PANEL)].T

    fig, ax = plt.subplots(figsize=(18, 5.6))
    im = ax.imshow(data.values, aspect="auto", cmap="YlGnBu", vmin=0, vmax=100)
    ax.set_yticks(range(data.shape[0]))
    ax.set_yticklabels([RECEPTOR_PANEL[g][1] for g in data.index], fontsize=9)
    ax.set_xticks(range(data.shape[1]))
    ax.set_xticklabels(order, rotation=60, ha="right", fontsize=7)
    ax.set_title("物理信号接收器在人体细胞类型中的表达百分位（HPA v25.1, nCPM）", fontsize=12)
    fig.colorbar(im, ax=ax, label="基因内百分位", shrink=0.8)
    fig.tight_layout()
    p = Path(path)
    p.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(p, dpi=150)
    plt.close(fig)
    return p
