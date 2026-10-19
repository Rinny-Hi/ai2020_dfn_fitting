"""Create a one-page summary figure for anode-OCP selection."""

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.font_manager as font_manager
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

import ocp_dynamic_combination_comparison as ocp_compare


ROOT = Path(__file__).resolve().parent
RESULT = ROOT / "results" / "261019_lowrate_anode_full_pipeline_comparison"
OUTPUT = RESULT / "anode_ocp_selection_onepage.png"

font_path = Path(r"C:\Windows\Fonts\malgun.ttf")
if font_path.exists():
    font_manager.fontManager.addfont(str(font_path))
    plt.rcParams["font.family"] = font_manager.FontProperties(fname=font_path).get_name()
plt.rcParams["axes.unicode_minus"] = False

AI_COLOR = "#D55E00"
LOW_COLOR = "#0072B2"
NAVY = "#163A5F"
GRID = "#D9DEE5"
MUTED = "#5E6875"
GREEN = "#1B7F5A"


def load_anode_ocp_curves() -> pd.DataFrame:
    """Load once, then cache the two anode OCP curves used in the comparison."""

    cache = RESULT / "anode_ocp_curve_comparison.csv"
    if cache.exists():
        return pd.read_csv(cache)
    anodes, _ = ocp_compare.load_ocp_sources({})
    rows = []
    for source, key in (("Ai2020", "nominal"), ("Low-rate", "old")):
        for x_value, voltage in zip(
            anodes[key]["grid"], anodes[key]["equilibrium"]
        ):
            rows.append(
                {"Source": source, "Stoichiometry_x": x_value, "OCP_V": voltage}
            )
    curves = pd.DataFrame(rows)
    curves.to_csv(cache, index=False, encoding="utf-8-sig")
    return curves


def branch_value(frame, scenario, rate, direction, column):
    row = frame[
        (frame.Scenario == scenario)
        & np.isclose(frame.C_rate.astype(float), rate)
        & (frame.Direction == direction)
    ].iloc[0]
    return float(row[column])


def main() -> None:
    branch = pd.read_csv(RESULT / "pipeline_branch_metrics.csv")
    pipeline = pd.read_csv(RESULT / "pipeline_summary.csv")
    fitted = pd.read_csv(RESULT / "s6_s7a_summary.csv")
    source_windows = pd.read_csv(
        ROOT
        / "results"
        / "260921_c50_ocp_source_combinations"
        / "ocp_source_combination_summary.csv"
    )
    curves = load_anode_ocp_curves()

    figure = plt.figure(figsize=(18, 10.3), facecolor="white")
    grid = figure.add_gridspec(
        2,
        3,
        width_ratios=(1.0, 1.25, 1.0),
        height_ratios=(1.0, 1.0),
        left=0.045,
        right=0.98,
        top=0.87,
        bottom=0.065,
        hspace=0.32,
        wspace=0.28,
    )
    ax_ocp = figure.add_subplot(grid[0, 0])
    ax_branch = figure.add_subplot(grid[0, 1:])
    ax_electrolyte = figure.add_subplot(grid[1, 0])
    ax_table = figure.add_subplot(grid[1, 1])
    ax_decision = figure.add_subplot(grid[1, 2])

    figure.suptitle(
        "음극 OCP 후보 비교 및 최종 선정",
        x=0.045,
        y=0.965,
        ha="left",
        fontsize=24,
        fontweight="bold",
        color=NAVY,
    )
    figure.text(
        0.045,
        0.91,
        "양극 GITT 고정  |  S5 → S6: ORegan2022 전해질 + $b_n=2.5$  |  S7a: $k_n+b_n$ 전압 fitting",
        ha="left",
        fontsize=12.5,
        color=MUTED,
    )

    # Panel 1: anode OCP shape and operational stoichiometry windows.
    window_rows = {
        "Ai2020": source_windows[
            (source_windows.anode_source == "nominal")
            & (source_windows.cathode_source == "old")
        ].iloc[0],
        "Low-rate": source_windows[
            (source_windows.anode_source == "old")
            & (source_windows.cathode_source == "old")
        ].iloc[0],
    }
    for source, color, marker, label in (
        ("Ai2020", AI_COLOR, "o", "Ai2020 graphite"),
        ("Low-rate", LOW_COLOR, "s", "실험 저율 OCP"),
    ):
        part = curves[curves.Source == source].sort_values("Stoichiometry_x")
        window = window_rows[source]
        x0, x100 = float(window.x0), float(window.x100)
        operational = part[
            (part.Stoichiometry_x >= x0) & (part.Stoichiometry_x <= x100)
        ]
        ax_ocp.plot(
            operational.Stoichiometry_x,
            operational.OCP_V,
            color=color,
            lw=2.4,
            label=label,
        )
        endpoint_x = np.array([x0, x100])
        endpoint_v = np.interp(
            endpoint_x, part.Stoichiometry_x, part.OCP_V
        )
        ax_ocp.scatter(
            endpoint_x,
            endpoint_v,
            marker=marker,
            s=55,
            facecolor="white",
            edgecolor=color,
            lw=1.8,
            zorder=4,
        )
    ai_window = window_rows["Ai2020"]
    low_window = window_rows["Low-rate"]
    ax_ocp.text(
        0.98,
        0.95,
        f"Ai2020: x={float(ai_window.x0):.3f}→{float(ai_window.x100):.3f},  Qn={float(ai_window.Qn_Ah):.2f} Ah\n"
        f"저율: x={float(low_window.x0):.3f}→{float(low_window.x100):.3f},  Qn={float(low_window.Qn_Ah):.2f} Ah",
        transform=ax_ocp.transAxes,
        ha="right",
        va="top",
        fontsize=9.2,
        color=MUTED,
        bbox=dict(facecolor="white", edgecolor="none", alpha=0.86, pad=2.5),
    )
    ax_ocp.set_xlim(0, 1.0)
    ax_ocp.set_ylim(0.04, 0.90)
    ax_ocp.set_xlabel("음극 stoichiometry, x")
    ax_ocp.set_ylabel("OCP vs Li/Li+ [V]")
    ax_ocp.set_title("① 음극 OCP 형상·window 차이", loc="left", fontweight="bold", color=NAVY)
    ax_ocp.grid(color=GRID, lw=0.8)
    ax_ocp.spines[["top", "right"]].set_visible(False)
    ax_ocp.legend(
        frameon=False,
        loc="upper right",
        bbox_to_anchor=(0.99, 0.79),
        fontsize=9.5,
    )
    ax_ocp.text(
        0.02,
        0.03,
        "저율 x축은 Ai2020 inventory 기반 조건부 환산",
        transform=ax_ocp.transAxes,
        fontsize=8.8,
        color=MUTED,
    )

    # Panel 2: S6 branch-level comparison.
    conditions = [
        (0.5, "Charge", "0.5C\n충전"),
        (0.5, "Discharge", "0.5C\n방전"),
        (1.0, "Charge", "1C\n충전"),
        (1.0, "Discharge", "1C\n방전"),
        (2.0, "Charge", "2C\n충전"),
        (2.0, "Discharge", "2C\n방전"),
    ]
    x = np.arange(len(conditions))
    ai_values = [
        branch_value(branch, "S6 Ai2020", rate, direction, "Full_time_RMSE_mV")
        for rate, direction, _ in conditions
    ]
    low_values = [
        branch_value(branch, "S6 Low-rate", rate, direction, "Full_time_RMSE_mV")
        for rate, direction, _ in conditions
    ]
    ax_branch.plot(
        x,
        ai_values,
        "o-",
        color=AI_COLOR,
        lw=2.5,
        ms=7,
        label="Ai2020 음극",
        zorder=3,
    )
    ax_branch.plot(
        x,
        low_values,
        "s--",
        color=LOW_COLOR,
        lw=2.2,
        ms=6.5,
        label="실험 저율 음극",
        zorder=3,
    )
    for index, (ai_value, low_value) in enumerate(zip(ai_values, low_values)):
        ax_branch.text(index, ai_value - 1.4, f"{ai_value:.1f}", ha="center", va="top", fontsize=9, color=AI_COLOR)
        ax_branch.text(index, low_value + 1.25, f"{low_value:.1f}", ha="center", va="bottom", fontsize=9, color=LOW_COLOR)
    ax_branch.set_xticks(x, [item[2] for item in conditions])
    ax_branch.set_ylim(14, 59)
    ax_branch.set_ylabel("전체 시간축 전압 RMSE [mV]")
    ax_branch.set_title("② 전해질 변경 후 S6: C-rate별 충·방전", loc="left", fontweight="bold", color=NAVY)
    ax_branch.grid(axis="y", color=GRID, lw=0.8)
    ax_branch.spines[["top", "right"]].set_visible(False)
    ax_branch.legend(frameon=False, ncol=2, loc="lower left")

    # Panel 3: electrolyte-set effect.
    labels = ["Ai2020 음극", "실험 저율 음극"]
    stage_pairs = [("Ai2020", AI_COLOR), ("Low-rate", LOW_COLOR)]
    y = np.array([1.0, 0.0])
    for y_value, (name, color) in zip(y, stage_pairs):
        s5 = float(pipeline[(pipeline.Stage == "S5") & (pipeline.Anode_OCP == name)].Overall_full_time_RMSE_mV.iloc[0])
        s6 = float(pipeline[(pipeline.Stage == "S6") & (pipeline.Anode_OCP == name)].Overall_full_time_RMSE_mV.iloc[0])
        cap5 = float(pipeline[(pipeline.Stage == "S5") & (pipeline.Anode_OCP == name)].Capacity_RMSE_pct.iloc[0])
        cap6 = float(pipeline[(pipeline.Stage == "S6") & (pipeline.Anode_OCP == name)].Capacity_RMSE_pct.iloc[0])
        ax_electrolyte.plot([s6, s5], [y_value, y_value], color=color, lw=3, alpha=0.55)
        ax_electrolyte.scatter(s5, y_value, s=75, facecolor="white", edgecolor=color, lw=2, zorder=3)
        ax_electrolyte.scatter(s6, y_value, s=90, facecolor=color, edgecolor=color, lw=1, zorder=3)
        ax_electrolyte.text(s5 + 0.45, y_value + 0.08, f"S5 {s5:.2f}", color=color, fontsize=10, va="bottom")
        ax_electrolyte.text(s6 - 0.45, y_value + 0.08, f"S6 {s6:.2f}", color=color, fontsize=10, va="bottom", ha="right")
        ax_electrolyte.text(
            55.3,
            y_value - 0.17,
            f"용량 RMSE {cap5:.2f}→{cap6:.2f}%",
            color=MUTED,
            fontsize=9.5,
            ha="right",
            va="top",
        )
    ax_electrolyte.set_yticks(y, labels)
    ax_electrolyte.set_xlim(34, 56)
    ax_electrolyte.set_ylim(-0.55, 1.45)
    ax_electrolyte.set_xlabel("전체 전압 RMSE [mV]  (낮을수록 우수)")
    ax_electrolyte.set_title("③ 전해질 변경 효과: S5 → S6", loc="left", fontweight="bold", color=NAVY)
    ax_electrolyte.grid(axis="x", color=GRID, lw=0.8)
    ax_electrolyte.spines[["top", "right", "left"]].set_visible(False)
    ax_electrolyte.tick_params(axis="y", length=0)
    ax_electrolyte.text(34.2, 1.32, "● S6 ORegan2022     ○ S5 Ai2020 보정 전해질", fontsize=9.5, color=MUTED)

    # Panel 4: final S7a metrics as a decision table.
    ai_fit = fitted[fitted.Scenario == "Ai2020 S7a"].iloc[0]
    low_fit = fitted[fitted.Scenario == "Low-rate S7a"].iloc[0]
    table_rows = [
        ("전체 전압 RMSE", ai_fit.Overall_full_time_RMSE_mV, low_fit.Overall_full_time_RMSE_mV, "mV", "Ai2020"),
        ("중앙 10–70% RMSE", ai_fit.Overall_center_time_RMSE_mV, low_fit.Overall_center_time_RMSE_mV, "mV", "저율"),
        ("용량 RMSE", ai_fit.Capacity_RMSE_pct, low_fit.Capacity_RMSE_pct, "%", "저율"),
        ("최대 용량 오차", ai_fit.Max_abs_capacity_error_mAh, low_fit.Max_abs_capacity_error_mAh, "mAh", "Ai2020"),
        ("최소 충전 음극전위", ai_fit.Minimum_charge_anode_potential_mV, low_fit.Minimum_charge_anode_potential_mV, "mV", "Ai2020"),
        ("qOCV RMSE", 6.2427, 8.2196, "mV", "Ai2020"),
    ]
    ax_table.axis("off")
    ax_table.set_title("④ 동일 S7a fitting 후 핵심 지표", loc="left", fontweight="bold", color=NAVY, pad=10)
    cell_text = []
    for metric, ai_value, low_value, unit, better in table_rows:
        suffix = "" if metric == "최소 충전 음극전위" else ""
        cell_text.append([metric, f"{ai_value:.2f} {unit}{suffix}", f"{low_value:.2f} {unit}{suffix}", better])
    table = ax_table.table(
        cellText=cell_text,
        colLabels=["지표", "Ai2020 음극", "실험 저율 음극", "우세"],
        cellLoc="center",
        colLoc="center",
        bbox=[0.0, 0.12, 1.0, 0.80],
        colWidths=[0.34, 0.22, 0.24, 0.20],
    )
    table.auto_set_font_size(False)
    table.set_fontsize(9.2)
    for (row, column), cell in table.get_celld().items():
        cell.set_edgecolor("#E2E6EB")
        cell.set_linewidth(0.8)
        if row == 0:
            cell.set_facecolor(NAVY)
            cell.get_text().set_color("white")
            cell.get_text().set_fontweight("bold")
        else:
            cell.set_facecolor("#F6F8FA" if row % 2 == 0 else "white")
            if column == 3:
                cell.get_text().set_color(GREEN if cell.get_text().get_text() == "Ai2020" else LOW_COLOR)
                cell.get_text().set_fontweight("bold")
    ax_table.text(
        0.0,
        0.0,
        "$k_n$ 변화: Ai2020 -19.7% / 저율 -38.3%     |     $b_n$ 변화: Ai2020 -0.2% / 저율 +4.3%",
        transform=ax_table.transAxes,
        fontsize=9.2,
        color=MUTED,
    )

    # Panel 5: recommendation.
    ax_decision.axis("off")
    ax_decision.set_title("⑤ 최종 선정", loc="left", fontweight="bold", color=NAVY, pad=10)
    ax_decision.text(
        0.0,
        0.82,
        "Ai2020 음극 OCP 선정",
        fontsize=18,
        fontweight="bold",
        color=AI_COLOR,
        transform=ax_decision.transAxes,
    )
    ax_decision.text(
        0.0,
        0.68,
        "실험 GITT 양극 + Ai2020 음극\n+ ORegan2022 전해질 + $b_n=2.5$",
        fontsize=12.5,
        fontweight="bold",
        color=NAVY,
        transform=ax_decision.transAxes,
        linespacing=1.45,
    )
    reasons = [
        "S6 전체 전압·용량 재현성이 우수",
        "2C 충·방전에서 더 안정적인 전압 오차",
        "S7a fitting 시 parameter 보상이 작음",
        "공개된 absolute stoichiometry x 축 확보",
    ]
    y_start = 0.48
    for index, reason in enumerate(reasons):
        ax_decision.text(
            0.01,
            y_start - 0.105 * index,
            f"•  {reason}",
            fontsize=10.8,
            color="#26323F",
            transform=ax_decision.transAxes,
        )
    ax_decision.text(
        0.0,
        0.02,
        "저율 OCP는 0.5–1C 충전에서 경쟁력 있으나,\n2C 열세·더 큰 parameter 보상·OCP 상한 외삽 때문에 보류",
        fontsize=9.5,
        color=MUTED,
        transform=ax_decision.transAxes,
        linespacing=1.4,
    )

    figure.savefig(OUTPUT, dpi=220, facecolor="white", bbox_inches="tight")
    print(OUTPUT)


if __name__ == "__main__":
    main()
