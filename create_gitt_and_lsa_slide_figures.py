"""Create presentation-ready Korean summary figures for GITT Ds and LSA."""

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import font_manager
from matplotlib.patches import FancyBboxPatch, Rectangle
import numpy as np


ROOT = Path(__file__).resolve().parent
OUT = ROOT / "results" / "260927_presentation_figures"
OUT.mkdir(parents=True, exist_ok=True)

FONT_PATH = Path(r"C:\Windows\Fonts\malgun.ttf")
FONT = font_manager.FontProperties(fname=str(FONT_PATH)) if FONT_PATH.exists() else None
FONT_BOLD_PATH = Path(r"C:\Windows\Fonts\malgunbd.ttf")
FONT_BOLD = (
    font_manager.FontProperties(fname=str(FONT_BOLD_PATH))
    if FONT_BOLD_PATH.exists()
    else FONT
)

NAVY = "#123B63"
BLUE = "#4F81BD"
TEAL = "#0B6B8A"
LIGHT = "#EAF1F8"
MID = "#C9D6E4"
TEXT = "#18222D"
MUTED = "#5B6570"
ORANGE = "#E67E22"
RED = "#C43D3D"
GREEN = "#2E7D5B"
WHITE = "#FFFFFF"


def tx(ax, x, y, text, size=15, color=TEXT, weight="normal", ha="left", va="center"):
    prop = FONT_BOLD if weight == "bold" else FONT
    ax.text(
        x,
        y,
        text,
        fontsize=size,
        color=color,
        fontproperties=prop,
        ha=ha,
        va=va,
        transform=ax.transAxes,
    )


def rounded(ax, xy, width, height, face, edge="none", radius=0.012, lw=1.0):
    box = FancyBboxPatch(
        xy,
        width,
        height,
        boxstyle=f"round,pad=0.008,rounding_size={radius}",
        transform=ax.transAxes,
        facecolor=face,
        edgecolor=edge,
        linewidth=lw,
    )
    ax.add_patch(box)
    return box


def slide_header(ax, title, subtitle, page):
    tx(ax, 0.045, 0.935, title, 31, TEXT, "bold")
    ax.add_patch(Rectangle((0.035, 0.885), 0.31, 0.004, transform=ax.transAxes, color="#233A9F"))
    ax.add_patch(Rectangle((0.345, 0.885), 0.61, 0.0015, transform=ax.transAxes, color="#7C838B"))
    tx(ax, 0.035, 0.845, subtitle, 20, TEAL, "bold")
    tx(ax, 0.955, 0.035, page, 12, MUTED, ha="right")


def create_gitt_slide():
    fig, ax = plt.subplots(figsize=(16, 9), dpi=120)
    fig.patch.set_facecolor(WHITE)
    ax.set_axis_off()
    slide_header(ax, "GITT 기반 Ds 파라미터", "Diffusion coefficient (Dsp, Dsn)", "9/27")

    # Formula and method flow.
    tx(
        ax,
        0.055,
        0.735,
        r"$D_{GITT}=\dfrac{4}{\pi\tau}\left(\dfrac{m_BV_M}{M_BS}\right)^2\left(\dfrac{\Delta E_s}{\Delta E_\tau}\right)^2$",
        22,
    )
    tx(ax, 0.055, 0.600, "pulse별 D 계산", 13.5, NAVY, "bold")
    tx(ax, 0.175, 0.600, "→", 17, MUTED, "bold", ha="center")
    tx(ax, 0.205, 0.600, "방향·셀별 log 평균", 13.5, NAVY, "bold")
    tx(ax, 0.375, 0.600, "→", 17, MUTED, "bold", ha="center")
    tx(ax, 0.405, 0.600, "전체 구간 log 대표값", 13.5, NAVY, "bold")

    # Input assumptions table.
    tx(ax, 0.62, 0.79, "계산 입력 및 가정", 16, NAVY, "bold")
    rows = [
        ("항목", "Anode", "Cathode", "근거"),
        ("mB", "0.01040 g", "0.01725 g", "단면 coating composite mass"),
        ("MB", "72.066", "97.87 g/mol", "C6, LiCoO2 화학식"),
        ("VM", "31.89", "19.38 cm³/mol", "MB / ρ"),
        ("S", "0.7854", "0.7854 cm²", "10 mm disk geometry"),
        ("τ", "600", "600 s", "GITT CC pulse"),
    ]
    x0, y0, w, h = 0.60, 0.60, 0.355, 0.16
    col = [0.055, 0.075, 0.085, 0.140]
    rh = h / len(rows)
    for i, row in enumerate(rows):
        yy = y0 + h - (i + 1) * rh
        ax.add_patch(
            Rectangle(
                (x0, yy),
                w,
                rh,
                transform=ax.transAxes,
                facecolor=LIGHT if i == 0 else WHITE,
                edgecolor=MID,
                linewidth=0.7,
            )
        )
        xx = x0
        for j, value in enumerate(row):
            cx = xx + col[j] / 2
            tx(ax, cx, yy + rh / 2, value, 10.5, NAVY if i == 0 else TEXT, "bold" if i == 0 else "normal", ha="center")
            xx += col[j]

    # Main result panel.
    rounded(ax, (0.045, 0.235), 0.575, 0.305, "#F5F8FC", MID, 0.012, 1.0)
    tx(ax, 0.065, 0.505, "대표값 산정 결과", 18, NAVY, "bold")
    tx(ax, 0.595, 0.505, r"단위: $10^{-14}\,\mathrm{m^2/s}$", 11, MUTED, ha="right")
    result_rows = [
        ("전극", "0% Charge", "0% DChg", "0% 결합", "10% 비교", "Fitting init", "Bound"),
        ("Anode Dsn", "6.79", "2.88", "4.42", "2.12–2.53", "4.42", "1.0–8.0"),
        ("Cathode Dsp", "4.56", "10.17", "6.81", "4.35–4.37", "6.81", "2.0–12.0"),
    ]
    rx, ry, rw, rrh = 0.065, 0.335, 0.535, 0.055
    widths = [0.105, 0.072, 0.078, 0.078, 0.078, 0.078, 0.068]
    for i, row in enumerate(result_rows):
        yy = ry + (len(result_rows) - 1 - i) * rrh
        xx = rx
        for j, value in enumerate(row):
            face = BLUE if i == 0 else ("#E0E8F2" if i == 1 else "#EDF2F7")
            ax.add_patch(Rectangle((xx, yy), widths[j], rrh, transform=ax.transAxes, facecolor=face, edgecolor=WHITE, linewidth=1.1))
            color = WHITE if i == 0 else (TEAL if j in (3, 5) else TEXT)
            tx(ax, xx + widths[j] / 2, yy + rrh / 2, value, 11.2, color, "bold" if i == 0 or j in (0, 3, 5) else "normal", ha="center")
            xx += widths[j]
    tx(ax, 0.065, 0.285, "0% 결합값 = √(Charge × Discharge)", 12.5, MUTED)
    tx(ax, 0.595, 0.285, "10% 값은 강건성 비교로 유지", 12.5, MUTED, ha="right")

    # Interpretation / model use panel.
    rounded(ax, (0.65, 0.235), 0.305, 0.305, "#F7FAFC", MID, 0.012, 1.0)
    tx(ax, 0.67, 0.505, "DFN 모델 적용", 18, NAVY, "bold")
    rounded(ax, (0.675, 0.42), 0.25, 0.055, LIGHT)
    tx(ax, 0.80, 0.447, "all-range apparent GITT reference", 12.2, NAVY, "bold", ha="center")
    tx(ax, 0.80, 0.39, "고정 물성값이 아니라 초기값·bound 근거", 12, TEXT, ha="center")
    rounded(ax, (0.675, 0.305), 0.115, 0.055, "#DDEFE8")
    rounded(ax, (0.81, 0.305), 0.115, 0.055, "#F8E7D5")
    tx(ax, 0.7325, 0.332, "Dsn: fitting", 12.5, GREEN, "bold", ha="center")
    tx(ax, 0.8675, 0.332, "Dsp: 6.81 고정", 12.5, ORANGE, "bold", ha="center")
    tx(ax, 0.80, 0.27, "향후 pulse + relaxation direct fitting으로 재검증", 11.5, MUTED, ha="center")

    # Caveat strip.
    rounded(ax, (0.045, 0.115), 0.91, 0.075, "#FFF7EA", "#E8C48A", 0.01, 1.0)
    tx(ax, 0.065, 0.152, "해석 한계", 13, ORANGE, "bold")
    tx(
        ax,
        0.15,
        0.152,
        "0%는 투명한 raw reference일 뿐 intrinsic D가 아님  |  600 s가 R²/D와 유사  |  composite mass를 100 wt% active로 가정",
        11.5,
        TEXT,
    )
    fig.savefig(OUT / "gitt_ds_slide_layout.png", dpi=160, facecolor=WHITE)
    plt.close(fig)


def create_fitting_design_slide():
    fig, ax = plt.subplots(figsize=(16, 9), dpi=120)
    fig.patch.set_facecolor(WHITE)
    ax.set_axis_off()
    slide_header(ax, "Charge fitting 설계", "0.5C / 1C / 2C voltage-only objective", "11/27")

    # Objective.
    rounded(ax, (0.045, 0.685), 0.91, 0.14, "#F5F8FC", MID, 0.012, 1.0)
    tx(ax, 0.065, 0.785, "목적함수", 17, NAVY, "bold")
    tx(
        ax,
        0.235,
        0.775,
        r"$J_V(\theta)=\sqrt{\dfrac{1}{3N_{cell}N_t}\sum_{c\in\{0.5,1,2\}}\sum_j\sum_i"
        r"\left[1000\left(V_{DFN}(t_i;\theta)-V_{exp,j}(t_i)\right)\right]^2}$",
        15.5,
    )
    tx(ax, 0.235, 0.710, "각 C-rate·셀 동일 가중  |  voltage L2/RMSE만 최소화  |  capacity·prior·형상항 없음", 11.7, MUTED)

    # Primary model.
    rounded(ax, (0.045, 0.35), 0.435, 0.285, "#EEF6F1", "#B9D7C8", 0.012, 1.0)
    tx(ax, 0.065, 0.600, "주 모델", 17, GREEN, "bold")
    tx(ax, 0.265, 0.600, "Dsn + kn + brugg_n", 18, NAVY, "bold", ha="center")
    primary_rows = [
        ("Parameter", "Initial", "Bound", "scale"),
        ("Dsn [m²/s]", "4.42e-14", "1.0e-14 – 8.0e-14", "log"),
        ("kn", "9.65e-7", "3.0e-7 – 3.0e-6", "log"),
        ("brugg_n", "2.914", "1.5 – 3.5", "linear"),
    ]
    px, py, ph = 0.065, 0.405, 0.043
    pwidths = [0.12, 0.09, 0.145, 0.055]
    for i, row in enumerate(primary_rows):
        yy = py + (len(primary_rows) - 1 - i) * ph
        xx = px
        for j, value in enumerate(row):
            ax.add_patch(Rectangle((xx, yy), pwidths[j], ph, transform=ax.transAxes, facecolor=GREEN if i == 0 else WHITE, edgecolor=WHITE, linewidth=1.0))
            tx(ax, xx + pwidths[j] / 2, yy + ph / 2, value, 10.4, WHITE if i == 0 else TEXT, "bold" if i == 0 else "normal", ha="center")
            xx += pwidths[j]
    tx(ax, 0.265, 0.372, "목적: 확산 + 반응속도 + 전해질 수송을 함께 분리", 11.2, MUTED, ha="center")

    # Conservative model.
    rounded(ax, (0.52, 0.35), 0.435, 0.285, "#FFF8ED", "#E8CAA0", 0.012, 1.0)
    tx(ax, 0.54, 0.600, "보수 모델", 17, ORANGE, "bold")
    tx(ax, 0.745, 0.600, "Dsn + kn", 18, NAVY, "bold", ha="center")
    conservative_rows = [
        ("Parameter", "Initial", "Bound", "scale"),
        ("Dsn [m²/s]", "4.42e-14", "1.0e-14 – 8.0e-14", "log"),
        ("kn", "9.65e-7", "3.0e-7 – 3.0e-6", "log"),
        ("brugg_n", "2.914", "고정", "—"),
    ]
    px = 0.54
    for i, row in enumerate(conservative_rows):
        yy = py + (len(conservative_rows) - 1 - i) * ph
        xx = px
        for j, value in enumerate(row):
            ax.add_patch(Rectangle((xx, yy), pwidths[j], ph, transform=ax.transAxes, facecolor=ORANGE if i == 0 else WHITE, edgecolor=WHITE, linewidth=1.0))
            tx(ax, xx + pwidths[j] / 2, yy + ph / 2, value, 10.4, WHITE if i == 0 else TEXT, "bold" if i == 0 else "normal", ha="center")
            xx += pwidths[j]
    tx(ax, 0.745, 0.372, "목적: 상관 위험을 낮춘 최소 식별 subset", 11.2, MUTED, ha="center")

    # Optimisation ladder.
    tx(ax, 0.055, 0.290, "동일 알고리즘 사다리", 16, NAVY, "bold")
    steps = [
        (0.055, "1", "Multi-start TRF", "raw/10%/중앙 시작점"),
        (0.275, "2", "Differential Evolution", "global → TRF polish"),
        (0.515, "3", "Dual Annealing", "global comparator → TRF"),
        (0.755, "4", "Discharge validation", "capacity는 사후 평가"),
    ]
    for i, (x, number, title, subtitle) in enumerate(steps):
        rounded(ax, (x, 0.155), 0.19, 0.095, LIGHT if i < 3 else "#DDEFE8", "none", 0.012, 0.0)
        rounded(ax, (x + 0.012, 0.178), 0.035, 0.045, BLUE if i < 3 else GREEN, "none", 0.018, 0.0)
        tx(ax, x + 0.0295, 0.200, number, 12, WHITE, "bold", ha="center")
        tx(ax, x + 0.055, 0.215, title, 12.2, NAVY, "bold")
        tx(ax, x + 0.055, 0.178, subtitle, 10.5, MUTED)
        if i < len(steps) - 1:
            tx(ax, x + 0.207, 0.202, "→", 18, MUTED, "bold", ha="center")

    tx(ax, 0.055, 0.095, "선택 규칙", 13, NAVY, "bold")
    tx(ax, 0.145, 0.095, "charge RMSE 단독 최저가 아니라, bound 비접촉 + 반복 수렴 + discharge 일반화가 확인된 해를 선택", 11.5, TEXT)
    tx(ax, 0.055, 0.055, "고정: Dsp=6.81e-14, kp=5.90e-7, brugg_p=1.83, brugg_s=1.5  |  조기 cutoff는 time-grid 미충족으로 실패 처리", 11.2, MUTED)

    fig.savefig(OUT / "charge_fitting_design_slide.png", dpi=160, facecolor=WHITE)
    plt.close(fig)


def create_fitting_bound_basis_slide():
    fig, ax = plt.subplots(figsize=(16, 9), dpi=120)
    fig.patch.set_facecolor(WHITE)
    ax.set_axis_off()
    slide_header(ax, "Fitting bound 설정 근거", "측정 envelope · 원논문 nominal · porous-electrode transport", "12/27")

    # Header row.
    cols = [0.16, 0.43, 0.18, 0.16]
    x0, y_top, row_h = 0.045, 0.76, 0.165
    headers = ["Parameter", "근거 데이터와 산정", "1차 fitting bound", "경계 도달 시"]
    xx = x0
    for width, label in zip(cols, headers):
        ax.add_patch(Rectangle((xx, y_top), width, 0.055, transform=ax.transAxes, facecolor=BLUE, edgecolor=WHITE, linewidth=1.2))
        tx(ax, xx + width / 2, y_top + 0.0275, label, 12.5, WHITE, "bold", ha="center")
        xx += width

    rows = [
        {
            "name": "Dsn",
            "sub": "apparent GITT",
            "basis1": "방향+결합 GITT envelope: 1.689 - 6.792e-14 m2/s",
            "basis2": "charge-only fitting이므로 0% charge 6.792e-14도 포함",
            "bound": "1.0e-14 -\n8.0e-14 m2/s",
            "edge": "상한 근접 시\n1.2e-13까지 확장",
            "face": "#EEF6F1",
            "accent": GREEN,
        },
        {
            "name": "kn",
            "sub": "reaction rate",
            "basis1": "Ai2020 kref=1.0e-11 m/s -> F*kref=9.65e-7",
            "basis2": "nominal +/-0.5 log10 decade; EIS 2.48/7.40e-7은 비교점",
            "bound": "3.0e-7 -\n3.0e-6",
            "edge": "경계 근접 시\n2.0e-7 - 5.0e-6",
            "face": "#EEF3FA",
            "accent": BLUE,
        },
        {
            "name": "brugg_n",
            "sub": "electrolyte exponent",
            "basis1": "이론 하한 1.5 · graphite EIS 환산 약 2.73",
            "basis2": "Ai2020: 1+alpha_B = 2.914",
            "bound": "1.5 – 3.5",
            "edge": "3.5 근접 시\n4.0까지 확장",
            "face": "#FFF8ED",
            "accent": ORANGE,
        },
    ]
    for idx, row in enumerate(rows):
        yy = y_top - (idx + 1) * row_h
        xx = x0
        for width in cols:
            ax.add_patch(Rectangle((xx, yy), width, row_h, transform=ax.transAxes, facecolor=row["face"], edgecolor=WHITE, linewidth=1.2))
            xx += width
        tx(ax, x0 + cols[0] / 2, yy + 0.105, row["name"], 18, row["accent"], "bold", ha="center")
        tx(ax, x0 + cols[0] / 2, yy + 0.060, row["sub"], 10.7, MUTED, ha="center")
        bx = x0 + cols[0] + 0.02
        tx(ax, bx, yy + 0.108, row["basis1"], 12.2, TEXT, "bold")
        tx(ax, bx, yy + 0.060, row["basis2"], 11.5, MUTED)
        bound_x = x0 + cols[0] + cols[1] + cols[2] / 2
        edge_x = x0 + cols[0] + cols[1] + cols[2] + cols[3] / 2
        tx(ax, bound_x, yy + row_h / 2, row["bound"], 14.2, NAVY, "bold", ha="center")
        tx(ax, edge_x, yy + row_h / 2, row["edge"], 11.5, TEXT, ha="center")

    rounded(ax, (0.045, 0.145), 0.91, 0.09, "#F5F8FC", MID, 0.01, 1.0)
    tx(ax, 0.065, 0.190, "해석 원칙", 13, NAVY, "bold")
    tx(ax, 0.155, 0.190, "bound는 95% 신뢰구간이 아니라 물리·실험 anchor를 포함하는 탐색 범위", 12.2, TEXT)
    tx(ax, 0.155, 0.158, "최적값이 경계에 붙으면 물성 확정이 아니라 bound/식별도 문제로 판정", 11.5, MUTED)

    tx(ax, 0.055, 0.095, "근거", 12.5, NAVY, "bold")
    tx(ax, 0.115, 0.095, "Ai2020 Table I · GITT endpoint sensitivity audit · Landesfeind et al. (2016) graphite tortuosity", 11.5, MUTED)
    tx(ax, 0.055, 0.055, "수정: Dsn 상한 5.0e-14 -> 8.0e-14 (0% charge 포함)  |  brugg_n 상한 4.0 -> 3.5", 11.5, RED, "bold")

    fig.savefig(OUT / "fitting_bound_basis_slide.png", dpi=160, facecolor=WHITE)
    plt.close(fig)


def create_lsa_slide():
    params = ["kn", "kp", "brugg_n", "Dsn", "Dsp"]
    rel = np.array([1.0, 0.9825486, 0.7709933, 0.283908, 0.0234306])
    corr_labels = ["Dsn", "Dsp", "kn", "kp", "brugg_n"]
    corr = np.array(
        [
            [1.000, 0.237, 0.411, 0.408, 0.033],
            [0.237, 1.000, 0.336, 0.367, 0.495],
            [0.411, 0.336, 1.000, 0.996, 0.815],
            [0.408, 0.367, 0.996, 1.000, 0.814],
            [0.033, 0.495, 0.815, 0.814, 1.000],
        ]
    )

    fig = plt.figure(figsize=(16, 9), dpi=120, facecolor=WHITE)
    canvas = fig.add_axes([0, 0, 1, 1])
    canvas.set_axis_off()
    slide_header(canvas, "민감도 / 식별도 분석", "Charge 0.5C / 1C / 2C · mixed perturbation LSA", "10/27")

    # Sensitivity ranking.
    ax1 = fig.add_axes([0.075, 0.31, 0.405, 0.48])
    colors = [BLUE, "#E2A84A", BLUE, BLUE, RED]
    bars = ax1.bar(params, rel, color=colors, width=0.62)
    ax1.axhline(0.10, color=RED, ls="--", lw=1.6, label="cutoff = 0.10")
    ax1.set_ylim(0, 1.12)
    ax1.set_ylabel("Relative sensitivity norm", fontproperties=FONT, fontsize=12)
    ax1.set_title("민감도 순위", fontproperties=FONT_BOLD, fontsize=18, color=NAVY, pad=12)
    ax1.grid(axis="y", alpha=0.22)
    ax1.legend(prop=FONT, frameon=False, loc="upper right")
    ax1.tick_params(labelsize=11)
    for label in ax1.get_xticklabels() + ax1.get_yticklabels():
        label.set_fontproperties(FONT)
    for bar, value in zip(bars, rel):
        ax1.text(bar.get_x() + bar.get_width() / 2, value + 0.025, f"{value:.3f}", ha="center", va="bottom", fontsize=11, fontproperties=FONT_BOLD, color=TEXT)

    # Correlation matrix.
    ax2 = fig.add_axes([0.555, 0.31, 0.37, 0.48])
    im = ax2.imshow(corr, vmin=0, vmax=1, cmap="Blues")
    ax2.set_xticks(range(5), corr_labels, rotation=25, ha="right")
    ax2.set_yticks(range(5), corr_labels)
    ax2.set_title("절대 Pearson correlation", fontproperties=FONT_BOLD, fontsize=18, color=NAVY, pad=12)
    for label in ax2.get_xticklabels() + ax2.get_yticklabels():
        label.set_fontproperties(FONT)
        label.set_fontsize(11)
    for i in range(5):
        for j in range(5):
            ax2.text(j, i, f"{corr[i, j]:.2f}", ha="center", va="center", fontsize=10.5, color=WHITE if corr[i, j] > 0.62 else TEXT, fontproperties=FONT_BOLD)
    cax = fig.add_axes([0.925, 0.37, 0.012, 0.36])
    cb = fig.colorbar(im, cax=cax)
    cb.ax.tick_params(labelsize=9)

    # Decisions.
    rounded(canvas, (0.055, 0.135), 0.27, 0.105, "#DDEFE8")
    tx(canvas, 0.19, 0.205, "1차 fitting subset", 14, GREEN, "bold", ha="center")
    tx(canvas, 0.19, 0.165, "Dsn + kn + brugg_n", 18, NAVY, "bold", ha="center")

    rounded(canvas, (0.365, 0.135), 0.27, 0.105, "#FFF1DA")
    tx(canvas, 0.50, 0.205, "kp 제외", 14, ORANGE, "bold", ha="center")
    tx(canvas, 0.50, 0.165, "민감도 낮음이 아니라 kn과 r = 0.996", 12.5, TEXT, ha="center")

    rounded(canvas, (0.675, 0.135), 0.27, 0.105, "#F8E2E2")
    tx(canvas, 0.81, 0.205, "Dsp 제외", 14, RED, "bold", ha="center")
    tx(canvas, 0.81, 0.165, "상대 민감도 0.023 < cutoff 0.10", 12.5, TEXT, ha="center")

    tx(canvas, 0.055, 0.085, "Robustness: eps = 2/5/10%, grid = 100/250/500 · β = 0.90에서 9/9 동일 subset", 12, MUTED)
    fig.savefig(OUT / "lsa_identifiability_slide.png", dpi=160, facecolor=WHITE)
    plt.close(fig)


if __name__ == "__main__":
    create_gitt_slide()
    create_lsa_slide()
    create_fitting_design_slide()
    create_fitting_bound_basis_slide()
