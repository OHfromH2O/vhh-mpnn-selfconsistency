#!/usr/bin/env python3
"""델타 B — 구조 그림 3종을 PyMOL 스크립트로 생성하고 헤드리스 렌더링한다.

수작업 스크린샷이 아니라 **데이터에서 스크립트를 만들어** 렌더링합니다.
생성된 .pml은 figures/scripts/ 에 남고 git에 커밋되므로, 그림이 어떤
선택을 거쳐 만들어졌는지 제3자가 읽어서 검증할 수 있습니다.

실행 (WSL2 / .venv-linux):
    python scripts/06_render.py              # 그림 3종
    python scripts/06_render.py --only 3     # Fig 3만
    python scripts/06_render.py --design M1_T0.1_s7   # Fig 3의 대상 교체
    python scripts/06_render.py --dry-run    # .pml만 만들고 렌더링은 생략

주의: PyMOL은 네이티브 Windows에서 돌지 않습니다(PyPI 휠이 의존 DLL을
번들하지 않음). WSL2 또는 Colab에서 실행하십시오.
"""
from __future__ import annotations

import argparse
import io
import json
import pathlib
import shutil
import subprocess
import sys

import pandas as pd

ROOT = pathlib.Path(__file__).resolve().parent.parent
FIGDIR = ROOT / "figures"
PMLDIR = FIGDIR / "scripts"

NATIVE = ROOT / "work/01_clean/8YVJ_clean.pdb"
S1DIR = ROOT / "work/03_refold/S1_esmfold"
DESIGNABLE = ROOT / "work/01_clean/designable.json"
METRICS = ROOT / "work/04_metrics/metrics.csv"

# 렌더 공통 설정. ray=1이므로 디스플레이가 없어도 됩니다.
#
# auto_zoom을 끄는 이유: 켜져 있으면 load 할 때마다 카메라가 리셋되어,
# 스크립트 끝에서 지정한 orient/zoom이 중간 load에 덮일 수 있습니다.
#
# ★ load 후에 hide everything을 호출해야 합니다 ★
#   PyMOL은 polymer를 로드하면 cartoon을 자동으로 켜고(auto_show_classified)
#   객체마다 다른 색을 자동 배정합니다. hide를 load보다 먼저 쓰면
#   그 뒤에 로드된 객체들이 전부 기본 표시로 떠버립니다.
PREAMBLE = """# 자동 생성 — scripts/06_render.py. 직접 수정하지 마십시오.
reinitialize
set ray_opaque_background, 0
set cartoon_transparency, 0
set ray_shadows, 0
set antialias, 2
set cartoon_fancy_helices, 1
set auto_zoom, off
bg_color white

# 투명 표면을 ray로 렌더링하면 잘린 안쪽 면이 기본값(검정)으로 칠해져
# 계면이 검은 덩어리로 덮입니다. 내부 색을 밝게 주고 양면 조명을 켜면
# 표면 안쪽이 자연스럽게 보입니다.
set ray_interior_color, grey80
set ray_interior_texture, 0
set two_sided_lighting, on
set transparency_mode, 2
"""


def load_regions() -> dict:
    """designable.json에서 재설계 영역을 읽는다.

    그림에 쓰는 잔기 번호를 코드에 박지 않고 산출물에서 읽습니다.
    정의가 바뀌면 그림도 자동으로 따라갑니다.
    """
    d = json.load(io.open(DESIGNABLE, encoding="utf-8"))
    lo, hi = d["cdr3_span_1based"]
    cdr3 = list(range(lo, hi + 1))
    designable = d["designable_1based"]
    # 접촉 잔기 = 재설계 영역에서 CDR3를 뺀 나머지 (designable = 접촉 ∪ CDR3)
    contacts = [r for r in designable if r not in cdr3]
    return {
        "chain": d["design_chain"],
        "partner": d["partner_chain"],
        "cdr3": cdr3,
        "contacts": contacts,
        "cys": d["cysteines_in_designable"],
    }


def sel(resnums: list[int]) -> str:
    """잔기 번호 리스트를 PyMOL selection 문자열로."""
    return "+".join(str(r) for r in resnums)


def pick_design(mode: str, explicit: str | None) -> tuple[str, float]:
    """Fig 3·5에 쓸 대표 설계를 metrics.csv에서 결정론적으로 고른다.

    기본은 **중앙값**입니다. 최저 scRMSD를 쓰면 '가장 잘된 것'만 보여주는
    선택 편향이 되고, 이 레포는 한계 5에서 이미 S2 선정 편향을 지적한 바
    있습니다. 같은 기준을 그림에도 적용합니다.
    """
    if explicit:
        return explicit, float("nan")
    m = pd.read_csv(METRICS)
    s1 = m[(m["stage"] == "S1") & (m["group"].str.startswith("design"))]
    s1 = s1.sort_values("scRMSD").reset_index(drop=True)
    if mode == "best":
        row = s1.iloc[0]
    elif mode == "worst":
        row = s1.iloc[-1]
    else:  # median
        row = s1.iloc[len(s1) // 2]
    return row["name"], float(row["scRMSD"])


# ---------------------------------------------------------------------
# Fig 3 — native vs 설계 superposition
# ---------------------------------------------------------------------
def pml_fig3(reg: dict, design: str, scrmsd: float) -> str:
    design_pdb = S1DIR / f"{design}.pdb"
    return PREAMBLE + f"""
# Fig 3 — native 사슬 A vs 설계 {design} (ESMFold S1 예측)
# 중첩 기준: 전체 Cα. 셀 [40]의 scRMSD와 같은 정의라 그림과 수치가 대응합니다.
# (참고: 프레임워크만으로 맞추면 CDR3 차이가 더 커 보이지만, 보고한 수치와
#  어긋납니다. Fig 5는 루프 비교가 목적이라 프레임워크 기준을 씁니다.)
load {NATIVE}, native
load {design_pdb}, design

# 항원 사슬은 이 그림의 주제가 아니므로 제거
remove native and not chain {reg['chain']}

hide everything
show cartoon
color grey70, native
color skyblue, design
color orange, native and resi {sel(reg['cdr3'])}
color salmon, design and resi {sel(reg['cdr3'])}

# 재설계 영역에 포함된 Cys — 나노바디 CDR1-CDR3 추가 이황화결합의 파트너.
# 보존 여부가 README 결과 (5)의 근거이므로 눈으로도 확인할 수 있게 둡니다.
show sticks, (native or design) and resi {sel(reg['cys'])} and not hydro
color yellow, (native or design) and resi {sel(reg['cys'])} and elem S

# scRMSD = {scrmsd:.2f} A (셀 [40] metrics.csv 기준)
super design, native
orient native
turn y, 20

png {FIGDIR}/fig3_superposition.png, width=1800, height=1350, dpi=300, ray=1
"""


# ---------------------------------------------------------------------
# Fig 4 — 계면 접촉 잔기
# ---------------------------------------------------------------------
def pml_fig4(reg: dict) -> str:
    return PREAMBLE + f"""
# Fig 4 — native 복합체의 계면. 접촉 정의는 designable.json과 동일한
# 4.5 A 기준이며, 잔기 번호를 코드에 박지 않고 그 파일에서 읽었습니다.
load {NATIVE}, cx

hide everything
show cartoon, cx and chain {reg['chain']}
color grey80, cx and chain {reg['chain']}

# 항원은 표면으로 — 포켓 형상과 리간드 사슬의 맞물림을 보기 위함
show surface, cx and chain {reg['partner']}
color palecyan, cx and chain {reg['partner']}
set transparency, 0.35

# 접촉 잔기(4.5 A 이내)와 CDR3를 구분해서 표시
show sticks, cx and chain {reg['chain']} and resi {sel(reg['contacts'])} and not hydro
color tv_orange, cx and chain {reg['chain']} and resi {sel(reg['contacts'])} and elem C
color salmon, cx and chain {reg['chain']} and resi {sel(reg['cdr3'])}

# 극성 접촉 — 수소결합 후보. 거리 3.5 A 기준은 PyMOL 기본값입니다.
distance polar_contacts, \\
    cx and chain {reg['chain']} and resi {sel(reg['contacts'])}, \\
    cx and chain {reg['partner']}, 3.5, mode=2
# 거리 라벨은 끕니다. 접촉이 15잔기라 숫자가 겹쳐 계면을 덮습니다.
# 개수와 위치만 보이면 되고, 정확한 거리는 Rosetta 쪽 hbonds_int로 셉니다.
hide labels, polar_contacts
color grey30, polar_contacts
set dash_width, 2.0
set dash_gap, 0.3

# 전체가 아니라 계면을 중심으로 잡습니다. 복합체 전체를 orient하면
# 항원의 beta-solenoid가 화면을 차지하고 계면은 작아집니다.
orient cx and chain {reg['chain']} and resi {sel(reg['contacts'])}
zoom cx and chain {reg['chain']} and resi {sel(reg['contacts'])}, 16
turn x, -10
png {FIGDIR}/fig4_interface_contacts.png, width=1800, height=1350, dpi=300, ray=1
"""


# ---------------------------------------------------------------------
# Fig 5 — CDR3 backbone 앙상블
# ---------------------------------------------------------------------
TEMP_COLORS = {"T0.1": "marine", "T0.2": "forest", "T0.3": "firebrick"}


def pml_fig5(reg: dict, per_temp: int = 4) -> str:
    """온도별 설계의 CDR3 백본을 겹쳐 그린다.

    README 결과 (2) "샘플링 온도의 효과는 scRMSD에서만 보인다"에 대한
    시각적 보강입니다. 온도가 올라가면 CDR3 구조 다양성이 커지는지를
    눈으로 확인할 수 있습니다.
    """
    m = pd.read_csv(METRICS)
    s1 = m[(m["stage"] == "S1") & (m["group"].str.startswith("design"))]

    lines = [PREAMBLE, f"""
# Fig 5 — 온도별 CDR3 백본 앙상블 (ESMFold S1 예측)
# 중첩 기준: **프레임워크(CDR3 제외) Cα**. 루프끼리 비교하는 것이 목적이므로
# 전체 Cα로 맞추면 루프 차이가 프레임워크 쪽으로 흡수되어 보이지 않습니다.
# Fig 3와 기준이 다른 이유가 이것입니다.
#
# 순서 주의: load -> super 를 모두 끝낸 뒤에 hide everything 을 합니다.
# 먼저 숨기면 이후 로드된 객체가 기본 cartoon으로 되살아납니다.
load {NATIVE}, native
remove native and not chain {reg['chain']}
"""]

    # 온도별로 scRMSD 중앙값 근처를 고르게 뽑습니다. 상위만 뽑으면
    # "가장 잘된 것들"의 앙상블이 되어 다양성을 과소평가합니다.
    shown: list[tuple[str, str]] = []
    for temp, color in TEMP_COLORS.items():
        grp = s1[s1["group"] == f"design_{temp}"].sort_values("scRMSD")
        if grp.empty:
            continue
        step = max(1, len(grp) // per_temp)
        picked = grp.iloc[::step].head(per_temp)
        for _, r in picked.iterrows():
            obj = r["name"].replace(".", "_")
            shown.append((obj, color))
            lines.append(f"""
load {S1DIR / (r['name'] + '.pdb')}, {obj}
super {obj} and not resi {sel(reg['cdr3'])}, native and not resi {sel(reg['cdr3'])}""")

    # --- 여기서부터 표시. 모든 load가 끝난 뒤여야 합니다 ---------------
    lines.append("\n\nhide everything\n")
    for obj, color in shown:
        lines.append(f"show ribbon, {obj} and resi {sel(reg['cdr3'])}\n"
                     f"color {color}, {obj} and resi {sel(reg['cdr3'])}\n")

    lines.append(f"""
# native: 프레임워크는 흐린 cartoon으로 맥락만, CDR3는 검은 굵은 선으로 기준선
show cartoon, native and not resi {sel(reg['cdr3'])}
color grey80, native and not resi {sel(reg['cdr3'])}
set cartoon_transparency, 0.75, native
show ribbon, native and resi {sel(reg['cdr3'])}
color black, native and resi {sel(reg['cdr3'])}
set ribbon_width, 5
set ribbon_sampling, 2

orient native and resi {sel(reg['cdr3'])}
zoom native and resi {sel(reg['cdr3'])}, 8
png {FIGDIR}/fig5_cdr3_ensemble.png, width=1800, height=1350, dpi=300, ray=1
""")
    return "".join(lines)


# ---------------------------------------------------------------------
def run_pymol(pml: pathlib.Path) -> None:
    """헤드리스로 .pml을 실행한다. -c = GUI 없음, -q = 배너 억제."""
    exe = shutil.which("pymol")
    cmd = [exe, "-cq", str(pml)] if exe else [sys.executable, "-m", "pymol", "-cq", str(pml)]
    r = subprocess.run(cmd, capture_output=True, text=True)
    if r.returncode != 0:
        print(r.stdout[-2000:], file=sys.stderr)
        print(r.stderr[-2000:], file=sys.stderr)
        raise SystemExit(f"PyMOL 실패: {pml.name}")
    # PyMOL은 실패해도 종료 코드가 0인 경우가 있어 경고를 따로 본다
    for line in r.stdout.splitlines():
        if "Error" in line or "error:" in line:
            print(f"  [warn] {line}")


def main() -> int:
    ap = argparse.ArgumentParser(description="델타 B — 구조 그림 렌더링")
    ap.add_argument("--only", type=int, choices=[3, 4, 5], help="그림 하나만")
    ap.add_argument("--design", help="Fig 3·5의 대표 설계를 직접 지정")
    ap.add_argument("--pick", default="median", choices=["median", "best", "worst"],
                    help="대표 설계 선택 기준 (기본: median — 선택 편향 회피)")
    ap.add_argument("--dry-run", action="store_true", help=".pml만 생성")
    args = ap.parse_args()

    PMLDIR.mkdir(parents=True, exist_ok=True)
    FIGDIR.mkdir(exist_ok=True)

    reg = load_regions()
    design, scrmsd = pick_design(args.pick, args.design)
    print(f"대표 설계: {design}  (기준: {args.pick}, scRMSD={scrmsd:.2f} A)")
    print(f"CDR3: {reg['cdr3'][0]}-{reg['cdr3'][-1]}  |  접촉 {len(reg['contacts'])}잔기")

    jobs = {
        3: ("fig3_superposition.pml", lambda: pml_fig3(reg, design, scrmsd)),
        4: ("fig4_interface_contacts.pml", lambda: pml_fig4(reg)),
        5: ("fig5_cdr3_ensemble.pml", lambda: pml_fig5(reg)),
    }
    for n, (fname, builder) in jobs.items():
        if args.only and args.only != n:
            continue
        pml = PMLDIR / fname
        pml.write_text(builder(), encoding="utf-8")
        print(f"  [{n}] {pml.relative_to(ROOT)}", end="")
        if args.dry_run:
            print("  (dry-run)")
            continue
        run_pymol(pml)
        out = FIGDIR / fname.replace(".pml", ".png")
        size = out.stat().st_size // 1024 if out.exists() else 0
        print(f"  ->  {out.name}  ({size} KB)")

    return 0


if __name__ == "__main__":
    sys.exit(main())
