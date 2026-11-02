# 그림

| 파일 | 내용 | 생성 |
|---|---|---|
| `fig1_scrmsd.png` | 그룹별 scRMSD 분포 | 셀 `[41]` (matplotlib) |
| `fig2_scrmsd_vs_plddt.png` | scRMSD vs pLDDT | 셀 `[41]` (matplotlib) |
| `fig3_superposition.png` | native vs 설계 중첩 | `scripts/06_render.py` (PyMOL) |
| `fig4_interface_contacts.png` | 계면 접촉 잔기 | `scripts/06_render.py` (PyMOL) |
| `fig5_cdr3_ensemble.png` | 온도별 CDR3 백본 앙상블 | `scripts/06_render.py` (PyMOL) |
| `fig6_rosetta_interface.png` | 계면 점수 분포·상관 | 셀 `[46]` (matplotlib) |

구조 그림(fig3–5)은 **수작업 스크린샷이 아닙니다.** `scripts/06_render.py`가
산출물에서 잔기 번호와 대상 설계를 읽어 `.pml` 스크립트를 생성하고,
PyMOL을 헤드리스로 돌려 렌더링합니다. 생성된 `.pml`은
[`figures/scripts/`](scripts/)에 커밋되어 있으므로 어떤 선택을 거쳤는지
제3자가 읽어서 확인할 수 있습니다.

```bash
# WSL2 / .venv-linux
python scripts/06_render.py                    # 3종 전부
python scripts/06_render.py --only 3           # 하나만
python scripts/06_render.py --pick best        # 대표 설계 기준 변경
python scripts/06_render.py --design M1_T0.1_s7  # 특정 설계 지정
python scripts/06_render.py --dry-run          # .pml만 생성
```

---

## 그림별 선택과 근거

### Fig 3 — native vs 설계 중첩

- **대상 설계는 scRMSD 중앙값**입니다(`M2_T0.3_s2`, 1.21 Å).
  최저값을 쓰면 "가장 잘된 것"만 보여주는 선택 편향이 됩니다.
  README 한계 5에서 S2 선정 편향을 이미 지적했으므로 같은 기준을 그림에도 적용했습니다.
  `--pick best` / `--pick worst`로 양 끝도 확인할 수 있습니다.
- **중첩 기준은 전체 Cα**로, 셀 `[40]`의 scRMSD 정의와 같습니다.
  그림과 보고한 수치가 대응하도록 맞췄습니다.
- 회색 = native, 파랑 = 설계. 주황/연분홍 = 각각의 CDR3.
- 노란 막대 = Cys110. 나노바디 CDR1–CDR3 추가 이황화결합의 파트너이며,
  README 결과 (5)의 보존율 집계 대상입니다.

### Fig 4 — 계면 접촉 잔기

- 접촉 잔기 정의(4.5 Å)와 번호는 `work/01_clean/designable.json`에서 읽습니다.
  코드에 박아 두지 않았으므로 정의가 바뀌면 그림도 따라갑니다.
- 연청록 표면 = 항원(사슬 B), 회색 cartoon = VHH(사슬 A).
- 주황 막대 = 접촉 잔기 15개, 연분홍 = CDR3.
- 회색 점선 = 3.5 Å 이내 극성 접촉. **거리 라벨은 껐습니다** — 15잔기라
  숫자가 겹쳐 계면을 덮습니다. 정확한 수소결합 수는 델타 A의 `hbonds_int`로 셉니다.

### Fig 5 — 온도별 CDR3 백본 앙상블

- README 결과 (2) "샘플링 온도의 효과는 scRMSD에서만 보인다"에 대한 시각적 보강입니다.
- 파랑 = T0.1, 초록 = T0.2, 빨강 = T0.3, 검정 = native. 온도당 4건.
- **중첩 기준이 Fig 3와 다릅니다 — 프레임워크(CDR3 제외) Cα입니다.**
  루프끼리 비교하는 것이 목적인데 전체 Cα로 맞추면 루프 차이가
  프레임워크 쪽으로 흡수되어 보이지 않습니다.
- 각 온도에서 scRMSD 순으로 등간격 추출합니다. 상위만 뽑으면
  다양성을 과소평가합니다.

---

## 한계

1. fig3–5는 **ESMFold S1 예측 구조**입니다. 실험 구조가 아닙니다.
   시각적 유사성은 scRMSD·TM-score를 대체하지 않습니다.
2. Fig 5의 앙상블은 온도당 4건으로, 분포가 아니라 예시입니다.
3. 그림은 단일 타깃(8YVJ) 단일 시점입니다. 회전각은 사람이 고른 값이며
   `.pml`의 `orient`/`turn`에 그대로 적혀 있습니다.

---

## 실행 환경

PyMOL은 **네이티브 Windows에서 돌지 않습니다.** PyPI 휠이 의존 DLL
(freetype, glew32, libpng16, libxml2, netcdf)을 번들하지 않아 `_cmd` import에서
실패합니다. WSL2 또는 Colab에서 실행하십시오.

검증본: pymol-open-source 3.2.0a0, Python 3.10.12, Ubuntu 22.04 (WSL2), 2026-11-02.
