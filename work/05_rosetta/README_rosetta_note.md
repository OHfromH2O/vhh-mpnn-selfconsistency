# 델타 A — PyRosetta 계면 점수 (탐색적 스테이지)

원 프로토콜(`protein_design_sc_protocol.md`, 2026-10-22)에 없던 스테이지입니다.
사후 추가이므로 **탐색적(exploratory)** 으로 표기하며, 원 가설 H1의 검정에는 쓰지 않습니다.

---

## 왜 추가했는가

이 레포는 계면에 대해 아무 주장도 하지 못한 상태였습니다.

> 양성 대조군인 native 서열의 ipTM이 **0.12**입니다. 실제로 결정화된 복합체인데
> 사전 등록 임계값 0.5의 1/4 수준입니다. (README 결과 6)

AlphaFold-Multimer가 실패했으므로 **그 예측기에 의존하지 않는 독립 경로**로
계면을 평가합니다. Rosetta는 학습된 모델이 아니라 물리·통계 기반 scoring function이므로,
ESMFold·AF2와 오차의 원천이 다릅니다.

---

## 입력을 어떻게 골랐는가

| 후보 | 채택 | 사유 |
|---|---|---|
| S1 ESMFold 단량체 65건 | ❌ | 항원 사슬이 없어 계면 자체가 존재하지 않음 |
| S2 AF2-multimer 복합체 13건 | ❌ | **ipTM 0.12로 신뢰할 수 없다고 이미 보고한 구조.** 이것을 ddG 입력으로 쓰면 자기모순 |
| **native 결정 복합체 + 서열 threading** | ✅ | 백본은 1.65 Å 실험 구조. inverse folding은 길이를 보존하므로 잔기 대응이 1:1 |

`work/01_clean/8YVJ_clean.pdb`의 사슬 A(VHH, 129 aa)를 설계 서열로 치환하고,
사슬 B(Toxin A CROPs 단편, 121 aa)는 그대로 둡니다.

---

## 이 스테이지가 가정하는 것

**설계 서열이 native 백본을 유지한다.**

고정 백본 threading이므로, 설계가 실제로 유발할 수 있는 백본 변화는 반영되지 않습니다.
이 가정은 이미 측정한 scRMSD가 간접적으로 지지합니다(S1 중앙값 1.4 Å).
즉 두 지표가 서로를 보강하는 구조입니다.

- scRMSD: "설계 서열이 이 백본으로 되돌아오는가" — 예측기 기반
- 계면 점수: "그 백본 위에서 계면이 성립하는가" — 물리 기반

---

## 무엇을 주장하지 않는가

1. **결합 친화도를 예측하지 않습니다.** `InterfaceAnalyzer`의 `dG_separated`는
   실험 ΔΔG와의 상관이 중간 수준으로 보고되어 있습니다.
   후보 **순위화의 보조 지표**이지, 결합 여부의 판정 근거가 아닙니다.
2. **설계 간 미세 순위를 주장하지 않습니다.** 단일 타깃·동일 백본·서열만 다른
   65건은 동적 범위가 좁습니다. 대조군과의 분리만 봅니다.
3. **실험 검증이 아닙니다.** 레포 전체의 한계 1번이 그대로 적용됩니다.

---

## 좌표계 주의 (중요)

`FastRelax`는 구조를 **변형합니다.** 따라서:

- 기존 `work/04_metrics/metrics.csv`의 **scRMSD는 재계산하지 않습니다.**
  그 값은 이완 전 좌표에서 나온 것이며, 두 지표는 서로 다른 좌표계에 속합니다.
- 이완 전후 구조를 모두 보존합니다 (`relaxed/*.pdb`).
- 두 지표를 같은 표에 놓고 비교할 때 이 사실을 각주로 달아야 합니다.

---

## 확률성

`FastRelax`는 몬테카를로 기반이라 같은 입력에서도 결과가 흔들립니다.
설계당 `N_REPEAT`회 반복 후 **중앙값**을 씁니다(평균이 아닌 이유는 기존 표들이
전부 중앙값 통계이기 때문입니다).

반복 간 산포는 `rosetta_interface_raw.csv`에 그대로 남겨 두어,
독자가 "이 차이가 노이즈보다 큰가"를 직접 판단할 수 있게 합니다.

---

## 라이선스 고지 (독자용)

**PyRosetta는 학술·비영리 사용만 무료입니다.**

- 상업적 사용은 University of Washington CoMotion에서 별도 유료 라이선스가 필요합니다.
- **학술 사용자의 fee-for-service 작업도 상업적 사용으로 분류됩니다.**
- 조건은 소속이 아니라 **용도**를 기준으로 쓰여 있습니다.

출처: <https://www.pyrosetta.org/home/licensing-pyrosetta>

본 레포는 공개 교육 자료로서 비상업적 용도에 해당한다고 판단하여 사용했습니다.
**이 판단은 본 저자의 것이며 독자에게 이전되지 않습니다.**
영리 목적이거나 판단이 애매하면 `license@uw.edu`에 직접 문의하십시오.

PyRosetta를 쓸 수 없는 독자를 위해, 이 스테이지의 기하 지표 일부
(매몰 표면적·계면 접촉 수)는 Biopython만으로도 근사할 수 있습니다.
다만 그 경우 `dG_separated`와 `packstat`은 얻을 수 없습니다.

---

## 실행 환경 (검증본)

PyRosetta·PyMOL 모두 네이티브 Windows를 지원하지 않습니다.
이 스테이지는 WSL2에서 실행했습니다.

| | |
|---|---|
| OS | Ubuntu 22.04 (WSL2) |
| Python | 3.10.12 (`.venv-linux`) |
| PyRosetta | 2026.39 (`PyRosetta4.Release.python310.ubuntu`) |
| pymol-open-source | 3.2.0a0 (델타 B용) |
| 검증 일자 | 2026-11-02 |

**흔한 함정:** conda `base`가 활성인 채 `python3 -m venv`를 하면 venv가
conda의 Python으로 만들어집니다. 본 환경에서는 base가 3.7이라
`pymol-open-source`(3.9 이상)의 설치 후보가 아예 없다는 오류가 났습니다.
`conda deactivate` 후 `/usr/bin/python3`를 명시해야 합니다.

---

## 산출물

| 파일 | 내용 |
|---|---|
| `rosetta_interface.csv` | 설계별 중앙값 + provenance |
| `rosetta_interface_raw.csv` | 반복별 원시값 (산포 확인용) |
| `relaxed/*.pdb` | 이완 후 복합체 — `.gitignore` 대상 (재생성 가능) |
