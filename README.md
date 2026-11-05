# VHH inverse-folding self-consistency benchmark

**한국어** | [English](README.en.md)

ProteinMPNN으로 재설계한 나노바디 서열이, 독립적인 구조 예측기로 되돌렸을 때
원본 backbone을 재현하는가(self-consistency)를 음성·양성 대조군과 함께 측정한 벤치마크 실습.

타깃: PDB 8YVJ — H5.2 나노바디(VHH, 129 aa) / *C. difficile* toxin A CROPs 단편(118 aa), X-ray 1.65 Å

---

## 한눈에 보기

| 단계 | 도구 | 핵심 결과 | 주요 한계 |
|---|---|---|---|
| **타깃 선정** | RCSB Search API | 8YVJ — VHH 129 aa / 항원 118 aa, X-ray 1.65 Å. 후보 3건 중 항원이 펩타이드가 아닌 유일한 건 | 단일 타깃 (한계 6) |
| **전처리 · 영역 정의** | PDBFixer, OpenMM, Biopython | 재설계 영역 **34/129자리** = 인터페이스(4.5 Å) ∪ CDR3 | ANARCI 미사용 → 정통 CDR 정의 아님 (한계 7) |
| **서열 설계** | **ProteinMPNN** | 설계 48건 (2모드 × 3온도 × 8) + 대조군 17건 | — |
| **단량체 재예측** | **ESMFold** | scRMSD **설계 1.20 Å** vs 셔플 3.01 / 무작위 3.90 → **H1 지지** | 조건당 14건 (한계 3) |
| **복합체 재예측** | **ColabFold** (AF2-multimer v3) | **ipTM 측정 실패** — 양성 대조군 native가 0.12 (임계값 0.5) | 축소 실행 (한계 4), 선정 편향 (한계 5) |
| **계면 교차검증** *(탐색적)* | **PyRosetta** ref2015 · InterfaceAnalyzer | 설계 **42/42** dG<0, 셔플 **0/5**. native가 전체 1위 (−49.6 REU).<br>단 **매몰 미충족 극성 원자**는 대조군과 구별 안 됨 (결과 8) | 고정 백본 가정 (한계 11), REU는 kcal/mol 아님 (한계 12) |
| **시각화** | **PyMOL** (스크립트 생성 + 헤드리스) | Fig 3–5 — `.pml`까지 커밋해 재현 가능 | ESMFold 예측 구조 기반 |

→ 상세: [무엇을 했는가](#무엇을-했는가) · [결과](#결과) · [한계](#한계) · [다음 단계](#다음-단계)

### 세 줄 요약

1. 되는 것 — inverse folding 서열은 원본 backbone을 재현하고(1.20 Å),
   그 재현도는 조성을 보존한 음성 대조군(3.01 Å)과 겹치지 않습니다.
2. 안 되는 것 — 복합체 계면은 평가조차 못 했습니다. 양성 대조군이 먼저 실패했고,
   그래서 그 축은 결과 대신 측정 실패로 보고하였습니다.
3. 그래서 한 것 — 물리 기반 경로(Rosetta)로 계면 축을 일부 복구했습니다.
   기준은 실행 전에 등록했고(`preregistration.md` §0.4), 선후는 커밋 순서로 확인됩니다.
   그 과정에서 설계가 복원하지 못한 축 하나를 찾았습니다. 계면에 묻힌
   극성 원자의 수소결합 짝(결과 8). 이 축에서는 설계가 무작위 대조군과
   구별되지 않습니다.

---

## 이 문서를 읽는 법

결과보다 한계를 먼저 읽도록 썼습니다. [무엇을 하지 않았는가](#무엇을-하지-않았는가)가
[무엇을 했는가](#무엇을-했는가) 바로 뒤에 오고, 판정 기준은 코드를 실행하기 전에
[`preregistration.md`](preregistration.md)에 등록했습니다. 양성 대조군이 통과하지
못한 지표(ipTM, native 0.12)는 측정 실패로 보고하였고, 그 축에 대해서는
아무것도 주장하지 않습니다.

이렇게 쓰는 이유는 저자가 계산 결과를 받는 쪽에 있어 봤기 때문입니다.
석사 과정에서 클로닝·단백질 정제·LC/MS를 직접 수행했고, iGEM에서는 드라이랩의
후보 선별 결과를 웻랩 실험 우선순위로 넘기는 역할을 맡았습니다. 실험으로 넘어가는
순간 중요해지는 것은 순위 자체보다 그 순위가 무엇을 보장하지 않는가였습니다.
in silico self-consistency는 발현성이나 결합력의 대리 지표가 아니며,
이 레포는 그 경계를 지표마다 명시합니다.

---

## 무엇을 했는가

1. RCSB Search API로 나노바디 복합체 후보를 조회하고, 사슬 구성으로 2차 선별해 타깃을 확정했습니다.
   해상도 ≤ 2.5 Å · 총 ≤ 300잔기 조건을 만족한 3건 중, 항원이 펩타이드가 아닌 단백질인 것은 8YVJ뿐이었습니다
   (나머지 2건은 14 aa Tau PHF6 펩타이드 복합체). → `data/chain_composition.csv`
2. PDBFixer로 전처리하고, 재설계 영역을 항원과 4.5 Å 이내 접촉하는 잔기 ∪ CDR3 (34/129자리)로 정의했습니다.
   → `work/01_clean/designable.json`, 정의 변경 사유는 [`work/01_clean/README_definition_note.md`](work/01_clean/README_definition_note.md)
3. ProteinMPNN으로 2모드 × 3온도 × 8서열 = 48 설계 서열을 생성했습니다.
   - M1: VHH 전체 129자리 재설계 / M2: 위 34자리만 재설계(프레임워크 고정)
   - temperature 0.1 / 0.2 / 0.3
4. 대조군 17건: 조성 보존 셔플 8, 완전 무작위 치환 8 (모두 재설계 영역 한정), native 1(양성 대조군)
5. 2단계로 재예측했습니다.
   - S1: ESMFold로 VHH 단량체 65건 전부
   - S2: ColabFold(AF2-multimer v3)로 복합체 13건 (S1 pLDDT 상위 설계 6 + 대조군 6 + native 1)
6. scRMSD · TM-score · pLDDT(전체/재설계 영역) · ipTM · seq_recovery를 산출했습니다.

---

## 무엇을 하지 않았는가

다음은 본 산출물로 입증되지 않습니다.

- 설계 서열이 실제로 발현·폴딩된다는 것 (실험 검증 없음)
- 설계 서열이 항원에 결합한다는 것 (결합 친화도 미측정)
- 복합체 계면을 재현한다는 것 — ipTM 측정에 실패했습니다 (결과 6).
  델타 A(결과 7)는 복합체를 예측하지 않고 native 결합 양식을 전제한 뒤
  서열이 그 기하에 맞는지만 봅니다. 결합 양식 자체는 여전히 가정입니다.
- 통계적 유의성 (N이 작아 기술통계 수준의 비교만 수행)

확실히 주장할 수 있는 사실: 주어진 backbone에 대해 inverse folding으로 생성한 서열이,
독립적인 구조예측기로 되돌렸을 때 원본 backbone을 재현하며, 그 재현도가 음성 대조군보다
명확히 높다.

---

## 사전 등록된 기준

판정 기준은 [`protein_design_sc_protocol.md`](protein_design_sc_protocol.md)(2026-10-22) §0.3에 코드 실행 전부터 적혀 있었고,
결과를 본 뒤 완화하거나 변경하지 않았습니다. 다만 이를 담은 [`preregistration.md`](preregistration.md) 파일은
분석을 마친 뒤에 작성되었으므로, 선후를 증명할 커밋 기록은 없습니다.
그 사실과 "데이터를 본 뒤에 내린 결정" 목록을 [`preregistration.md`](preregistration.md)에 전부 적었습니다.

```
성공 판정 (설계 1건 단위)
  scRMSD (재예측 구조 vs 원본 backbone, Cα) < 2.0 Å
  재예측 pLDDT 평균 > 80
  (복합체 수행 시) ipTM > 0.5

H1: 설계 서열의 scRMSD 중앙값 < 셔플 대조군의 scRMSD 중앙값
H0: 차이 없음
```

기준은 사후에 변경하지 않았습니다. 아래 통과 집계는 위 기준을 그대로 적용한 결과입니다.

---

## 방법 (버전 고정)

| 항목 | 값 |
|---|---|
| 실행 환경 | Google Colab, Tesla T4 (15,360 MiB), Linux 6.6.122, Python 3.13.16 |
| 실행 시각 | 2026-10-26 05:52 UTC |
| seed | 42 (전 단계 고정) |
| ProteinMPNN | `dauparas/ProteinMPNN` (MIT License), `vanilla_model_weights/v_48_020.pt` |
| ESMFold | `facebook/esmfold_v1` (HuggingFace transformers), 언어모델 fp16 / 폴딩 트렁크 fp32 |
| ColabFold | `alphafold2_multimer_v3`, `--num-recycle 3 --num-models 1 --random-seed 42` |
| 전처리 | PDBFixer 1.x + OpenMM 8.6.1 (측쇄 보완, 말단 결측은 보완하지 않음) |
| 구조 비교 | Biopython 1.88 Cα 중첩 (US-align 미사용 — 아래 한계 참조) |
| 기타 | numpy 2.1.3, pandas 2.2.3, matplotlib 3.10.0, torch 2.11.0+cu130 |
| 계면 점수 (탐색적) | PyRosetta 2026.39 `python310.ubuntu`, ref2015, WSL2 Ubuntu 22.04 / Python 3.10.12 |
| 구조 렌더링 | pymol-open-source 3.2.0a0, 헤드리스(`pymol -cq`), 같은 WSL2 환경 |

셀 `[45]`·`[46]`과 `scripts/06_render.py`는 네이티브 Windows에서 돌지 않습니다.
PyRosetta는 Linux 휠만 제공하고, PyMOL은 PyPI Windows 휠이 의존 DLL을 번들하지
않아 import에서 실패합니다. 해당 단계는 WSL2에서 실행했습니다.

전 단계의 seed·버전·플랫폼은 `work/04_metrics/provenance.json`과 각 출력 CSV의 컬럼에 기록됩니다.

---

## 결과

### 표 1a — 예측기 × 그룹 요약 (scRMSD 중앙값, Å)

|  | n | scRMSD | TM-score | pLDDT(전체) | **pLDDT(재설계 영역)** | ipTM | seq_recovery |
|---|---:|---:|---:|---:|---:|---:|---:|
| S1 design_T0.1 | 14 | **0.90** | 0.97 | 90.4 | **85.0** | — | 0.52 |
| S1 design_T0.2 | 14 | **1.10** | 0.96 | 89.2 | **83.1** | — | 0.50 |
| S1 design_T0.3 | 14 | **1.36** | 0.94 | 89.4 | **82.1** | — | 0.50 |
| S1 shuffle | 5 | 3.01 | 0.83 | 78.9 | 56.3 | — | 0.06 |
| S1 random | 5 | 3.90 | 0.77 | 79.9 | 58.4 | — | 0.00 |
| S2 native | 1 | 1.69 | 0.92 | 78.9 | 60.7 | **0.12** | 1.00 |
| S2 design_T0.1 | 2 | 1.35 | 0.93 | 85.7 | 72.5 | 0.48 | 0.52 |
| S2 design_T0.2 | 2 | 1.38 | 0.93 | 85.1 | 73.3 | 0.40 | 0.52 |
| S2 design_T0.3 | 2 | 1.59 | 0.94 | 81.6 | 67.5 | 0.18 | 0.54 |
| S2 shuffle | 3 | 3.20 | 0.80 | 76.2 | 51.3 | 0.14 | 0.12 |
| S2 random | 3 | 4.81 | 0.70 | 78.6 | 53.5 | 0.14 | 0.03 |

S1 = ESMFold 단량체, S2 = AF2-multimer 복합체. 두 예측기는 입력도 모델도 다르므로 한 축에서 비교하지 않습니다.
S2 13건의 개별 값은 [`work/04_metrics/table1.md`](work/04_metrics/table1.md) 표 1b 참조.

### Fig 1 — 그룹별 scRMSD 분포

![Fig 1 — scRMSD by group, split by predictor](figures/fig1_scrmsd.png)

예측기별 2패널(왼쪽 S2 AF2-multimer 복합체, 오른쪽 S1 ESMFold 단량체).
점선은 사전 등록 임계값 2.0 Å. N이 작으므로 개별 점을 모두 표시했습니다.

### Fig 2 — scRMSD vs pLDDT

![Fig 2 — scRMSD versus pLDDT](figures/fig2_scrmsd_vs_plddt.png)

색 = 그룹, 모양 = 예측기(`o` S2 / `^` S1).
오른쪽 위에 외따로 찍힌 `random [S2]` 점이 random_s3입니다.
pLDDT가 85.4로 높은데 scRMSD는 5.50 Å으로 가장 나쁩니다.
신뢰도가 높다고 정확한 것은 아니라는 사실을, 이 한 점이 보여줍니다.

### Fig 3 — native vs 설계 중첩

![Fig 3 — native versus design superposition](figures/fig3_superposition.png)

회색 = native 사슬 A, 파랑 = 설계 `M2_T0.3_s2`(ESMFold S1 예측).
주황·연분홍 = 각각의 CDR3. 노란 막대 = Cys110.

중첩 대상으로는 중앙값에 해당하는 설계(1.21 Å)를 골랐습니다. 최저값을 쓰면
가장 잘 된 것만 보여주는 셈이 되기 때문입니다. 따라서 한계 5에서 지적한
S2 선정 편향과 같은 기준을 그림에도 적용하였습니다. 중첩은 전체 Cα로, 위에 보고한 scRMSD와 같은 정의입니다.

β-샌드위치 프레임워크는 거의 겹치고 CDR3에서 벌어집니다. 재설계 영역이
설계 자유도를 받은 자리이므로 예상되는 양상이며, 수치로는 (3)의 pLDDT 차이와 대응합니다.

### Fig 4 — 계면 접촉 잔기

![Fig 4 — interface contact residues](figures/fig4_interface_contacts.png)

연청록 표면 = 항원(사슬 B, Toxin A CROPs), 회색 cartoon = VHH(사슬 A).
주황 막대 = 접촉 잔기 16개 중 CDR3와 겹치지 않는 15개, 연분홍 = CDR3.
회색 점선 = 3.5 Å 이내 극성 접촉.

이 그림이 [`work/01_clean/README_definition_note.md`](work/01_clean/README_definition_note.md)의
정의 변경 사유를 시각적으로 보여줍니다. 항원과 맞닿는 면(주황)과 CDR3(연분홍)가
서로 다른 쪽을 향하고 있고, 16개 접촉 잔기 중 CDR3에 드는 것은 108번 하나뿐입니다.
인터페이스 정의만 썼다면 M2가 CDR3를 거의 재설계하지 않았을 것이고,
그래서 CDR3 전체를 합집합으로 더했습니다.

접촉 잔기를 IMGT CDR1·CDR2에 귀속시키지는 않았습니다. ANARCI를 쓰지 않은
탓에 numbering scheme이 검증되지 않은 상태이기 때문입니다(한계 7).

### (1) H1 지지 — 두 예측기 모두

| stage | design | shuffle | random |
|---|---|---|---|
| S1 (n=52) | **1.20 Å** | 3.01 Å | 3.90 Å |
| S2 (n=13) | **1.55 Å** | 3.20 Å | 4.81 Å |

조성을 보존한 셔플 대조군이 3 Å 이상이라는 점이 핵심입니다. 설계가 1.2 Å을 받은 것은
아미노산 조성이 그럴듯해서가 아닙니다. 위치별로 어떤 잔기를 고르느냐에서
나온 결과입니다.

### (2) 샘플링 온도의 효과는 scRMSD에서만 보인다

| temperature | scRMSD (S1) | seq_recovery | pLDDT |
|---|---|---|---|
| 0.1 | **0.90** | 0.52 | 90.4 |
| 0.2 | **1.10** | 0.50 | 89.2 |
| 0.3 | **1.36** | 0.50 | 89.4 |

recovery와 pLDDT는 평평한데 scRMSD만 단조 증가합니다. 다양성–정확도 트레이드오프가
실재하되 서열 복원율로는 보이지 않는다는 것으로, seq_recovery를 품질 지표로 쓰면 안 되는
이유를 본 데이터가 직접 보여줍니다.

![Fig 5 — CDR3 backbone ensemble by temperature](figures/fig5_cdr3_ensemble.png)

Fig 5 — 온도별 CDR3 백본 앙상블. 파랑 T0.1 / 초록 T0.2 / 빨강 T0.3 / 검정 native,
온도당 4건. 프레임워크(CDR3 제외) Cα로 중첩했습니다. 루프끼리 비교하는 것이 목적이라
Fig 3의 전체 Cα 기준과 다릅니다.

위 표의 scRMSD 증가가 백본에서 어떤 모습인지 보여줍니다. 다만 온도당 4건이므로
분포가 아니라 예시입니다. 이 그림으로 다양성 증가를 정량화하지는 않습니다.

### (3) pLDDT는 재설계 영역에서 재야 판별력이 있다

전체 평균으로는 설계 90 vs 대조군 79로 10점 차이지만, 재설계 34자리만 보면
85 vs 56으로 26점 차이입니다. 무작위 서열의 전체 pLDDT가 79점이나 되는 이유는
129잔기 중 95잔기(프레임워크)를 건드리지 않아 평균이 희석되기 때문입니다.

### (4) M1(전체 재설계)과 M2(CDR 한정)는 차이가 없다

| mode | n | scRMSD | pLDDT(재설계 영역) | seq_recovery |
|---|---:|---:|---:|---:|
| M1 (129자리) | 19 | 1.19 | 83.8 | 0.50 |
| M2 (34자리) | 23 | 1.21 | 84.9 | 0.53 |

M2는 프레임워크 95자리를 원본 그대로 물려받는데도 사슬 전체 scRMSD가 M1과 같습니다.
M1이 새로 쓴 95자리에서 오차가 유의하게 발생했다면 M1이 나빠야 하므로, 이 결과는
프레임워크 재현은 쉽고 잔여 오차가 재설계 루프에 집중됨을 시사합니다.

### (5) 이황화결합 보존 — 설계 100%, 셔플 12.5%, 무작위 0%

재설계 영역에 포함된 Cys110(나노바디 CDR1–CDR3 추가 이황화결합의 파트너)을
설계 48건이 전부 보존했습니다. → `work/04_metrics/cys_retention.csv`

### (6) ipTM 측정 실패 — 음성 결과

양성 대조군인 native 서열의 ipTM이 0.12입니다. 실제로 결정화된 복합체인데
사전 등록 임계값 0.5의 1/4 수준입니다.

양성 대조군이 통과하지 못하는 지표로는 설계를 평가할 수 없습니다. 설계 2건이 0.82·0.67로
native보다 높게 나왔지만, 같은 실행에서 무작위 치환 대조군 random_s3가 0.59를
기록했습니다(scRMSD 5.50 Å, 최악의 구조). 이 분포는 신호가 아니라 잡음입니다.

따라서 본 산출물은 ipTM을 결과로 보고하지 않고 측정 실패로 보고하며, 계면 재현에
대해서는 어떤 주장도 하지 않습니다.

### (7) 계면 점수 — Rosetta 교차검증 (탐색적)

(6)에서 ipTM이 실패해 계면에 대해 아무 주장도 하지 못한 상태였습니다.
AlphaFold-Multimer와 오차의 원천이 다른 경로로 다시 평가했습니다.
native 결정 백본에 설계 서열을 threading하고 PyRosetta로 계면을 분해했습니다.
판정 기준은 실행 전에 [`preregistration.md`](preregistration.md) §0.4에 고정했습니다.

> "독립적인 검증"이라고 쓰지 않는 이유. Rosetta는 AF2로부터 독립이지만
> ProteinMPNN으로부터 독립은 아닙니다. ProteinMPNN은 PDB 구조로 학습된
> 신경망이고 ref2015는 PDB 통계가 들어간 에너지 함수라, 둘 다 "이 백본에
> 이 서열이 어울리는가"를 비슷한 근거로 판단합니다. 아래에서 설계가 셔플보다
> 좋게 나온 데에는 이 공통 전제가 일부 기여했을 수 있습니다.

| 그룹 | dG_separated (REU) | hbonds_int | unsat_hbonds | packstat |
|---|---|---|---|---|
| **native** | **−49.55** | **13** | **5.0** | 0.73 |
| design_T0.1 | −30.52 | 8 | 14.5 | 0.69 |
| design_T0.2 | −30.42 | 8 | 14.0 | 0.70 |
| design_T0.3 | −31.23 | 8 | 12.5 | 0.70 |
| random | +37.46 | 6 | 15.5 | 0.70 |
| shuffle | +105.16 | 5 | 13.0 | 0.66 |

그룹 중앙값, `work/05_rosetta/rosetta_interface.csv` (65건 × 3반복의 중앙값).

![Fig 6 — Rosetta interface score](figures/fig6_rosetta_interface.png)

부호가 갈립니다. 설계는 계면 에너지가 음수(유리)이고 대조군은 양수입니다.

| 그룹 | dG < 0 비율 |
|---|---|
| 설계 | **100%** (42/42) |
| random | 20% (1/5) |
| shuffle | **0%** (0/5) |

설계 최대 −22.72와 대조군 최소 −22.54 사이에 겹침이 없습니다(간격 0.18 REU).
scRMSD(설계 ≤ 2.11 Å, 대조군 ≥ 2.51 Å)와 마찬가지로 완전 분리입니다.

판정 (n=42 설계 / 대조군 각 5, S1 병합 기준)

| 기준 | 값 | 결과 |
|---|---|---|
| 1. 설계 < shuffle | −31.18 vs +116.25 | **통과** |
| 2. native가 설계 상위 50% 이내 | −49.55 = **상위 0%** | **통과** |

ipTM과 달리 양성 대조군이 작동합니다. native가 설계 42건 전부보다 좋게 나왔습니다.
(6)에서 native ipTM이 0.12로 임계값의 1/4이었던 것과 대조됩니다.

#### 물리 점수와 구조 재현은 같은 것을 보는가

| Spearman ρ | dG ~ scRMSD | dG ~ pLDDT_cdr | hbonds ~ scRMSD |
|---|---|---|---|
| 전체 (대조군 포함, n=52) | **+0.50** | **−0.42** | **−0.49** |
| 설계만 (n=42) | +0.06 | +0.09 | −0.01 |

두 축은 큰 질문에는 같은 답을, 작은 질문에는 독립된 답을 줍니다.

- 설계 vs 대조군 — 둘 다 같은 방향을 가리킵니다(ρ ≈ 0.5). 부호도 맞습니다:
  dG가 높을수록(불리) scRMSD가 크고(나쁨), pLDDT_cdr은 낮습니다.
- 설계끼리 — 상관이 0입니다. 어느 설계가 더 나은지에 대해 두 지표는 서로 다른 말을 합니다.

후자는 §0.4 기준 3에서 실행 전에 예상해 둔 바입니다. "동적 범위가 좁아 판별력이 없다".
파일럿에서 설계 하나의 반복 간 산포가 13.8 REU로 설계–native 차이와 같은 자릿수였던 것과
일관됩니다. 그래서 이 지표를 개별 설계 순위화에 쓰지 않습니다.

그리고 상관이 중간 수준(0.5)에 그친 것 자체가 교차검증으로서의 가치입니다.
1.0이면 scRMSD를 다시 잰 것에 불과하고, 0이면 둘 중 하나가 신호가 없다는 뜻입니다.

#### 그래도 계면 재현을 주장하지 않습니다

(6)의 실패는 "복합체를 예측하지 못했다"는 것인데, 이 스테이지는 복합체를
예측하지 않습니다. native 결합 양식을 전제로 놓고 "이 서열이 그 기하에
들어맞는가"를 물었을 뿐입니다. 결합 양식 자체는 여전히 가정입니다.

즉 복구된 것은 *계면 예측 능력*이 아닙니다. 고정 백본에서의 서열–기하
적합성이라는 별개의 축입니다.

#### 설계가 native에 못 미치는 지점

| | native | 설계 |
|---|---|---|
| dG_separated | −49.6 | −30.7 (약 62%) |
| 계면 수소결합 | 13 | 8 |
| 매몰 미충족 극성 원자 | 5.0 | 12.5~14.5 |

설계가 native보다 낫다고 읽어서는 안 됩니다. 세 지표 모두 native가 앞섭니다.
마지막 줄은 대조군과도 구별되지 않아 따로 다룹니다. 아래 (8) 참조.

#### 보조 지표 두 개는 정보가 없습니다

- dSASA_int: 1204~1280 Å²로 전 그룹이 비슷합니다. 고정 백본 threading이라
  계면의 모양이 애초에 동일하기 때문입니다. 설계의 성질이라기보다 측정 설계가
그렇게 되어 있어서 나온 값입니다.
  이 값을 "설계가 계면을 잘 만들었다"의 근거로 쓸 수 없습니다.
- packstat: 0.66~0.73으로 분리가 약합니다.

정보가 있는 것은 dG_separated와 hbonds_int 둘뿐입니다.

#### 관측: shuffle이 random보다 나쁩니다

| | scRMSD | dG_separated |
|---|---|---|
| shuffle | 3.01 Å | **+105** |
| random | 3.90 Å | +37 |

두 지표가 음성 대조군의 순위를 반대로 매깁니다. 고정 백본에서는 조성이 보존된
shuffle이 원래의 큰 잔기를 엉뚱한 자리에 놓아 입체 충돌을 만드는 반면, 구조를 다시
접을 수 있는 scRMSD에서는 그 벌점이 완화되기 때문일 수 있습니다.
이는 아직 검증되지 않은 가설입니다. 확인하려면 반발 항(`fa_rep`)을 분해해야
하는데, 본 산출물의 범위 밖에 두었습니다.

#### 이 지표가 아닌 것

`dG_separated`는 Rosetta Energy Unit이며 kcal/mol이 아닙니다. 절대값을 해석하지
않습니다. 결합 친화도 예측이 아니라 후보 순위화의 보조 지표이며, 실험 ΔΔG와의 상관은
중간 수준으로 보고되어 있습니다.

상세는 [`work/05_rosetta/README_rosetta_note.md`](work/05_rosetta/README_rosetta_note.md).

### (8) 설계가 복원하지 못한 축 — 매몰 미충족 극성 원자

(7)의 보조 지표 중 하나가 다른 지표들과 반대 방향을 가리킵니다.
설계가 대조군과 구별되지 않는 유일한 축이라 따로 적습니다.

| 지표 | native | 설계 | shuffle | random | 설계가 대조군과 구별되는가 |
|---|---|---|---|---|---|
| dG_separated (REU) | −49.6 | −30.7 | +105.2 | +37.5 | ✅ 완전 분리 |
| 계면 수소결합 수 | 13 | 8 | 5 | 6 | ✅ 단조 |
| **매몰 미충족 극성 원자** | **5.0** | **12.5~14.5** | 13.0 | 15.5 | ❌ **겹침** |

매몰 미충족 극성 원자란 계면에 묻혔는데 수소결합 짝을 찾지 못한 극성
원자를 말합니다. 극성 원자가 용매에 노출되어 있으면 물과 수소결합을 하지만,
계면에 묻히면 상대 분자에서 짝을 찾아야 합니다. 못 찾으면 결합에 불리합니다.

native는 5.0인데 설계는 12.5~14.5로 약 3배이고, 이 값이 셔플·무작위
대조군(13.0~15.5)과 겹칩니다. 다른 두 지표에서는 설계가 대조군과 완전히
갈렸으므로, 저는 이 겹침을 우연보다 축의 성격 차이로 읽었습니다.

읽는 법. ProteinMPNN은 "이 백본에 어울리는 서열"을 찾도록 학습된 모델입니다.
백본 적합성은 scRMSD 1.20 Å과 dG −30.7 REU로 확인됩니다. 그러나 계면 건너편의
수소결합 상대까지 맞추는 것은 다른 문제이고, 이 지표는 그쪽을 재고 있습니다.
설계가 그 축에서는 무작위 서열과 구별되지 않았습니다.

이것이 본 산출물에서 기존 문헌의 재현이 아닌 유일한 관찰입니다.
다만 다음 조건 아래에서만 성립합니다.

- 단일 타깃(8YVJ), 단일 계면
- 고정 백본 threading — 설계가 백본을 바꿔 짝을 만들 가능성은 배제됨
- 통계 검정 없음. 그룹 중앙값의 차이만 봄

일반화할 수 없습니다. 타깃을 늘려 재현되는지 확인하는 것이
[다음 단계](#다음-단계) B-3의 동기 중 하나입니다.

### 사전 등록 기준 통과 집계

| 적용 | 통과 |
|---|---|
| S1 (scRMSD + pLDDT) | 40 / 52 |
| S2 (scRMSD + pLDDT) | 4 / 13 |
| **S2 (scRMSD + pLDDT + ipTM 전부)** | **설계 6건 중 1건** (M1_T0.1_s2) |
| S2 native (양성 대조군) | **불통과** (pLDDT 78.9, ipTM 0.12) |
| §0.4 기준 1 (설계 < shuffle, 탐색적) | **통과** (−31.2 vs +116.3) |
| §0.4 기준 2 (native 상위 50%, 탐색적) | **통과** (−49.6, 상위 0%) |

---

## 한계

1. 실험 검증 없음. in silico self-consistency는 발현성·결합력의 대리 지표가 아닙니다.
   델타 A의 계면 점수도 계산값이며 결합을 입증하지 않습니다.
2. ipTM 측정 실패. 양성 대조군이 임계값에 미달하여 계면 재현을 평가할 수 없었습니다.
   결과 7의 Rosetta 점수는 이 공백을 부분적으로만 메웁니다. 복합체를 예측하는
   방법이 아니고, native 결합 양식을 전제한 고정 백본 평가이기 때문입니다.
   AlphaFold-Multimer의 항체–항원 복합체 예측 성능이 일반 복합체보다 크게 낮다는 것은
   보고된 바입니다. 427개 복합체 벤치마크에서 상위 1순위 예측이 medium 이상 정확도에
   도달한 비율이 나노바디–항원 27%(n=132), 기존 항체 14%(n=295)로,
   본 결과는 실패하는 다수에 해당할 수 있습니다 (Yin & Pierce 2024).
   다만 본 실행이 그 이유로 실패했다고 단정하지는 않습니다. 축소 설정(모델 1개·recycle 3·
   단일 seed)과 항원의 β-solenoid 반복 구조도 후보이며, 이를 분리하는 실험은 하지 않았습니다.
3. n이 작습니다. S1은 조건당 14건, S2는 조건당 1–3건. 통계적 검정은 수행하지 않았습니다.
4. AF2 모델 1개 · recycle 3회 · 단일 seed로 축소 실행 — 예측 불확실성을 과소평가합니다.
5. S2 선정 편향. 13건은 ESMFold pLDDT 상위로 골랐으므로 설계 분포의 대표 표본이 아니라
   "최선의 설계들"입니다. 또한 M1 5건 / M2 1건으로 치우쳐 복합체 수준의 모드 비교는 불가능합니다.
6. 단일 타깃. 일반화를 주장할 수 없습니다.
7. 재설계 영역 정의를 원문 프로토콜에서 변경했습니다. ANARCI 대신 "인터페이스(4.5 Å) ∪ CDR3"를
   사용했으며, IMGT CDR1/CDR2 전체를 포함하지 않습니다. 사유와 영향은
   [`work/01_clean/README_definition_note.md`](work/01_clean/README_definition_note.md)에 기록했습니다.
8. US-align 미사용. Biopython Cα 중첩으로 대체했습니다. inverse folding은 길이를 보존하므로
   Cα 대응은 1:1로 정확하지만, TM-score는 최적 중첩이 아닌 RMSD-최적 중첩에서 계산한
   근사값(`tm_score_approx`)입니다.
9. 재현성의 범위. seed 42로 고정했으나 PyTorch 난수 스트림은 CUDA와 CPU가 다릅니다.
   본 결과는 Colab T4(GPU) 실행본이며, 재현성은 "동일 장치·동일 seed"로 한정됩니다.
10. ProteinMPNN 설계는 ESMFold가 더 높은 pLDDT를 주었습니다 (설계 89–90 vs native 79, S1 기준).
    이를 "설계가 원본보다 낫다"로 읽어서는 안 됩니다. pLDDT가 재는 것은 정확도가
    아닙니다. 예측기의 자신감일 뿐이고, random_s3(pLDDT 85.4, scRMSD 5.50 Å)가 그
    반례입니다.
11. 델타 A는 고정 백본 가정 위에 있습니다. 설계 서열을 native 백본에 threading했으므로,
    설계가 실제로 유발할 수 있는 백본 변화는 반영되지 않습니다. 이 가정은 scRMSD
    중앙값 1.2 Å이 간접적으로 지지할 뿐 입증하지는 않습니다.
12. dG_separated는 순위화 보조 지표입니다. 단위는 REU이고 kcal/mol이 아니며,
    `InterfaceAnalyzer`의 ddG는 실험 ΔΔG와의 상관이 중간 수준으로 보고되어 있습니다.
    개별 설계 간 순위에는 사용하지 않았습니다(반복 간 산포 13.8 REU).
13. 델타 A는 사후 추가된 탐색적 스테이지입니다. 원 프로토콜(2026-10-22)에 없었고,
    원 가설 H1의 검정에 사용하지 않았습니다. 기준은 실행 전에 §0.4로 등록했으며
    선후는 커밋 순서로 확인할 수 있습니다.

---

## 다음 단계

위 한계 중 무엇을 먼저 고칠 것인가를 비용과 효과로 정리합니다.
번호는 한계 목록의 번호입니다.

### A. 지금 자원으로 할 수 있는 것

추가 GPU 없이, 이미 구성한 WSL2 환경에서 가능합니다.

| 한계 | 할 일 | 비용 | 얻는 것 |
|---|---|---|---|
| 8 | **US-align 설치** | 수 시간 | TM-score를 근사값(`tm_score_approx`)이 아닌 정식값으로 |
| 7 | **ANARCI 설치** | 반나절 | IMGT 정통 CDR 정의. M2를 "CDR 한정 설계"로 정확히 부를 수 있게 됨 |
| — | **shuffle > random 역전 원인 규명** | CPU 1시간 | 결과 (7)의 가설을 검증. Rosetta 에너지를 항별로 분해해 `fa_rep`(반발) 기여를 확인 |
| 11 | **백본 유연 FastRelax** | CPU 5시간 | 고정 백본 가정의 영향 측정. `BB_FLEX=True`로 다시 돌려 결과가 뒤집히는지 확인 |
| 3 | **설계 수 확대** | S1이 병목 | ProteinMPNN 자체는 CPU로도 빠름. 48 → 200건은 ESMFold 시간만 늘어남 |

첫 두 개는 "도구를 못 써서 우회한 것"이고, 나머지는 "결론이 바뀔 수 있는지 확인"입니다.
후자가 더 중요합니다. 특히 백본 유연 재실행은 결과 (7) 전체의 전제를 검증합니다.

---

### B. GPU가 있다면 — 어떤 설계로 갈 것인가

규모를 키우는 것보다 먼저 할 일이 있습니다. 아래는 우선순위 순입니다.

#### B-1. ipTM 실패의 원인 분리 (최우선)

한계 2는 지금 "실패했다"까지만 말하고 원인은 셋 중 무엇인지 모릅니다.

```
(a) AF2-multimer가 항체–항원에 약하다  (Yin & Pierce 2024: 나노바디 27%)
(b) 축소 설정 (모델 1개 · recycle 3 · 단일 seed)
(c) 항원의 β-솔레노이드 반복 구조
```

GPU가 있으면 요인 분해 실험으로 갈라낼 수 있습니다. native 서열 하나만 가지고:

| 요인 | 수준 |
|---|---|
| 모델 | 1개 vs 5개 전부 |
| recycle | 3 vs 20 |
| seed | 1개 vs 5개 |

= 전체 조합 약 50회 실행.

- native ipTM이 모든 조합에서 낮으면 → (b)가 아니라 (a)/(c). 한계 2의 서술이 확정됩니다.
- 설정을 올렸을 때 ipTM이 올라가면 → 축소 실행이 원인이었고, 사전 등록 기준을 다시 적용할 수 있습니다.

저는 이것이 가장 가치가 크다고 봅니다. 데이터를 늘리는 작업이 아니라, 모르는
것을 아는 것으로 바꾸는 실험이기 때문입니다.

#### B-2. 예측기 다중화

| 현재 | 확장 |
|---|---|
| S1: ESMFold | + ESMFold (유지) |
| S2: AF2-multimer | + **Boltz-2**, **Chai-1** |

목적은 성능 향상이 아니라 실패의 귀속입니다.

- 세 예측기가 모두 native 복합체를 못 맞히면 → 타깃이 어려운 것
- AF2만 못 맞히면 → AF2 고유의 한계

AlphaFold3는 가중치에 비상업 제한이 있어 직접 실행 대신 공개 가중치 모델로 대체합니다.
이 대체 자체를 README에 근거와 함께 적는 것이 맞습니다.

#### B-3. n과 타깃 확대 — 통계를 쓸 수 있게

현재 한계 3·6이 동시에 걸려 있습니다. 조건당 14건, 타깃 1개.

```
목표:  타깃 5~10개 (나노바디–단백질 항원 복합체)
       타깃당 설계 100건 이상
       → 조건당 수백 건
```

여기까지 와야 통계적 검정이 정당해집니다. 지금은 기술통계만 보고하고 있고,
그건 n이 작아서 내린 올바른 판단이지만 영구적인 해법은 아닙니다.

타깃 선정은 이미 `data/chain_composition.csv`의 2차 선별 로직을 그대로 확장하면 됩니다.

#### B-4. 고정 백본 가정 해소 (한계 11)

```
현재:   설계 서열 → native 백본에 threading → Rosetta
확장:   설계 서열 → Boltz-2로 복합체 예측 → Rosetta
```

후자는 결합 양식을 가정하지 않습니다. 다만 B-1에서 복합체 예측이 믿을 만하다는 것이
먼저 확인되어야 합니다. 순서가 중요합니다. B-1을 건너뛰고 B-4로 가면
믿을 수 없는 구조 위에 계면 점수를 쌓는 것이 됩니다(결과 7에서 S2 복합체를 쓰지 않은 것과 같은 이유).

#### B-5. de novo 설계로 확장

현재는 inverse folding입니다. 백본이 주어지고 서열을 채웁니다.
RFdiffusion을 앞에 붙이면 백본 자체를 생성하는 파이프라인이 됩니다.

```
RFdiffusion (백본 생성)
   → ProteinMPNN (서열 설계)      ← 현재 레포는 여기만
   → ESMFold (단량체 필터)
   → Rosetta (계면 점수)           ← 델타 A 스테이지 재사용
   → Boltz-2 / Chai-1 (복합체 검증)
```

이 레포의 스테이지가 그대로 중간 모듈이 됩니다. `work/05_rosetta/`의 계면 점수와
`scripts/06_render.py`의 렌더링은 파이프라인이 길어져도 재사용됩니다.

#### B-6. 실험 가능성 필터 추가

in silico 점수가 좋아도 발현·정제가 안 되면 실험으로 넘길 수 없습니다.

- 응집 경향(aggregation propensity), 등전점, 소수성 패치
- 발현 숙주에 맞는 코돈·잔기 제약

어느 단계에서 거르느냐가 설계 판단입니다. 맨 뒤에 붙이면 좋은 설계를 늦게 버리고,
앞에 붙이면 탐색 공간을 미리 좁힙니다. ESMFold 필터 직후가 합리적인 위치로 보입니다.

---

### C. 규모가 커져도 바꾸지 않을 것

위를 전부 해도 아래는 유지합니다. 이 레포의 결론이 믿을 만하다면 그 이유는
규모보다 설계에 있기 때문입니다.

1. 판정 기준을 실행 전에 등록한다. 델타 A에서는 커밋 순서로 선후를 증명했습니다
   (`e5ce22f` 기준 등록 → 결과 커밋).
2. 대조군을 같은 파이프라인에 태운다. 셔플·무작위가 없었다면 scRMSD 1.2 Å도
   계면 −31 REU도 해석할 수 없습니다.
3. 양성 대조군이 실패하면 측정 실패로 보고한다. ipTM이 그랬고, 델타 A는 그 조항을
   미리 적어 둔 뒤 통과했습니다. 통과했기 때문에 의미가 있는 것은 아닙니다.
   실패했더라도 그렇게 보고했을 것이기 때문에 의미가 있습니다.
4. 무엇을 하지 않았는지 먼저 적는다.

---

## 재현 방법

```bash
git clone <this-repo> && cd vhh-mpnn-selfconsistency
# Colab에서 protein_design_sc.ipynb 를 열고 런타임을 T4 GPU로 설정
# 셀 01 → 50 을 순서대로 실행 (10·11은 타깃 선정 기록용, 생략 가능)
```

노트북은 자신이 있는 폴더를 저장소 루트로 인식하므로(`LOCAL_ROOT = cwd`),
clone 위치에 관계없이 그대로 동작합니다. 실행 결과를 담은 노트북은
[`protein_design_sc_executed.ipynb`](protein_design_sc_executed.ipynb)에 출력과 함께 보존돼 있습니다.

노트북은 Colab과 로컬 Windows 양쪽에서 동작합니다. 로컬 환경 구축과 단계별 실행 환경은
[`COLAB_RUNBOOK.md`](COLAB_RUNBOOK.md) 참조. GPU 없이 로컬에서 돌릴 경우 ProteinMPNN(CPU, 8서열/27초)까지만
가능하며, 재예측 단계는 GPU가 필요합니다.

---

## 산출물

```
protein_design_sc.ipynb         실행용 노트북 (Colab / 로컬 공용)
protein_design_sc_executed.ipynb  출력이 포함된 실행본 — 본 결과의 원본 기록
protein_design_sc_protocol.md   사전 작성한 기획·프로토콜 문서 (2026-10-22)
data/
  target_candidates.csv       RCSB 조회 결과
  chain_composition.csv       사슬 구성 2차 선별 — 타깃 선정 근거
work/01_clean/
  8YVJ_clean.pdb              전처리 구조
  designable.json             재설계 영역 정의 + 변경 사유
  README_definition_note.md   정의 변경 상세
work/02_mpnn/
  all_sequences.csv           65 서열 (설계 48 + 대조군 16 + native 1)
  all_monomer.fasta / all_complex.fasta / selected.fasta
work/03_refold/
  S1_esmfold/                 ESMFold 구조 65건
  S2_colabfold/               AF2-multimer 구조 13건 + MSA/PAE/점수
work/04_metrics/
  screening.csv               S1 pLDDT
  metrics.csv                 전 지표 + provenance
  cys_retention.csv           이황화결합 보존율
  table1.md                   표 1a / 1b
  provenance.json             seed·버전·플랫폼
work/05_rosetta/                (탐색적 — 사전 등록 §0.4)
  rosetta_interface.csv       설계별 계면 점수 중앙값 + provenance
  rosetta_interface_raw.csv   반복별 원시값 (산포 확인용)
  README_rosetta_note.md      가정·한계·라이선스 고지
figures/
  fig1_scrmsd.png  fig2_scrmsd_vs_plddt.png
  fig3_superposition.png  fig4_interface_contacts.png  fig5_cdr3_ensemble.png
  fig6_rosetta_interface.png
  scripts/*.pml               그림을 만든 PyMOL 스크립트 (생성물이자 검증 대상)
  README.md                   그림별 선택 근거
scripts/
  06_render.py                .pml 생성 + 헤드리스 렌더링
```

상세 분석과 해석은 [`RESULTS_ANALYSIS.md`](RESULTS_ANALYSIS.md) 참조.

---

## 참고문헌

한계 2번의 근거입니다. 상세 논의와 추가 문헌은 [`RESULTS_ANALYSIS.md`](RESULTS_ANALYSIS.md) §11.

- Yin R, Pierce BG (2024) Evaluation of AlphaFold antibody–antigen modeling with implications
  for improving predictive accuracy. *Protein Science* 33(1):e4865.
  [doi:10.1002/pro.4865](https://doi.org/10.1002/pro.4865)
- Yin R, Feng BY, Varshney A, Pierce BG (2022) Benchmarking AlphaFold for protein complex
  modeling reveals accuracy determinants. *Protein Science* 31(8):e4379.
  [doi:10.1002/pro.4379](https://doi.org/10.1002/pro.4379)

---

## 라이선스 고지

ProteinMPNN은 MIT License (Copyright (c) 2022 Justas Dauparas) 하에 사용했습니다.
