# VHH inverse-folding self-consistency benchmark

ProteinMPNN으로 재설계한 나노바디 서열이, 독립적인 구조 예측기로 되돌렸을 때
원본 backbone을 재현하는가(self-consistency)를 음성·양성 대조군과 함께 측정한 벤치마크 실습.

타깃: **PDB 8YVJ** — H5.2 나노바디(VHH, 129 aa) / *C. difficile* toxin A CROPs 단편(118 aa), X-ray 1.65 Å

---

## 무엇을 했는가

1. RCSB Search API로 나노바디 복합체 후보를 조회하고, **사슬 구성으로 2차 선별**해 타깃을 확정했습니다.
   해상도 ≤ 2.5 Å · 총 ≤ 300잔기 조건을 만족한 3건 중, 항원이 펩타이드가 아닌 단백질인 것은 8YVJ뿐이었습니다
   (나머지 2건은 14 aa Tau PHF6 펩타이드 복합체). → `data/chain_composition.csv`
2. PDBFixer로 전처리하고, 재설계 영역을 **항원과 4.5 Å 이내 접촉하는 잔기 ∪ CDR3** (34/129자리)로 정의했습니다.
   → `work/01_clean/designable.json`, 정의 변경 사유는 [`work/01_clean/README_definition_note.md`](work/01_clean/README_definition_note.md)
3. ProteinMPNN으로 **2모드 × 3온도 × 8서열 = 48 설계 서열**을 생성했습니다.
   - M1: VHH 전체 129자리 재설계 / M2: 위 34자리만 재설계(프레임워크 고정)
   - temperature 0.1 / 0.2 / 0.3
4. **대조군 17건**: 조성 보존 셔플 8, 완전 무작위 치환 8 (모두 재설계 영역 한정), native 1(양성 대조군)
5. 2단계로 재예측했습니다.
   - **S1**: ESMFold로 VHH 단량체 65건 전부
   - **S2**: ColabFold(AF2-multimer v3)로 복합체 13건 (S1 pLDDT 상위 설계 6 + 대조군 6 + native 1)
6. scRMSD · TM-score · pLDDT(전체/재설계 영역) · ipTM · seq_recovery를 산출했습니다.

---

## 무엇을 하지 않았는가

다음은 본 산출물로 **입증되지 않습니다.**

- 설계 서열이 실제로 발현·폴딩된다는 것 (실험 검증 없음)
- 설계 서열이 항원에 결합한다는 것 (결합 친화도 미측정)
- **복합체 계면을 재현한다는 것** — ipTM 측정에 실패했습니다 (아래 결과 참조)
- 통계적 유의성 (N이 작아 기술통계 수준의 비교만 수행)

**확실히 주장할 수 있는 사실:** 주어진 backbone에 대해 inverse folding으로 생성한 서열이,
독립적인 구조예측기로 되돌렸을 때 원본 backbone을 재현하며, 그 재현도가 음성 대조군보다
명확히 높다.

---

## 사전 등록된 기준

판정 기준은 [`protein_design_sc_protocol.md`](protein_design_sc_protocol.md)(2026-10-22) §0.3에 코드 실행 전부터 적혀 있었고,
결과를 본 뒤 완화하거나 변경하지 않았습니다. 다만 이를 담은 [`preregistration.md`](preregistration.md) 파일은
**분석을 마친 뒤에 작성**되었으므로, 선후를 증명할 커밋 기록은 없습니다.
그 사실과 "데이터를 본 뒤에 내린 결정" 목록을 [`preregistration.md`](preregistration.md)에 전부 적었습니다.

```
성공 판정 (설계 1건 단위)
  scRMSD (재예측 구조 vs 원본 backbone, Cα) < 2.0 Å
  재예측 pLDDT 평균 > 80
  (복합체 수행 시) ipTM > 0.5

H1: 설계 서열의 scRMSD 중앙값 < 셔플 대조군의 scRMSD 중앙값
H0: 차이 없음
```

**기준은 사후에 변경하지 않았습니다.** 아래 통과 집계는 위 기준을 그대로 적용한 결과입니다.

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

셀 `[45]`·`[46]`과 `scripts/06_render.py`는 **네이티브 Windows에서 돌지 않습니다.**
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

S1 = ESMFold 단량체, S2 = AF2-multimer 복합체. **두 예측기는 입력도 모델도 다르므로 한 축에서 비교하지 않습니다.**
S2 13건의 개별 값은 [`work/04_metrics/table1.md`](work/04_metrics/table1.md) 표 1b 참조.

### Fig 1 — 그룹별 scRMSD 분포

![Fig 1 — scRMSD by group, split by predictor](figures/fig1_scrmsd.png)

예측기별 2패널(왼쪽 S2 AF2-multimer 복합체, 오른쪽 S1 ESMFold 단량체).
점선은 사전 등록 임계값 2.0 Å. N이 작으므로 개별 점을 모두 표시했습니다.

### Fig 2 — scRMSD vs pLDDT

![Fig 2 — scRMSD versus pLDDT](figures/fig2_scrmsd_vs_plddt.png)

색 = 그룹, 모양 = 예측기(`o` S2 / `^` S1).
오른쪽 위에 외따로 찍힌 `random [S2]` 점이 random_s3입니다 —
pLDDT 85.4로 높은데 scRMSD는 5.50 Å으로 최악입니다. **신뢰도 지표가
정확도의 대리 지표가 아니라는 것**을 한 점으로 보여줍니다.

### Fig 3 — native vs 설계 중첩

![Fig 3 — native versus design superposition](figures/fig3_superposition.png)

회색 = native 사슬 A, 파랑 = 설계 `M2_T0.3_s2`(ESMFold S1 예측).
주황·연분홍 = 각각의 CDR3. 노란 막대 = Cys110.

**대상은 최저 scRMSD가 아니라 중앙값(1.21 Å)입니다.** 최저값을 쓰면
"가장 잘된 것"만 보여주는 선택 편향이 되므로, 한계 5에서 지적한 S2 선정 편향과
같은 기준을 그림에도 적용했습니다. 중첩은 전체 Cα로, 위에 보고한 scRMSD와 같은 정의입니다.

β-샌드위치 프레임워크는 거의 겹치고 CDR3에서 벌어집니다. 재설계 영역이
설계 자유도를 받은 자리이므로 예상되는 양상이며, 수치로는 (3)의 pLDDT 차이와 대응합니다.

### Fig 4 — 계면 접촉 잔기

![Fig 4 — interface contact residues](figures/fig4_interface_contacts.png)

연청록 표면 = 항원(사슬 B, Toxin A CROPs), 회색 cartoon = VHH(사슬 A).
주황 막대 = 접촉 잔기 16개 중 CDR3와 겹치지 않는 15개, 연분홍 = CDR3.
회색 점선 = 3.5 Å 이내 극성 접촉.

이 그림이 [`work/01_clean/README_definition_note.md`](work/01_clean/README_definition_note.md)의
정의 변경 사유를 시각적으로 보여줍니다. **항원과 맞닿는 면(주황)과 CDR3(연분홍)가
서로 다른 쪽을 향하고 있고, 16개 접촉 잔기 중 CDR3에 드는 것은 108번 하나뿐입니다.**
인터페이스 정의만 썼다면 M2가 CDR3를 거의 재설계하지 않았을 것이고,
그래서 CDR3 전체를 합집합으로 더했습니다.

접촉 잔기를 IMGT CDR1·CDR2로 귀속시키지는 않습니다. ANARCI를 쓰지 않아
numbering scheme을 검증하지 않았기 때문입니다(한계 7).

### (1) H1 지지 — 두 예측기 모두

| stage | design | shuffle | random |
|---|---|---|---|
| S1 (n=52) | **1.20 Å** | 3.01 Å | 3.90 Å |
| S2 (n=13) | **1.55 Å** | 3.20 Å | 4.81 Å |

조성을 보존한 셔플 대조군이 3 Å 이상이라는 점이 핵심입니다. 설계의 1.2 Å은
"그럴듯한 아미노산 조성"이 아니라 **위치별 선택**에서 나옵니다.

### (2) 샘플링 온도의 효과는 scRMSD에서만 보인다

| temperature | scRMSD (S1) | seq_recovery | pLDDT |
|---|---|---|---|
| 0.1 | **0.90** | 0.52 | 90.4 |
| 0.2 | **1.10** | 0.50 | 89.2 |
| 0.3 | **1.36** | 0.50 | 89.4 |

recovery와 pLDDT는 평평한데 scRMSD만 단조 증가합니다. 다양성–정확도 트레이드오프가
실재하되 **서열 복원율로는 보이지 않는다**는 것으로, seq_recovery를 품질 지표로 쓰면 안 되는
이유를 본 데이터가 직접 보여줍니다.

![Fig 5 — CDR3 backbone ensemble by temperature](figures/fig5_cdr3_ensemble.png)

**Fig 5 — 온도별 CDR3 백본 앙상블.** 파랑 T0.1 / 초록 T0.2 / 빨강 T0.3 / 검정 native,
온도당 4건. 프레임워크(CDR3 제외) Cα로 중첩했습니다 — 루프끼리 비교하는 것이 목적이라
Fig 3의 전체 Cα 기준과 다릅니다.

위 표의 scRMSD 증가가 백본에서 어떤 모습인지 보여줍니다. 다만 **온도당 4건이므로
분포가 아니라 예시입니다.** 이 그림으로 다양성 증가를 정량화하지는 않습니다.

### (3) pLDDT는 재설계 영역에서 재야 판별력이 있다

전체 평균으로는 설계 90 vs 대조군 79로 **10점** 차이지만, 재설계 34자리만 보면
85 vs 56으로 **26점** 차이입니다. 무작위 서열의 전체 pLDDT가 79점이나 되는 이유는
129잔기 중 95잔기(프레임워크)를 건드리지 않아 평균이 희석되기 때문입니다.

### (4) M1(전체 재설계)과 M2(CDR 한정)는 차이가 없다

| mode | n | scRMSD | pLDDT(재설계 영역) | seq_recovery |
|---|---:|---:|---:|---:|
| M1 (129자리) | 19 | 1.19 | 83.8 | 0.50 |
| M2 (34자리) | 23 | 1.21 | 84.9 | 0.53 |

M2는 프레임워크 95자리를 원본 그대로 물려받는데도 사슬 전체 scRMSD가 M1과 같습니다.
M1이 새로 쓴 95자리에서 오차가 유의하게 발생했다면 M1이 나빠야 하므로, 이 결과는
**프레임워크 재현은 쉽고 잔여 오차가 재설계 루프에 집중**됨을 시사합니다.

### (5) 이황화결합 보존 — 설계 100%, 셔플 12.5%, 무작위 0%

재설계 영역에 포함된 Cys110(나노바디 CDR1–CDR3 추가 이황화결합의 파트너)을
설계 48건이 **전부** 보존했습니다. → `work/04_metrics/cys_retention.csv`

### (6) ipTM 측정 실패 — 음성 결과

**양성 대조군인 native 서열의 ipTM이 0.12**입니다. 실제로 결정화된 복합체인데
사전 등록 임계값 0.5의 1/4 수준입니다.

양성 대조군이 통과하지 못하는 지표로는 설계를 평가할 수 없습니다. 설계 2건이 0.82·0.67로
native보다 높게 나왔지만, 같은 실행에서 **무작위 치환 대조군 random_s3가 0.59**를
기록했습니다(scRMSD 5.50 Å, 최악의 구조). 이 분포는 신호가 아니라 잡음입니다.

따라서 **본 산출물은 ipTM을 결과로 보고하지 않고 측정 실패로 보고하며, 계면 재현에
대해서는 아무것도 주장하지 않습니다.**

### (7) 계면 점수 — Rosetta 교차검증 ⟨A 대기⟩

<!-- ────────────────────────────────────────────────────────────────
작성 지침 (채우고 나면 이 주석 블록을 지울 것)

셀 [46] 출력의 두 기준에 따라 서술이 갈립니다.

[기준 2 통과 = native가 설계 분포 상위 50% 이내]
  → "ipTM과 달리 양성 대조군이 작동하는 지표를 확보했다"로 서술.
     '무엇을 하지 않았는가'의 계면 항목과 한계 2번을 함께 수정할 것.

[기준 2 불통과]
  → ipTM과 똑같이 **측정 실패로 보고**. 순위화에 쓰지 않는다.
     사전 등록 §0.4에 그렇게 하기로 미리 적어 두었으므로 그대로 따른다.

[기준 1: 설계 중앙값 < shuffle 중앙값]
  → 통과면 "조성이 같은 셔플보다 계면이 낫다"로, (1)의 scRMSD 결과와 같은 구조.
  → 불통과면 그대로 적고, 동적 범위가 좁다는 §0.4의 사전 판단과 연결.

어느 쪽이든 아래 세 문장은 유지할 것:
  - dG_separated는 REU이며 kcal/mol이 아니다. 절대값을 해석하지 않는다.
  - 결합 친화도 예측이 아니라 순위화 보조 지표다.
  - 개별 설계 순위에는 쓰지 않는다 (파일럿에서 반복 간 산포 13.8 REU 관측).
──────────────────────────────────────────────────────────────── -->

ipTM이 실패해 계면에 대해 아무 주장도 하지 못한 상태였습니다((6) 참조).
AlphaFold-Multimer에 의존하지 않는 경로로 계면을 다시 평가했습니다.
native 결정 백본에 설계 서열을 threading하고 PyRosetta로 계면을 분해했습니다.

| 그룹 | dG_separated (REU) | dSASA_int (Å²) | hbonds_int |
|---|---|---|---|
| native | ⟨A 대기⟩ | ⟨A 대기⟩ | ⟨A 대기⟩ |
| design_T0.1 | ⟨A 대기⟩ | ⟨A 대기⟩ | ⟨A 대기⟩ |
| design_T0.2 | ⟨A 대기⟩ | ⟨A 대기⟩ | ⟨A 대기⟩ |
| design_T0.3 | ⟨A 대기⟩ | ⟨A 대기⟩ | ⟨A 대기⟩ |
| shuffle | ⟨A 대기⟩ | ⟨A 대기⟩ | ⟨A 대기⟩ |
| random | ⟨A 대기⟩ | ⟨A 대기⟩ | ⟨A 대기⟩ |

![Fig 6 — Rosetta interface score](figures/fig6_rosetta_interface.png)

**판정** ⟨A 대기⟩

**물리 점수와 예측기 신뢰도의 상관** ⟨A 대기⟩
(Spearman ρ: dG vs scRMSD = ⟨A 대기⟩, dG vs pLDDT_cdr = ⟨A 대기⟩)

상관이 약하게 나오더라도 그것이 음성 결과는 아닙니다. 두 축이 독립된 정보를
담는다는 뜻이고, 교차검증 지표로서는 오히려 그쪽이 가치가 큽니다.

**이 지표가 아닌 것.** `dG_separated`는 Rosetta Energy Unit이며 kcal/mol이
아닙니다. 절대값을 해석하지 않습니다. 결합 친화도 예측이 아니라 후보 순위화의
보조 지표이며, 실험 ΔΔG와의 상관은 중간 수준으로 보고되어 있습니다.
**개별 설계 간 순위에는 사용하지 않습니다** — 파일럿에서 설계 하나의 반복 간
산포가 13.8 REU로, 설계와 native의 차이와 같은 자릿수였습니다.
이는 사전 등록 §0.4에 실행 전부터 적어 둔 판단과 일치합니다.

상세는 [`work/05_rosetta/README_rosetta_note.md`](work/05_rosetta/README_rosetta_note.md).

### 사전 등록 기준 통과 집계

| 적용 | 통과 |
|---|---|
| S1 (scRMSD + pLDDT) | 40 / 52 |
| S2 (scRMSD + pLDDT) | 4 / 13 |
| **S2 (scRMSD + pLDDT + ipTM 전부)** | **설계 6건 중 1건** (M1_T0.1_s2) |
| S2 native (양성 대조군) | **불통과** (pLDDT 78.9, ipTM 0.12) |
| §0.4 기준 1 (설계 < shuffle, 탐색적) | ⟨A 대기⟩ |
| §0.4 기준 2 (native 상위 50%, 탐색적) | ⟨A 대기⟩ |

---

## 한계

1. **실험 검증 없음.** in silico self-consistency는 발현성·결합력의 대리 지표가 아닙니다.
   델타 A의 계면 점수도 계산값이며 결합을 입증하지 않습니다.
2. **ipTM 측정 실패.** 양성 대조군이 임계값에 미달하여 계면 재현을 평가할 수 없었습니다.
   AlphaFold-Multimer의 항체–항원 복합체 예측 성능이 일반 복합체보다 크게 낮다는 것은
   보고된 바입니다 — 427개 복합체 벤치마크에서 상위 1순위 예측이 medium 이상 정확도에
   도달한 비율이 **나노바디–항원 27%(n=132)**, 기존 항체 14%(n=295)로,
   본 결과는 실패하는 다수에 해당할 수 있습니다 (Yin & Pierce 2024).
   다만 본 실행이 그 이유로 실패했다고 단정하지는 않습니다. 축소 설정(모델 1개·recycle 3·
   단일 seed)과 항원의 β-solenoid 반복 구조도 후보이며, 이를 분리하는 실험은 하지 않았습니다.
3. **n이 작습니다.** S1은 조건당 14건, S2는 조건당 1–3건. 통계적 검정은 수행하지 않았습니다.
4. **AF2 모델 1개 · recycle 3회 · 단일 seed**로 축소 실행 — 예측 불확실성을 과소평가합니다.
5. **S2 선정 편향.** 13건은 ESMFold pLDDT 상위로 골랐으므로 설계 분포의 대표 표본이 아니라
   "최선의 설계들"입니다. 또한 M1 5건 / M2 1건으로 치우쳐 복합체 수준의 모드 비교는 불가능합니다.
6. **단일 타깃.** 일반화를 주장할 수 없습니다.
7. **재설계 영역 정의를 원문 프로토콜에서 변경했습니다.** ANARCI 대신 "인터페이스(4.5 Å) ∪ CDR3"를
   사용했으며, IMGT CDR1/CDR2 전체를 포함하지 않습니다. 사유와 영향은
   [`work/01_clean/README_definition_note.md`](work/01_clean/README_definition_note.md)에 기록했습니다.
8. **US-align 미사용.** Biopython Cα 중첩으로 대체했습니다. inverse folding은 길이를 보존하므로
   Cα 대응은 1:1로 정확하지만, TM-score는 최적 중첩이 아닌 RMSD-최적 중첩에서 계산한
   **근사값**(`tm_score_approx`)입니다.
9. **재현성의 범위.** seed 42로 고정했으나 PyTorch 난수 스트림은 CUDA와 CPU가 다릅니다.
   본 결과는 Colab T4(GPU) 실행본이며, 재현성은 "동일 장치·동일 seed"로 한정됩니다.
10. **ProteinMPNN 설계는 ESMFold가 더 높은 pLDDT를 주었습니다** (설계 89–90 vs native 79, S1 기준).
    이를 "설계가 원본보다 낫다"로 읽어서는 안 됩니다. pLDDT는 정확도가 아니라 예측기의
    자신감이며, random_s3(pLDDT 85.4, scRMSD 5.50 Å)가 그 반례입니다.

---

## 재현 방법

```bash
git clone <this-repo> && cd vhh-mpnn-selfconsistency
# Colab에서 protein_design_sc.ipynb 를 열고 런타임을 T4 GPU로 설정
# 셀 01 → 50 을 순서대로 실행 (10·11은 타깃 선정 기록용, 생략 가능)
```

노트북은 **자신이 있는 폴더를 저장소 루트로 인식**하므로(`LOCAL_ROOT = cwd`),
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
