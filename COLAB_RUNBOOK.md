# Colab 실행 런북 — S1(ESMFold) · S2(ColabFold)

`protein_design_sc.ipynb` 32셀 기준. 로컬에서 셀 22까지 끝낸 상태에서 이어서 하는 절차입니다.

---

## 0. 왜 Colab인가 (실측 근거)

| 단계 | 로컬 (Windows, RAM 15.7 GB) | 조치 |
|---|---|---|
| 셀 20–22 ProteinMPNN | 조건당 27초 × 6 = 약 3분 (CPU) | **로컬 완료** |
| 셀 30 ESMFold | **불가** — 모델 fp32 8.4 GB, 가용 RAM 7 GB → segfault(exit 139) | Colab GPU |
| 셀 31 ColabFold | **불가** — `colabfold_batch`는 Linux 전용 | Colab GPU |

셀 30은 로딩 전에 가용 RAM을 점검해, 조용히 죽는 대신 원인과 선택지를 출력하고 멈춥니다.

---

## 1. 올릴 파일

**`protein_design_sc.ipynb` 하나입니다.**

`vhh_upload.zip`(144 KB)은 **올리지 않아도 됩니다.** 셀 22까지의 산출물은 Colab에서 그대로
재생성되며(RCSB 다운로드 + ProteinMPNN 클론, seed 42 고정), 샌드박스와 프로젝트 폴더에서
각각 돌렸을 때 결과가 완전히 일치하는 것을 확인했습니다.
재계산 3분을 아끼고 싶으면 zip도 함께 올리십시오 (§5 참조).

`third_party/`(180 MB)는 **절대 올리지 마십시오.** 셀 20이 다시 클론합니다.

### 올리는 법

1. <https://drive.google.com> → `4. ML DL archive` → `3. Exercises` → `ProteinMPNN`
2. `protein_design_sc.ipynb`를 드래그
3. 우클릭 → 연결 앱 → Google Colaboratory
4. **런타임 → 런타임 유형 변경 → T4 GPU** ← 빠뜨리면 셀 30이 RAM 가드에 걸립니다

---

## 2. 실행 순서

| 셀 | 하는 일 | 비고 |
|---|---|---|
| 01 → 04 | 마운트·경로·의존성·provenance | `IN_COLAB=True` 자동 전환 |
| 05 | Drive → `/content/vhh_work` 복원 | zip을 안 올렸으면 "이어받을 작업물이 없습니다" 출력 후 통과 (정상) |
| 10, 11 | 타깃 조회 | **건너뛰기.** 8YVJ로 확정됨. 100회 넘는 API 호출만 낭비 |
| 12 → 15 | 전처리 · designable · README 안내문 | |
| 20 → 22 | ProteinMPNN | 20이 레포 클론, 21이 설계 실행, 22가 65서열 조립 |
| **30** | ESMFold 65건 | 목적지. 첫 줄 `device = cuda` 확인 |
| 31 | ColabFold S2 13건 | 30 완료 후 |
| 40 → 41 | 지표 · 그림 2종 | 로컬에서 해도 됨 |

**20–22를 건너뛰고 30으로 갈 수 없습니다.** 셀 30이 쓰는 `ALL` · `NATIVE` · `DESIGNABLE` ·
`cys_in_design`이 전부 셀 22에서 만들어집니다.

---

## 3. 셀 30에서 확인할 것

### (1) `device = cuda` 인지

`cpu`로 뜨면 GPU 미할당입니다. 런타임 유형을 다시 확인하십시오.

### (2) `S1_LIMIT`을 먼저 3으로

```python
S1_LIMIT = 3      # 동작·소요시간 측정
```

서열당 실제 소요시간과 전체 추정치가 출력됩니다. 괜찮으면 `None`으로 바꿔 전체를 돌리십시오.

부분 실행 결과는 `screening.csv`가 아니라 **`screening_partial.csv`** 로 저장됩니다.
셀 31이 `screening.csv`를 읽어 S2 상위 6건을 고르는 구조라, 3건짜리 부분 결과를 전체로
오인하면 **에러 없이 틀린 13건을 선정**하게 되기 때문입니다.

### (3) 중간 저장

세션이 끊기면 셀 01 → 05를 다시 돌리면 이어집니다 (셀 30은 기존 `.pdb`를 보고 건너뜀).
단 끊기기 전에 `sync_to_drive("work/03_refold")`가 실행됐어야 합니다 — 셀 30 끝에 있지만
중간에 끊기면 그 세션 분량은 날아갑니다. **`S1_LIMIT`을 20씩 끊어 돌리는 편이 안전합니다.**

Colab 무료 티어는 유휴 90분 · 최대 12시간에 끊깁니다.

---

## 4. 알려진 함정

| 증상 | 원인 | 대응 |
|---|---|---|
| ESMFold가 메시지 없이 죽음 (exit 139) | RAM 부족 | 셀 30 RAM 가드가 미리 차단. GPU 런타임 확인 |
| `hf-xet` Rust 패닉, OSError 1450 | HuggingFace Rust 다운로더, Windows 한정 | 셀 01의 `HF_HUB_DISABLE_XET=1`이 우회 |
| `protein_mpnn_run.p/vanilla_model_weights/...` 없음 | ProteinMPNN이 `rfind("/")`로 경로 계산, Windows 한정 | 셀 21의 `--path_to_model_weights`가 우회 |
| 출력 `.fa` 파일명이 전체 경로 | `parse_multiple_chains.py`도 `/` 기준 basename, Windows 한정 | 셀 20이 `parsed.jsonl`의 `name` 정규화 |
| AF2-multimer OOM (셀 31) | T4 VRAM 부족 | `--num-recycle 1` → 그래도 안 되면 단량체 폴딩 (ipTM 포기) |

---

## 5. zip을 올려 재계산을 건너뛰려면

`vhh_upload.zip`을 `ProteinMPNN` 폴더에 함께 올린 뒤, **셀 02와 셀 05 사이에** 새 셀을 만들어
한 번만 실행하십시오.

```python
!cd "/content/drive/MyDrive/4. ML DL archive/3. Exercises/ProteinMPNN" && unzip -o vhh_upload.zip
```

그 다음 셀 05가 Drive에서 복원하고, 셀 12는 기존 `8YVJ_clean.pdb`를 재사용하며
(Colab에서 pdbfixer 설치가 실패해 다른 backbone이 만들어지는 것을 막습니다),
셀 21은 기존 FASTA를 보고 건너뜁니다.

---

## 6. 끝난 뒤 남은 숙제 — 코드가 해주지 않습니다

1. **`preregistration.md` 커밋** (§0.3). 아직이면 지금 커밋하고 시점을 정직하게 적을 것
2. **README "무엇을 하지 않았는가"** — §0.2를 그대로
3. **재설계 영역 정의 변경** — `work/01_clean/README_definition_note.md`를 README 방법 섹션에 그대로 옮길 것 (인터페이스 ∪ CDR3, ANARCI 미사용)
4. **한계 항목**
   - 실험 검증 없음 / 조건당 n=8 / 단일 타깃
   - AF2 모델 1개 · recycle 3회로 축소
   - US-align 대신 Biopython Cα 중첩을 썼다면 `tm_score_approx` 표기 유지
   - **온도 0.1–0.3 구간에서 recovery 곡선이 거의 평평했음** (M2 0.529 → 0.515 → 0.496, n=8). §10.2 질문에 이대로 답하는 편이 강함
5. **§9.2 한계 5번 (AbMPNN 문헌 주장)** — 1차 문헌으로 교차 확인하고 인용을 달거나, 어려우면 **문장 삭제**

### 보고할 만한 자체 관찰 (인용 불필요, 본인 데이터)

- **M2 고정 메커니즘 정확성:** designable 34자리 외 변경 0건
- **이황화결합 보존:** 설계 48건 전부가 Cys110 보존 (1.000) vs shuffle 0.125, random 0.000.
  나노바디 CDR1–CDR3 추가 이황화결합의 파트너이며, 조성 보존 대조군 대비 명확한 차이
- **재현성:** seed 42로 서로 다른 경로에서 두 번 실행해 동일 결과
- **A2 가정 실측:** ProteinMPNN CPU 추론 8서열/27초 (247잔기 복합체)
