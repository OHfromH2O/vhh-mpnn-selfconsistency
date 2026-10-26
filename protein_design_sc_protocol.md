# 최소 설계 산출물 생산 프로토콜 — ProteinMPNN 기반 in silico self-consistency 검증

작성일: 2026-10-22 / 목표 마감: 2026-10-28 (아론티어 Computational Protein Design Scientist 지원)

---

## 0. 이 문서의 전제와 범위 제한

### 0.1 명시적 전제

| 번호 | 전제 | 검증 필요 여부 |
|---|---|---|
| A1 | Colab 무료 티어에서 T4 GPU(16GB VRAM) 할당을 받을 수 있다 | 세션마다 변동. 할당 실패 시 §7 대체 경로 |
| A2 | ProteinMPNN은 CPU만으로도 수 초~수십 초 내 추론이 가능하다 | 모델 파라미터 수가 작아 통상 성립. 실측으로 확인할 것 |
| A3 | ColabFold(AF2-multimer)는 T4에서 총 잔기 수 300 이하 복합체를 처리할 수 있다 | 실측 필요. 초과 시 단량체로 축소 |
| A4 | Boltz-2는 T4 무료 티어에서 복합체 추론 시 VRAM 부족 가능성이 있다 | **미검증.** 선택 단계로 분리(§6.3), 실패해도 산출물이 성립하도록 설계 |

### 0.2 이 산출물이 **주장하지 않는 것**

다음은 본 프로토콜로 입증되지 않습니다. README와 면접 양쪽에서 동일하게 명시해야 합니다.

- 설계 서열이 실제로 발현·폴딩된다는 것 (실험 검증 없음)
- 설계 서열이 항원에 결합한다는 것 (결합 친화도 미측정)
- 통계적 유의성 (N이 작아 분포 비교는 기술통계 수준)

**본 산출물이 주장하는 것은 단 하나입니다:** "주어진 backbone에 대해 inverse folding으로 생성한 서열이, 독립적인 구조예측기로 되돌렸을 때 원본 backbone을 재현하는가(self-consistency)" — 그리고 그 재현도가 음성 대조군 대비 유의하게 높은가.

### 0.3 사전 등록(pre-registration) — 분석 전에 확정할 것

결과를 본 뒤 기준을 바꾸면 산출물의 신뢰도가 사라집니다. 코드 작성 전에 다음을 고정하고 `preregistration.md`로 커밋하십시오.

```
성공 판정 기준 (설계 1건 단위)
  - scRMSD (설계 서열 재예측 구조 vs 원본 backbone, Cα) < 2.0 Å
  - 재예측 pLDDT 평균 > 80
  - (복합체 수행 시) ipTM > 0.5

가설
  H1: ProteinMPNN 설계 서열의 scRMSD 중앙값 < 셔플 대조군의 scRMSD 중앙값
  H0: 차이 없음

실패해도 산출물로 제출한다. 음성 결과는 결함이 아니라 결과다.
```

---

## 1. 타깃 선정

### 1.1 선정 기준 (우선순위 순)

1. **본인이 도메인 맥락을 설명할 수 있을 것** — 면접 방어 가능성이 최우선
2. 단백질–단백질/항체 계열일 것 — JD 축(단백질·항체 디자인)과 일치
3. X-ray 구조, 해상도 ≤ 2.5 Å, 결측 잔기 최소
4. 총 잔기 수 ≤ 300 (복합체 기준) — Colab 제약

### 1.2 1순위: anti-CD20 VHH (나노바디)

KBioHealth 과정에서 수행한 워크플로우와 직접 연결됩니다. 현재 이력서에서 "SBDD 실습"으로 묻혀 있는 항목을 전면에 세울 수 있다는 점이 전략적으로 중요합니다.

단, **실제 PDB에 적합한 VHH–항원 복합체 구조가 존재하는지 직접 확인해야 합니다.** 제가 특정 PDB ID를 단정하지 않는 이유는, 검증하지 않은 ID를 코드에 박아 넣으면 재현성이 깨지기 때문입니다. 아래 스크립트로 조회하십시오.

### 1.3 타깃 조회 스크립트

```python
# scripts/00_select_target.py
"""
RCSB Search API로 나노바디(VHH) 복합체 후보를 조회한다.
출력: data/target_candidates.csv
"""
import requests, pandas as pd, json, pathlib

RCSB_SEARCH = "https://search.rcsb.org/rcsbsearch/v2/query"
RCSB_DATA   = "https://data.rcsb.org/rest/v1/core/entry/"

def search_candidates(keyword: str, max_res: float = 2.5, limit: int = 50):
    query = {
        "query": {
            "type": "group", "logical_operator": "and",
            "nodes": [
                {"type": "terminal", "service": "full_text",
                 "parameters": {"value": keyword}},
                {"type": "terminal", "service": "text",
                 "parameters": {
                     "attribute": "rcsb_entry_info.resolution_combined",
                     "operator": "less_or_equal", "value": max_res}},
                {"type": "terminal", "service": "text",
                 "parameters": {
                     "attribute": "rcsb_entry_info.polymer_entity_count_protein",
                     "operator": "greater_or_equal", "value": 2}},
            ],
        },
        "return_type": "entry",
        "request_options": {"paginate": {"start": 0, "rows": limit}},
    }
    r = requests.post(RCSB_SEARCH, json=query, timeout=30)
    r.raise_for_status()
    return [h["identifier"] for h in r.json().get("result_set", [])]

def fetch_meta(pdb_id: str) -> dict:
    r = requests.get(RCSB_DATA + pdb_id, timeout=30)
    r.raise_for_status()
    j = r.json()
    info = j.get("rcsb_entry_info", {})
    return {
        "pdb_id": pdb_id,
        "title": j.get("struct", {}).get("title", ""),
        "resolution": (info.get("resolution_combined") or [None])[0],
        "method": ",".join(info.get("experimental_method", []) or []),
        "n_protein_entities": info.get("polymer_entity_count_protein"),
        "deposited_residues": info.get("deposited_polymer_monomer_count"),
    }

if __name__ == "__main__":
    pathlib.Path("data").mkdir(exist_ok=True)
    ids = search_candidates("nanobody VHH complex")
    rows = []
    for i in ids:
        try:
            rows.append(fetch_meta(i))
        except Exception as e:
            print(f"[skip] {i}: {e}")
    df = pd.DataFrame(rows)
    df = df[(df.deposited_residues <= 300)].sort_values("resolution")
    df.to_csv("data/target_candidates.csv", index=False)
    print(df.head(20).to_string())
```

**선택 후 반드시 수동 확인할 것:** 구조를 PyMOL/ChimeraX로 열어 ① VHH 사슬과 항원 사슬이 실제로 접촉하는지, ② CDR 루프에 결측 잔기가 없는지. 결측 루프가 있으면 재설계 대상 영역이 비어 있다는 뜻이므로 해당 구조는 탈락입니다.

### 1.4 2순위(대체): iGEM DHCR7

VHH 복합체 확보가 지연되면 전환하십시오. 단량체 효소 재설계가 되므로 ipTM 항목이 사라지고, JD의 "항체 디자인" 축과의 거리가 멀어집니다. 시간 손실 시의 차선책입니다.

---

## 2. 디렉토리 구조 (재현성 요건)

```
vhh_mpnn_selfconsistency/
├── README.md
├── preregistration.md
├── environment.yml
├── data/
│   ├── target_candidates.csv
│   └── raw/              # RCSB 원본 .cif
├── work/
│   ├── 01_clean/         # 전처리된 PDB
│   ├── 02_mpnn/          # 설계 서열 FASTA
│   ├── 03_refold/        # 재예측 구조
│   └── 04_metrics/       # 지표 CSV
├── scripts/
│   ├── 00_select_target.py
│   ├── 01_prepare_structure.py
│   ├── 02_run_mpnn.py
│   ├── 03_refold.py
│   ├── 04_metrics.py
│   └── 05_report.py
└── figures/
```

**모든 스크립트는 `--seed` 인자를 받고, 실행 시 사용한 seed·버전·커밋 해시를 출력 CSV에 컬럼으로 기록할 것.** 이것이 "재현 가능한 산출물"과 "노트북 한 장"의 차이입니다.

---

## 3. Step 1 — 구조 전처리

### 3.1 수행 내용

1. RCSB에서 비대칭 단위 다운로드
2. 물·이온·리간드 제거 (단, 결합에 관여하는 이온이 있으면 기록 후 제거 사유 명시)
3. PDBFixer로 결측 측쇄 보완 (결측 **주쇄**가 있는 구조는 §1.3대로 애초에 탈락시킬 것)
4. VHH 사슬 / 항원 사슬 식별 및 라벨링

```python
# scripts/01_prepare_structure.py
"""
입력: PDB ID
출력: work/01_clean/{pdb_id}_clean.pdb, work/01_clean/{pdb_id}_chains.json
"""
import argparse, json, pathlib, requests
from pdbfixer import PDBFixer
from openmm.app import PDBFile
from Bio.PDB import PDBParser, PPBuilder

def download(pdb_id: str, out_dir: pathlib.Path) -> pathlib.Path:
    out_dir.mkdir(parents=True, exist_ok=True)
    dst = out_dir / f"{pdb_id}.pdb"
    if not dst.exists():
        r = requests.get(f"https://files.rcsb.org/download/{pdb_id}.pdb", timeout=60)
        r.raise_for_status()
        dst.write_bytes(r.content)
    return dst

def clean(src: pathlib.Path, dst: pathlib.Path) -> None:
    fixer = PDBFixer(filename=str(src))
    fixer.removeHeterogens(keepWater=False)
    fixer.findMissingResidues()
    # 말단 결측은 보완하지 않는다(인위적 꼬리 생성 방지)
    chains = list(fixer.topology.chains())
    keys = list(fixer.missingResidues.keys())
    for key in keys:
        ch = chains[key[0]]
        if key[1] == 0 or key[1] == len(list(ch.residues())):
            del fixer.missingResidues[key]
    fixer.findMissingAtoms()
    fixer.addMissingAtoms()
    dst.parent.mkdir(parents=True, exist_ok=True)
    PDBFile.writeFile(fixer.topology, fixer.positions, open(dst, "w"))

def summarize_chains(pdb_path: pathlib.Path) -> dict:
    st = PDBParser(QUIET=True).get_structure("s", str(pdb_path))
    ppb = PPBuilder()
    out = {}
    for ch in st[0]:
        seq = "".join(str(pp.get_sequence()) for pp in ppb.build_peptides(ch))
        if seq:
            out[ch.id] = {"length": len(seq), "sequence": seq}
    return out

if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--pdb_id", required=True)
    a = ap.parse_args()
    raw = download(a.pdb_id, pathlib.Path("data/raw"))
    clean_path = pathlib.Path(f"work/01_clean/{a.pdb_id}_clean.pdb")
    clean(raw, clean_path)
    chains = summarize_chains(clean_path)
    json.dump(chains, open(f"work/01_clean/{a.pdb_id}_chains.json", "w"), indent=2)
    for cid, v in chains.items():
        print(f"chain {cid}: {v['length']} aa")
```

### 3.2 CDR 위치 지정

ProteinMPNN에 "어디를 재설계할지" 알려야 합니다. 두 가지 방식이 있고, **둘 다 수행해 비교하는 것을 권합니다.**

| 설계 모드 | 고정(fixed) 영역 | 재설계 영역 | 의미 |
|---|---|---|---|
| **M1 전체 재설계** | 항원 전체 | VHH 전체 | inverse folding 기본 성능 측정 |
| **M2 CDR 한정** | 항원 + VHH 프레임워크 | CDR1/2/3만 | 항체 공학에 가까운 현실적 설정 |

CDR 경계는 ANARCI(IMGT numbering)로 자동 산출하는 것이 수작업보다 재현성이 높습니다. ANARCI 설치가 막히면, 구조상 항원과 4.5 Å 이내 접촉하는 VHH 잔기를 "인터페이스 잔기"로 정의해 대체하고, **CDR이 아니라 인터페이스 정의를 썼다는 사실을 README에 명시**하십시오. 정의를 바꾼 것을 숨기는 것이 가장 나쁜 선택입니다.

```python
# scripts/01b_define_designable.py (발췌)
from Bio.PDB import PDBParser, NeighborSearch

def interface_residues(pdb, design_chain, partner_chain, cutoff=4.5):
    st = PDBParser(QUIET=True).get_structure("s", pdb)[0]
    partner_atoms = [a for a in st[partner_chain].get_atoms()]
    ns = NeighborSearch(partner_atoms)
    hits = []
    for res in st[design_chain]:
        if any(ns.search(a.coord, cutoff) for a in res):
            hits.append(res.id[1])
    return sorted(hits)
```

---

## 4. Step 2 — ProteinMPNN 서열 재설계

### 4.1 설계 매트릭스

온도(temperature)는 서열 다양성을 좌우하는 핵심 파라미터입니다. 단일 값만 쓰면 "파라미터를 이해하고 썼는가"라는 면접 질문에 답할 근거가 없습니다.

| 조건 | temperature | 설계 서열 수 | 목적 |
|---|---|---|---|
| T01 | 0.1 | 8 | 보수적 설계, 높은 복원율 기대 |
| T02 | 0.2 | 8 | 중간 |
| T03 | 0.3 | 8 | 다양성 확대, scRMSD 악화 여부 관찰 |

모드(M1/M2) × 온도 3수준 × 8서열 = **48 서열**. 이 중 재예측은 비용 때문에 일부만 수행합니다(§5.1).

### 4.2 실행

ProteinMPNN은 공개 레포(`dauparas/ProteinMPNN`)를 클론해 사용합니다. 라이선스(MIT 계열로 알려져 있으나 **클론 후 LICENSE 파일을 직접 확인**할 것)를 README에 기재하십시오.

```bash
# scripts/02_run_mpnn.sh
set -euo pipefail
PDB=${1:?pdb path}
OUT=${2:?out dir}
DESIGN_CHAIN=${3:?e.g. B}
FIXED_JSON=${4:-""}      # M2일 때만 사용
SEED=${5:-42}

MPNN=third_party/ProteinMPNN
mkdir -p "$OUT"

python $MPNN/helper_scripts/parse_multiple_chains.py \
  --input_path "$(dirname "$PDB")" --output_path "$OUT/parsed.jsonl"

python $MPNN/helper_scripts/assign_fixed_chains.py \
  --input_path "$OUT/parsed.jsonl" --output_path "$OUT/chains.jsonl" \
  --chain_list "$DESIGN_CHAIN"

EXTRA=""
if [ -n "$FIXED_JSON" ]; then EXTRA="--fixed_positions_jsonl $FIXED_JSON"; fi

for T in 0.1 0.2 0.3; do
  python $MPNN/protein_mpnn_run.py \
    --jsonl_path "$OUT/parsed.jsonl" \
    --chain_id_jsonl "$OUT/chains.jsonl" \
    $EXTRA \
    --out_folder "$OUT/T${T}" \
    --num_seq_per_target 8 \
    --sampling_temp "$T" \
    --seed "$SEED" \
    --batch_size 1
done
```

### 4.3 음성 대조군 (생략 불가)

self-consistency 수치는 **대조군 없이는 해석 불가능**합니다. 설계 서열의 scRMSD가 1.5 Å라고 해도, 아무 서열이나 넣어도 1.5 Å가 나온다면 아무것도 입증하지 못합니다.

최소 2종을 생성하십시오.

```python
# scripts/02b_controls.py (발췌)
import random

def shuffled_control(seq: str, designable: list[int], seed: int) -> str:
    """재설계 대상 위치의 잔기만 서로 섞는다 (조성 보존, 배치 무작위)"""
    rng = random.Random(seed)
    chars = list(seq)
    idx = [i - 1 for i in designable]          # 1-based -> 0-based
    pool = [chars[i] for i in idx]
    rng.shuffle(pool)
    for i, aa in zip(idx, pool):
        chars[i] = aa
    return "".join(chars)

def random_control(seq: str, designable: list[int], seed: int) -> str:
    """재설계 대상 위치를 20종 균등 무작위로 치환"""
    AA = "ACDEFGHIKLMNPQRSTVWY"
    rng = random.Random(seed)
    chars = list(seq)
    for i in designable:
        chars[i - 1] = rng.choice(AA)
    return "".join(chars)
```

- **대조군 1 (shuffle):** 아미노산 조성은 보존되고 배치만 무작위. 가장 엄격한 대조군.
- **대조군 2 (random):** 완전 무작위. 상한 참조값.
- **양성 대조군:** 원본(native) 서열 그대로 재예측. 예측기 자체의 오차 하한을 제공합니다. **이것을 빠뜨리면 scRMSD가 예측기 오차인지 설계 실패인지 구분할 수 없습니다.**

각 대조군 8개씩.

---

## 5. Step 3 — 재예측 (refolding)

### 5.1 수행 대상 축소

총 서열: 설계 48 + 대조군 24 + native 1 = 73. 전부 AF2-multimer로 돌리면 무료 티어에서 시간이 부족합니다.

**2단계 전략을 쓰십시오.**

| 단계 | 예측기 | 대상 | 목적 |
|---|---|---|---|
| S1 (스크리닝) | ESMFold — VHH 단량체만 | 73개 전부 | 빠른 1차 선별. 단일 서열 입력, MSA 불필요 |
| S2 (정밀) | ColabFold AF2-multimer | S1 상위 설계 6 + 대조군 6 + native 1 = 13 | 복합체 ipTM 산출 |

ESMFold는 MSA를 쓰지 않으므로 설계 서열(자연계에 상동체가 없음)에 특히 적합합니다. AF2는 MSA 의존도가 높아 설계 서열에서 성능이 떨어질 수 있는데, **이것 자체가 보고할 만한 관찰 결과입니다.**

### 5.2 ColabFold 실행 (S2)

```python
# scripts/03_refold.py (ColabFold batch 모드)
"""
입력: work/02_mpnn/selected.fasta  (복합체는 사슬을 ':'로 연결)
출력: work/03_refold/{name}/*.pdb + *_scores.json
"""
import subprocess, pathlib, argparse

def run_colabfold(fasta: pathlib.Path, out: pathlib.Path,
                  num_recycle: int = 3, models: int = 1, seed: int = 42):
    out.mkdir(parents=True, exist_ok=True)
    cmd = [
        "colabfold_batch", str(fasta), str(out),
        "--num-recycle", str(num_recycle),
        "--num-models", str(models),
        "--random-seed", str(seed),
        "--model-type", "alphafold2_multimer_v3",
    ]
    subprocess.run(cmd, check=True)
```

**복합체 FASTA 형식:** `>design_T01_s1` 다음 줄에 `항원서열:VHH설계서열` — 콜론이 사슬 구분자입니다. 사슬 순서를 전 샘플에서 동일하게 유지하십시오.

`--num-models 1`은 시간 절약을 위한 타협입니다. 모델 1개만 쓰면 예측 불확실성을 과소평가하게 되므로, **이 타협을 README 한계 항목에 적으십시오.**

---

## 6. Step 4 — 지표 산출

### 6.1 측정 항목

| 지표 | 정의 | 출처 |
|---|---|---|
| `scRMSD` | 재예측 구조 vs 원본 backbone, VHH 사슬 Cα | US-align / TM-align |
| `TM-score` | 길이 정규화 구조 유사도 | US-align |
| `pLDDT_mean` | 재예측 구조 전체 평균 | ColabFold scores.json |
| `pLDDT_cdr` | 재설계 영역만의 평균 | B-factor 컬럼에서 추출 |
| `ipTM` | 사슬 간 계면 신뢰도 | ColabFold scores.json (multimer 전용) |
| `seq_recovery` | 원본 서열과의 동일 잔기 비율 (재설계 영역 한정) | 직접 계산 |

`seq_recovery`는 **설계 품질 지표가 아닙니다.** 원본을 그대로 복원하면 recovery는 100%지만 설계로서의 가치는 0입니다. "복원율이 높을수록 좋다"고 쓰면 면접에서 바로 지적당합니다. recovery는 온도 파라미터가 의도대로 작동했는지 확인하는 **진단 지표**로만 쓰십시오.

### 6.2 산출 스크립트

```python
# scripts/04_metrics.py
"""
출력: work/04_metrics/metrics.csv
컬럼: name, group, mode, temperature, seed, scRMSD, tm_score,
      pLDDT_mean, pLDDT_cdr, ipTM, seq_recovery, tool_version, git_commit
"""
import json, re, subprocess, pathlib, pandas as pd
from Bio.PDB import PDBParser

def us_align(pred: str, ref: str, usalign_bin="USalign") -> dict:
    out = subprocess.run([usalign_bin, pred, ref], capture_output=True,
                         text=True, check=True).stdout
    rmsd = float(re.search(r"RMSD=\s*([\d.]+)", out).group(1))
    tm   = float(re.search(r"TM-score=\s*([\d.]+)", out).group(1))
    return {"scRMSD": rmsd, "tm_score": tm}

def plddt_from_pdb(pdb: str, residues: list[int] | None = None) -> float:
    st = PDBParser(QUIET=True).get_structure("s", pdb)[0]
    vals = [a.get_bfactor() for ch in st for res in ch
            if residues is None or res.id[1] in residues
            for a in res if a.get_id() == "CA"]
    return sum(vals) / len(vals) if vals else float("nan")

def seq_recovery(design: str, native: str, positions: list[int]) -> float:
    same = sum(design[i-1] == native[i-1] for i in positions)
    return same / len(positions)

def git_commit() -> str:
    try:
        return subprocess.run(["git", "rev-parse", "--short", "HEAD"],
                              capture_output=True, text=True).stdout.strip()
    except Exception:
        return "NA"
```

### 6.3 Boltz-2 (선택 단계 — 실패를 허용할 것)

AF2-multimer와 다른 아키텍처로 교차 검증하면 산출물의 설득력이 올라갑니다. 다만 §0.1 A4대로 무료 티어 VRAM이 불확실합니다.

**중단 규칙:** Boltz-2 환경 구축에 **90분**을 넘기면 중단하고 "AF2-multimer 단일 예측기로만 검증함 — 교차 검증은 미수행"을 한계로 기록하십시오. 산출물 없이 마감을 넘기는 것이 최악의 결과입니다.

---

## 7. 계산 자원 실패 시 대체 경로

| 장애 | 대응 |
|---|---|
| Colab GPU 미할당 | ProteinMPNN은 CPU로 진행(A2). 재예측만 GPU 세션을 재시도 |
| AF2-multimer OOM | ① `--num-recycle 1`로 축소 → ② VHH 단량체 폴딩으로 전환(ipTM 포기, scRMSD/pLDDT만 보고) |
| ESMFold 설치 실패 | S1 스크리닝을 생략하고 S2 대상을 13개로 직접 선정 |
| 타깃 구조 확보 실패 | §1.4 DHCR7 단량체로 전환 |

**어떤 경로를 탔든 README에 실제 경로를 적으십시오.** 축소된 설계를 원래 계획인 것처럼 쓰면 코드와 대조했을 때 드러납니다.

---

## 8. 3일 작업 분배

| 일자 | 작업 | 산출 | 중단 기준 |
|---|---|---|---|
| **D1 오전** | 타깃 조회·선정, 구조 육안 검증 | `target_candidates.csv`, PDB 1건 확정 | 3시간 내 미확정 시 DHCR7 전환 |
| **D1 오후** | 전처리, CDR/인터페이스 정의, `preregistration.md` 커밋 | `*_clean.pdb`, `designable.json` | — |
| **D2 오전** | ProteinMPNN M1/M2 × 3온도 실행, 대조군 생성 | 73 서열 FASTA | — |
| **D2 오후** | S1 ESMFold 스크리닝 → S2 대상 13건 선정 | `screening.csv` | ESMFold 2시간 초과 시 생략 |
| **D3 오전** | ColabFold S2 실행 | 재예측 구조 13건 | — |
| **D3 오후** | 지표 산출, 그림 2종, README 작성 | `metrics.csv`, `figures/`, README | — |

Boltz-2는 일정에 넣지 않았습니다. D3 오후에 여유가 생길 때만 착수하십시오.

---

## 9. 보고 형식

### 9.1 그림 2종 (그 이상은 불필요)

- **Fig 1:** 그룹별 scRMSD 분포 — native / design(T01,T02,T03) / shuffle / random. 박스플롯 + 개별 점 오버레이. N이 작으므로 점을 반드시 보일 것.
- **Fig 2:** scRMSD vs pLDDT 산점도, 그룹별 색상. 두 지표가 일치하는지/어긋나는지를 보여줍니다.

```python
# scripts/05_report.py (발췌)
import pandas as pd, matplotlib.pyplot as plt, numpy as np

df = pd.read_csv("work/04_metrics/metrics.csv")
order = ["native", "design_T0.1", "design_T0.2", "design_T0.3", "shuffle", "random"]

fig, ax = plt.subplots(figsize=(7, 4))
data = [df.loc[df.group == g, "scRMSD"].dropna().values for g in order]
ax.boxplot(data, labels=order, showfliers=False)
for i, v in enumerate(data, start=1):
    ax.scatter(np.random.normal(i, 0.05, len(v)), v, s=18, alpha=0.7, zorder=3)
ax.axhline(2.0, ls="--", lw=1, color="grey")   # 사전 등록 임계값
ax.set_ylabel("scRMSD to native backbone (Å)")
ax.set_xlabel(f"group (n per group shown as points)")
plt.xticks(rotation=20); plt.tight_layout()
plt.savefig("figures/fig1_scrmsd.png", dpi=200)
```

### 9.2 README 필수 섹션

```markdown
# VHH inverse-folding self-consistency benchmark

## 무엇을 했는가
## 무엇을 하지 않았는가   ← §0.2를 그대로
## 방법 (버전 고정 표)
## 사전 등록된 기준       ← preregistration.md 링크
## 결과 (표 1 + Fig 1,2)
## 한계
## 재현 방법 (명령어 3줄)
```

**한계 섹션에 반드시 포함:**

1. 실험 검증 없음. in silico self-consistency는 발현성·결합력의 대리 지표가 아님
2. N이 조건당 8로 작아 기술통계 수준의 비교만 가능
3. AF2 모델 1개·recycle 3회로 축소 — 예측 불확실성 과소평가 가능
4. 단일 타깃. 일반화 주장 불가
5. ProteinMPNN은 일반 단백질 코퍼스로 학습되었고, 항체 CDR 루프는 서열-구조 대응이 약해 범용 inverse folding 모델의 성능이 떨어진다는 보고가 있음(AbMPNN 등 항체 특화 모델이 별도로 제안된 배경). **이 서술은 제가 인용을 교차 확인하지 않았으므로, 제출 전 1차 문헌으로 직접 확인하고 인용을 달거나, 확인이 어려우면 문장을 삭제하십시오.** 다만 본인 결과에서 M2(CDR 한정)가 M1보다 나쁘게 나온다면, 그것은 본인 데이터로 직접 말할 수 있는 관찰입니다

---

## 10. 이력서·면접 반영

### 10.1 이력서 기재 문안 (수행 결과 확정 후 수치 삽입)

```
단백질 설계 self-consistency 벤치마크 | 2026.10 | 개인 프로젝트
• VHH–항원 복합체를 대상으로 ProteinMPNN inverse folding 설계 서열 48종 생성
  (전체 재설계 / CDR 한정 2모드 × 샘플링 온도 3수준)
• AF2-multimer 재예측으로 scRMSD·pLDDT·ipTM 산출, 조성 보존 셔플 및
  무작위 치환 음성 대조군, native 양성 대조군과 비교
• 결과: [실측값 삽입]
• 한계: 실험 검증 없음, 조건당 n=8, 단일 타깃, 예측 모델 1개
• 기술: ProteinMPNN, ColabFold(AF2-multimer), ESMFold, US-align, PDBFixer, Biopython
• GitHub: [링크]
```

스킬 섹션에는 **실제로 실행한 것만** 추가하십시오: ProteinMPNN, ColabFold, ESMFold, US-align. RFdiffusion은 실행하지 않았다면 넣지 말고, 자기소개서에서 "학습 중"으로 분리 서술하십시오.

### 10.2 예상 기술질문과 답변 방향

| 질문 | 답변의 핵심 |
|---|---|
| 왜 self-consistency가 설계 품질의 지표가 되나? | 완전한 지표가 아님을 먼저 인정. inverse folding이 생성한 서열이 독립 예측기에서 원본 fold로 되돌아오는지를 보는 필요조건이며, 충분조건이 아님 |
| 대조군을 왜 두 종 뒀나? | shuffle은 조성 효과를 통제, random은 상한 참조. native는 예측기 자체 오차 하한 제공 |
| temperature를 왜 3수준으로? | 다양성–복원율 트레이드오프를 한 점이 아니라 곡선으로 보기 위함. seq_recovery로 파라미터 작동 확인 |
| 이 결과로 실제 결합을 주장할 수 있나? | 불가. 다음 단계는 ΔΔG 계산 또는 display 기반 실험 검증 |
| RFdiffusion은 왜 안 썼나? | 6일 일정에서 유의미한 결과를 낼 수 없다고 판단. 백본 생성은 평가 지표 설계가 더 복잡해 별도 과제로 분리 |

마지막 질문에 대한 답변은 **시간 제약을 솔직히 말하는 것이 가장 강한 답**입니다. 범위를 스스로 잘라낸 판단은 연구자로서의 감점 요소가 아닙니다.
