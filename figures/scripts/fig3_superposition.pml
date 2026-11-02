# 자동 생성 — scripts/06_render.py. 직접 수정하지 마십시오.
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

# Fig 3 — native 사슬 A vs 설계 M2_T0.3_s2 (ESMFold S1 예측)
# 중첩 기준: 전체 Cα. 셀 [40]의 scRMSD와 같은 정의라 그림과 수치가 대응합니다.
# (참고: 프레임워크만으로 맞추면 CDR3 차이가 더 커 보이지만, 보고한 수치와
#  어긋납니다. Fig 5는 루프 비교가 목적이라 프레임워크 기준을 씁니다.)
load /mnt/c/Users/PC/Downloads/exercise/work/01_clean/8YVJ_clean.pdb, native
load /mnt/c/Users/PC/Downloads/exercise/work/03_refold/S1_esmfold/M2_T0.3_s2.pdb, design

# 항원 사슬은 이 그림의 주제가 아니므로 제거
remove native and not chain A

hide everything
show cartoon
color grey70, native
color skyblue, design
color orange, native and resi 100+101+102+103+104+105+106+107+108+109+110+111+112+113+114+115+116+117+118
color salmon, design and resi 100+101+102+103+104+105+106+107+108+109+110+111+112+113+114+115+116+117+118

# 재설계 영역에 포함된 Cys — 나노바디 CDR1-CDR3 추가 이황화결합의 파트너.
# 보존 여부가 README 결과 (5)의 근거이므로 눈으로도 확인할 수 있게 둡니다.
show sticks, (native or design) and resi 110 and not hydro
color yellow, (native or design) and resi 110 and elem S

# scRMSD = 1.21 A (셀 [40] metrics.csv 기준)
super design, native
orient native
turn y, 20

png /mnt/c/Users/PC/Downloads/exercise/figures/fig3_superposition.png, width=1800, height=1350, dpi=300, ray=1
