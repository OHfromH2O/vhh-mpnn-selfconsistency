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

# Fig 5 — 온도별 CDR3 백본 앙상블 (ESMFold S1 예측)
# 중첩 기준: **프레임워크(CDR3 제외) Cα**. 루프끼리 비교하는 것이 목적이므로
# 전체 Cα로 맞추면 루프 차이가 프레임워크 쪽으로 흡수되어 보이지 않습니다.
# Fig 3와 기준이 다른 이유가 이것입니다.
#
# 순서 주의: load -> super 를 모두 끝낸 뒤에 hide everything 을 합니다.
# 먼저 숨기면 이후 로드된 객체가 기본 cartoon으로 되살아납니다.
load /mnt/c/Users/PC/Downloads/exercise/work/01_clean/8YVJ_clean.pdb, native
remove native and not chain A

load /mnt/c/Users/PC/Downloads/exercise/work/03_refold/S1_esmfold/M2_T0.1_s8.pdb, M2_T0_1_s8
super M2_T0_1_s8 and not resi 100+101+102+103+104+105+106+107+108+109+110+111+112+113+114+115+116+117+118, native and not resi 100+101+102+103+104+105+106+107+108+109+110+111+112+113+114+115+116+117+118
load /mnt/c/Users/PC/Downloads/exercise/work/03_refold/S1_esmfold/M2_T0.1_s7.pdb, M2_T0_1_s7
super M2_T0_1_s7 and not resi 100+101+102+103+104+105+106+107+108+109+110+111+112+113+114+115+116+117+118, native and not resi 100+101+102+103+104+105+106+107+108+109+110+111+112+113+114+115+116+117+118
load /mnt/c/Users/PC/Downloads/exercise/work/03_refold/S1_esmfold/M1_T0.1_s4.pdb, M1_T0_1_s4
super M1_T0_1_s4 and not resi 100+101+102+103+104+105+106+107+108+109+110+111+112+113+114+115+116+117+118, native and not resi 100+101+102+103+104+105+106+107+108+109+110+111+112+113+114+115+116+117+118
load /mnt/c/Users/PC/Downloads/exercise/work/03_refold/S1_esmfold/M2_T0.1_s2.pdb, M2_T0_1_s2
super M2_T0_1_s2 and not resi 100+101+102+103+104+105+106+107+108+109+110+111+112+113+114+115+116+117+118, native and not resi 100+101+102+103+104+105+106+107+108+109+110+111+112+113+114+115+116+117+118
load /mnt/c/Users/PC/Downloads/exercise/work/03_refold/S1_esmfold/M2_T0.2_s7.pdb, M2_T0_2_s7
super M2_T0_2_s7 and not resi 100+101+102+103+104+105+106+107+108+109+110+111+112+113+114+115+116+117+118, native and not resi 100+101+102+103+104+105+106+107+108+109+110+111+112+113+114+115+116+117+118
load /mnt/c/Users/PC/Downloads/exercise/work/03_refold/S1_esmfold/M1_T0.2_s4.pdb, M1_T0_2_s4
super M1_T0_2_s4 and not resi 100+101+102+103+104+105+106+107+108+109+110+111+112+113+114+115+116+117+118, native and not resi 100+101+102+103+104+105+106+107+108+109+110+111+112+113+114+115+116+117+118
load /mnt/c/Users/PC/Downloads/exercise/work/03_refold/S1_esmfold/M1_T0.2_s1.pdb, M1_T0_2_s1
super M1_T0_2_s1 and not resi 100+101+102+103+104+105+106+107+108+109+110+111+112+113+114+115+116+117+118, native and not resi 100+101+102+103+104+105+106+107+108+109+110+111+112+113+114+115+116+117+118
load /mnt/c/Users/PC/Downloads/exercise/work/03_refold/S1_esmfold/M2_T0.2_s1.pdb, M2_T0_2_s1
super M2_T0_2_s1 and not resi 100+101+102+103+104+105+106+107+108+109+110+111+112+113+114+115+116+117+118, native and not resi 100+101+102+103+104+105+106+107+108+109+110+111+112+113+114+115+116+117+118
load /mnt/c/Users/PC/Downloads/exercise/work/03_refold/S1_esmfold/M2_T0.3_s7.pdb, M2_T0_3_s7
super M2_T0_3_s7 and not resi 100+101+102+103+104+105+106+107+108+109+110+111+112+113+114+115+116+117+118, native and not resi 100+101+102+103+104+105+106+107+108+109+110+111+112+113+114+115+116+117+118
load /mnt/c/Users/PC/Downloads/exercise/work/03_refold/S1_esmfold/M1_T0.3_s4.pdb, M1_T0_3_s4
super M1_T0_3_s4 and not resi 100+101+102+103+104+105+106+107+108+109+110+111+112+113+114+115+116+117+118, native and not resi 100+101+102+103+104+105+106+107+108+109+110+111+112+113+114+115+116+117+118
load /mnt/c/Users/PC/Downloads/exercise/work/03_refold/S1_esmfold/M2_T0.3_s4.pdb, M2_T0_3_s4
super M2_T0_3_s4 and not resi 100+101+102+103+104+105+106+107+108+109+110+111+112+113+114+115+116+117+118, native and not resi 100+101+102+103+104+105+106+107+108+109+110+111+112+113+114+115+116+117+118
load /mnt/c/Users/PC/Downloads/exercise/work/03_refold/S1_esmfold/M1_T0.3_s7.pdb, M1_T0_3_s7
super M1_T0_3_s7 and not resi 100+101+102+103+104+105+106+107+108+109+110+111+112+113+114+115+116+117+118, native and not resi 100+101+102+103+104+105+106+107+108+109+110+111+112+113+114+115+116+117+118

hide everything
show ribbon, M2_T0_1_s8 and resi 100+101+102+103+104+105+106+107+108+109+110+111+112+113+114+115+116+117+118
color marine, M2_T0_1_s8 and resi 100+101+102+103+104+105+106+107+108+109+110+111+112+113+114+115+116+117+118
show ribbon, M2_T0_1_s7 and resi 100+101+102+103+104+105+106+107+108+109+110+111+112+113+114+115+116+117+118
color marine, M2_T0_1_s7 and resi 100+101+102+103+104+105+106+107+108+109+110+111+112+113+114+115+116+117+118
show ribbon, M1_T0_1_s4 and resi 100+101+102+103+104+105+106+107+108+109+110+111+112+113+114+115+116+117+118
color marine, M1_T0_1_s4 and resi 100+101+102+103+104+105+106+107+108+109+110+111+112+113+114+115+116+117+118
show ribbon, M2_T0_1_s2 and resi 100+101+102+103+104+105+106+107+108+109+110+111+112+113+114+115+116+117+118
color marine, M2_T0_1_s2 and resi 100+101+102+103+104+105+106+107+108+109+110+111+112+113+114+115+116+117+118
show ribbon, M2_T0_2_s7 and resi 100+101+102+103+104+105+106+107+108+109+110+111+112+113+114+115+116+117+118
color forest, M2_T0_2_s7 and resi 100+101+102+103+104+105+106+107+108+109+110+111+112+113+114+115+116+117+118
show ribbon, M1_T0_2_s4 and resi 100+101+102+103+104+105+106+107+108+109+110+111+112+113+114+115+116+117+118
color forest, M1_T0_2_s4 and resi 100+101+102+103+104+105+106+107+108+109+110+111+112+113+114+115+116+117+118
show ribbon, M1_T0_2_s1 and resi 100+101+102+103+104+105+106+107+108+109+110+111+112+113+114+115+116+117+118
color forest, M1_T0_2_s1 and resi 100+101+102+103+104+105+106+107+108+109+110+111+112+113+114+115+116+117+118
show ribbon, M2_T0_2_s1 and resi 100+101+102+103+104+105+106+107+108+109+110+111+112+113+114+115+116+117+118
color forest, M2_T0_2_s1 and resi 100+101+102+103+104+105+106+107+108+109+110+111+112+113+114+115+116+117+118
show ribbon, M2_T0_3_s7 and resi 100+101+102+103+104+105+106+107+108+109+110+111+112+113+114+115+116+117+118
color firebrick, M2_T0_3_s7 and resi 100+101+102+103+104+105+106+107+108+109+110+111+112+113+114+115+116+117+118
show ribbon, M1_T0_3_s4 and resi 100+101+102+103+104+105+106+107+108+109+110+111+112+113+114+115+116+117+118
color firebrick, M1_T0_3_s4 and resi 100+101+102+103+104+105+106+107+108+109+110+111+112+113+114+115+116+117+118
show ribbon, M2_T0_3_s4 and resi 100+101+102+103+104+105+106+107+108+109+110+111+112+113+114+115+116+117+118
color firebrick, M2_T0_3_s4 and resi 100+101+102+103+104+105+106+107+108+109+110+111+112+113+114+115+116+117+118
show ribbon, M1_T0_3_s7 and resi 100+101+102+103+104+105+106+107+108+109+110+111+112+113+114+115+116+117+118
color firebrick, M1_T0_3_s7 and resi 100+101+102+103+104+105+106+107+108+109+110+111+112+113+114+115+116+117+118

# native: 프레임워크는 흐린 cartoon으로 맥락만, CDR3는 검은 굵은 선으로 기준선
show cartoon, native and not resi 100+101+102+103+104+105+106+107+108+109+110+111+112+113+114+115+116+117+118
color grey80, native and not resi 100+101+102+103+104+105+106+107+108+109+110+111+112+113+114+115+116+117+118
set cartoon_transparency, 0.75, native
show ribbon, native and resi 100+101+102+103+104+105+106+107+108+109+110+111+112+113+114+115+116+117+118
color black, native and resi 100+101+102+103+104+105+106+107+108+109+110+111+112+113+114+115+116+117+118
set ribbon_width, 5
set ribbon_sampling, 2

orient native and resi 100+101+102+103+104+105+106+107+108+109+110+111+112+113+114+115+116+117+118
zoom native and resi 100+101+102+103+104+105+106+107+108+109+110+111+112+113+114+115+116+117+118, 8
png /mnt/c/Users/PC/Downloads/exercise/figures/fig5_cdr3_ensemble.png, width=1800, height=1350, dpi=300, ray=1
