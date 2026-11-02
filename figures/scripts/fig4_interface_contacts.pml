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

# Fig 4 — native 복합체의 계면. 접촉 정의는 designable.json과 동일한
# 4.5 A 기준이며, 잔기 번호를 코드에 박지 않고 그 파일에서 읽었습니다.
load /mnt/c/Users/PC/Downloads/exercise/work/01_clean/8YVJ_clean.pdb, cx

hide everything
show cartoon, cx and chain A
color grey80, cx and chain A

# 항원은 표면으로 — 포켓 형상과 리간드 사슬의 맞물림을 보기 위함
show surface, cx and chain B
color palecyan, cx and chain B
set transparency, 0.35

# 접촉 잔기(4.5 A 이내)와 CDR3를 구분해서 표시
show sticks, cx and chain A and resi 29+30+31+32+33+34+35+56+57+58+75+77+78+79+80 and not hydro
color tv_orange, cx and chain A and resi 29+30+31+32+33+34+35+56+57+58+75+77+78+79+80 and elem C
color salmon, cx and chain A and resi 100+101+102+103+104+105+106+107+108+109+110+111+112+113+114+115+116+117+118

# 극성 접촉 — 수소결합 후보. 거리 3.5 A 기준은 PyMOL 기본값입니다.
distance polar_contacts, \
    cx and chain A and resi 29+30+31+32+33+34+35+56+57+58+75+77+78+79+80, \
    cx and chain B, 3.5, mode=2
# 거리 라벨은 끕니다. 접촉이 15잔기라 숫자가 겹쳐 계면을 덮습니다.
# 개수와 위치만 보이면 되고, 정확한 거리는 Rosetta 쪽 hbonds_int로 셉니다.
hide labels, polar_contacts
color grey30, polar_contacts
set dash_width, 2.0
set dash_gap, 0.3

# 전체가 아니라 계면을 중심으로 잡습니다. 복합체 전체를 orient하면
# 항원의 beta-solenoid가 화면을 차지하고 계면은 작아집니다.
orient cx and chain A and resi 29+30+31+32+33+34+35+56+57+58+75+77+78+79+80
zoom cx and chain A and resi 29+30+31+32+33+34+35+56+57+58+75+77+78+79+80, 16
turn x, -10
png /mnt/c/Users/PC/Downloads/exercise/figures/fig4_interface_contacts.png, width=1800, height=1350, dpi=300, ray=1
