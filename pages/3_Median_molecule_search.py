import streamlit as st
from streamlit_ketcher import st_ketcher

from app_utils import (
    get_default_selfies_alphabet,
    median_molecule_search,
    smiles_list_to_grid_image,
    smiles_to_mol,
)

st.title("3) Median molecule search (multiple references)")

st.markdown(
    """
여러 개(최소 2개)의 참조 분자에 동시에 유사한 median-like 분자를 SELFIES 변이로 탐색합니다.
"""
)

random_seed = st.sidebar.number_input("Random seed", min_value=0, value=0, step=1)
allow_charged = st.sidebar.checkbox("생성된 분자에서 이온(+ / -) 허용", value=False)

alphabet = get_default_selfies_alphabet()
st.write(f"SELFIES alphabet size: **{len(alphabet)}** (기본 alphabet)")

median_input_mode = st.radio(
    "Reference 분자 입력 방식",
    ["SMILES 텍스트", "구조 그리기"],
    horizontal=True,
)

if "median_ref_smiles_drawn" not in st.session_state:
    st.session_state["median_ref_smiles_drawn"] = []
if "median_df" not in st.session_state:
    st.session_state["median_df"] = None

ref_smiles = []

if median_input_mode == "SMILES 텍스트":
    ref_smiles_input = st.text_area(
        "Reference SMILES 목록 (한 줄에 하나씩)",
        value="c1ccccc1\nc1ccccc1O\nc1ccccc1N",
        height=120,
    )
    ref_smiles = [s.strip() for s in ref_smiles_input.splitlines() if s.strip()]
else:
    st.caption("구조 수정 후 Ketcher 안의 Apply(✔)를 눌러야 값이 반영됩니다.")

    smiles_ref_drawn = st_ketcher(key="median_ref_ketcher")

    col_d1, col_d2 = st.columns(2)
    with col_d1:
        if st.button("현재 그린 구조를 ref 목록에 추가"):
            if smiles_ref_drawn and smiles_ref_drawn.strip():
                smi = smiles_ref_drawn.strip()
                if smiles_to_mol(smi) is None:
                    st.error("그려진 구조(SMILES)를 RDKit에서 해석하지 못했습니다.")
                else:
                    st.session_state["median_ref_smiles_drawn"].append(smi)
            else:
                st.warning("먼저 구조를 그리고 Ketcher 안에서 Apply(✔)를 누른 뒤 다시 시도해 주세요.")
    with col_d2:
        if st.button("ref 목록 비우기"):
            st.session_state["median_ref_smiles_drawn"] = []

    if st.session_state["median_ref_smiles_drawn"]:
        st.write("현재 reference SMILES 목록:")
        for i, smi in enumerate(st.session_state["median_ref_smiles_drawn"], start=1):
            st.write(f"{i}. `{smi}`")
    else:
        st.info("아직 추가된 reference가 없습니다.")

    ref_smiles = list(st.session_state["median_ref_smiles_drawn"])

c4, c5, c6 = st.columns(3)
with c4:
    n_starts = st.number_input(
        "시작점 개수 (ref에서 선택)",
        min_value=1,
        max_value=10,
        value=3,
        step=1,
        key="median_n_starts",
    )
with c5:
    steps_per_start = st.number_input(
        "각 시작점당 최대 스텝 수",
        min_value=1,
        max_value=100,
        value=20,
        step=1,
        key="median_steps_per_start",
    )
with c6:
    neighbors_per_step = st.number_input(
        "각 스텝당 후보 개수",
        min_value=8,
        max_value=512,
        value=128,
        step=8,
        key="median_neighbors_per_step",
    )

use_hill_climb = st.checkbox(
    "Hill-climb 모드 사용 (개선되지 않으면 해당 시작점에서 탐색 중단)",
    value=False,
)

if st.button("Median molecule 탐색", type="primary"):
    if len(ref_smiles) < 2:
        st.error("최소 2개 이상의 reference 분자가 필요합니다.")
    else:
        with st.spinner("SELFIES 공간에서 median molecules 탐색 중..."):
            st.session_state["median_df"] = median_molecule_search(
                smiles_list=ref_smiles,
                alphabet=alphabet,
                n_starts=int(n_starts),
                steps_per_start=int(steps_per_start),
                neighbors_per_step=int(neighbors_per_step),
                use_hill_climb=use_hill_climb,
                random_seed=int(random_seed),
                allow_charged=allow_charged,
            )

if st.session_state["median_df"] is not None:
    df_med = st.session_state["median_df"]
    if df_med.empty:
        st.error("Median 후보를 찾지 못했습니다. 파라미터를 조정해보세요.")
    else:
        st.success(f"총 {len(df_med)} 개의 median 후보를 찾았습니다.")

        st.subheader("Median candidates (상위)")
        st.dataframe(df_med.head(50))

        top_n = st.slider(
            "이미지로 볼 median 후보 개수",
            min_value=4,
            max_value=32,
            value=12,
            step=4,
            key="median_top_n",
        )
        top_df = df_med.head(top_n)
        smiles_for_image = top_df["smiles"].tolist()
        legends = [f"score:{row.median_score_sum:.2f}" for _, row in top_df.iterrows()]
        img = smiles_list_to_grid_image(
            smiles_for_image, legends=legends, mols_per_row=4
        )
        st.image(img, caption="Median molecule candidates", use_column_width=True)
else:
    st.info("탐색 버튼을 눌러 결과를 확인하세요.")
