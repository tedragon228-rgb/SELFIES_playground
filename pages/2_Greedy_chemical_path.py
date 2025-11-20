import streamlit as st
from streamlit_ketcher import st_ketcher

from app_utils import (
    get_default_selfies_alphabet,
    greedy_path_between_two,
    smiles_list_to_grid_image,
    smiles_to_mol,
)

st.title("2) Greedy chemical path between two molecules")

st.markdown(
    """
시작/목표 분자를 SMILES로 입력하거나 직접 그려서 SELFIES 변이 기반 경로를 만듭니다.
"""
)

random_seed = st.sidebar.number_input("Random seed", min_value=0, value=0, step=1)
allow_charged = st.sidebar.checkbox("생성된 분자에서 이온(+ / -) 허용", value=False)

alphabet = get_default_selfies_alphabet()
st.write(f"SELFIES alphabet size: **{len(alphabet)}** (기본 alphabet)")

path_input_mode = st.radio(
    "Path 분자 입력 방식",
    ["SMILES 텍스트", "구조 그리기"],
    horizontal=True,
)

if "path_df" not in st.session_state:
    st.session_state["path_df"] = None

start_smi = None
target_smi = None

if path_input_mode == "SMILES 텍스트":
    col1, col2 = st.columns(2)
    with col1:
        start_smi = st.text_input("시작 SMILES", value="CCO")
    with col2:
        target_smi = st.text_input("목표 SMILES", value="c1ccccc1O")
else:
    st.caption("구조 수정 후 Ketcher 안의 Apply(✔)를 눌러야 값이 반영됩니다.")
    colp1, colp2 = st.columns(2)
    with colp1:
        st.write("시작 분자")
        smiles_start_drawn = st_ketcher(key="path_start_ketcher")
    with colp2:
        st.write("목표 분자")
        smiles_target_drawn = st_ketcher(key="path_target_ketcher")

    start_smi = smiles_start_drawn.strip() if smiles_start_drawn else None
    target_smi = smiles_target_drawn.strip() if smiles_target_drawn else None

c3, c4 = st.columns(2)
with c3:
    max_steps = st.number_input(
        "최대 스텝 수",
        min_value=1,
        max_value=100,
        value=20,
        step=1,
        key="path_max_steps",
    )
with c4:
    neighbors_per_step = st.number_input(
        "각 스텝당 후보 생성 개수",
        min_value=8,
        max_value=512,
        value=128,
        step=8,
        key="path_neighbors_per_step",
    )

if st.button("Chemical path 생성", type="primary"):
    if not start_smi or not target_smi:
        st.error(
            "시작/목표 분자를 올바르게 입력 또는 그려 주세요. 그림 모드에서는 Apply(✔) 버튼을 눌러야 합니다."
        )
    elif smiles_to_mol(start_smi) is None or smiles_to_mol(target_smi) is None:
        st.error("시작/목표 SMILES를 RDKit에서 해석하지 못했습니다.")
    else:
        with st.spinner("SELFIES 공간에서 greedy path 탐색 중..."):
            st.session_state["path_df"] = greedy_path_between_two(
                start_smiles=start_smi,
                target_smiles=target_smi,
                alphabet=alphabet,
                max_steps=int(max_steps),
                neighbors_per_step=int(neighbors_per_step),
                random_seed=int(random_seed),
                allow_charged=allow_charged,
            )

if st.session_state["path_df"] is not None:
    df_path = st.session_state["path_df"]
    if df_path.empty:
        st.error("경로를 생성하지 못했습니다. SMILES나 파라미터를 조정해보세요.")
    else:
        st.success(f"경로 길이: {len(df_path)} (step 0 포함)")

        st.subheader("경로 테이블")
        st.dataframe(df_path)

        st.subheader("경로 분자 이미지")
        smiles_for_image = df_path["smiles"].tolist()
        legends = [
            f"step {row.step}, sim:{row.similarity_to_target:.2f}"
            for _, row in df_path.iterrows()
        ]
        img = smiles_list_to_grid_image(smiles_for_image, legends=legends, mols_per_row=4)
        st.image(img, caption="Greedy path", use_column_width=True)
else:
    st.info("생성 버튼을 눌러 경로를 확인하세요.")
