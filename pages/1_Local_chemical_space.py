import streamlit as st
from streamlit_ketcher import st_ketcher

from app_utils import (
    generate_local_chemical_space,
    get_selfies_alphabet,
    smiles_list_to_grid_image,
    smiles_to_mol,
)

st.title("1) Local chemical space")

st.markdown(
    """
Seed 분자를 중심으로 SELFIES 변이를 적용해 주변 분자를 생성합니다.
좌측 패널에서 Seed를 직접 입력하거나 그려서 추가한 뒤, 생성 버튼을 눌러 탐색하세요.
"""
)

random_seed = st.sidebar.number_input("Random seed", min_value=0, value=0, step=1)
allow_charged = st.sidebar.checkbox("생성된 분자에서 이온(+ / -) 허용", value=False)

seed_input_mode = st.radio(
    "Seed 분자 입력 방식 (alphabet 구성용)",
    ["SMILES 텍스트", "구조 그리기"],
    horizontal=True,
)

if "local_seed_smiles_drawn" not in st.session_state:
    st.session_state["local_seed_smiles_drawn"] = []
if "local_space_df" not in st.session_state:
    st.session_state["local_space_df"] = None

seed_smiles = []

if seed_input_mode == "SMILES 텍스트":
    seed_smiles_input = st.text_area(
        "Seed SMILES (한 줄에 하나씩)",
        value="CCO\nCCN\nc1ccccc1",
        height=100,
    )
    seed_smiles = [s.strip() for s in seed_smiles_input.splitlines() if s.strip()]
else:
    st.caption("Ketcher에서 구조를 수정 후 Apply(✔)를 누른 뒤 버튼을 사용하세요.")

    smiles_seed_drawn = st_ketcher(key="local_seed_ketcher")

    col_seed1, col_seed2 = st.columns(2)
    with col_seed1:
        if st.button("현재 구조를 Seed 목록에 추가"):
            if smiles_seed_drawn and smiles_seed_drawn.strip():
                smi = smiles_seed_drawn.strip()
                if smiles_to_mol(smi) is None:
                    st.error("그려진 구조(SMILES)를 RDKit에서 해석하지 못했습니다.")
                else:
                    st.session_state["local_seed_smiles_drawn"].append(smi)
            else:
                st.warning("먼저 구조를 그리고 Ketcher 안에서 Apply(✔)를 누른 뒤 다시 시도해 주세요.")
    with col_seed2:
        if st.button("Seed 목록 비우기"):
            st.session_state["local_seed_smiles_drawn"] = []

    if st.session_state["local_seed_smiles_drawn"]:
        st.write("현재 Seed SMILES 목록:")
        for i, smi in enumerate(st.session_state["local_seed_smiles_drawn"], start=1):
            st.write(f"{i}. `{smi}`")
    else:
        st.info("아직 추가된 Seed가 없습니다. 기본 SELFIES alphabet만 사용됩니다.")

    seed_smiles = list(st.session_state["local_seed_smiles_drawn"])

alphabet = get_selfies_alphabet(seed_smiles)
st.write(f"SELFIES alphabet size: **{len(alphabet)}** (seed 기반 + 기본 alphabet)")

col1, col2 = st.columns(2)
with col1:
    n_mutations_per_seed = st.number_input(
        "각 seed 당 생성 시도 수",
        min_value=1,
        max_value=5000,
        value=256,
        step=1,
        key="local_mutations_per_seed",
    )
    max_total = st.number_input(
        "전체 최대 생성 분자 수",
        min_value=10,
        max_value=10000,
        value=1000,
        step=10,
        key="local_max_total",
    )
with col2:
    min_edits = st.number_input(
        "최소 편집 횟수 (SELFIES token)",
        min_value=1,
        max_value=10,
        value=1,
        step=1,
        key="local_min_edits",
    )
    max_edits = st.number_input(
        "최대 편집 횟수 (SELFIES token)",
        min_value=min_edits,
        max_value=10,
        value=3,
        step=1,
        key="local_max_edits",
    )

if st.button("Local chemical space 생성", type="primary"):
    with st.spinner("SELFIES 공간 탐색 중..."):
        st.session_state["local_space_df"] = generate_local_chemical_space(
            seed_smiles=seed_smiles if seed_smiles else ["CCO"],
            alphabet=alphabet,
            n_mutations_per_seed=int(n_mutations_per_seed),
            min_edits=int(min_edits),
            max_edits=int(max_edits),
            max_total=int(max_total),
            random_seed=int(random_seed),
            allow_charged=allow_charged,
        )

if st.session_state["local_space_df"] is not None:
    df_space = st.session_state["local_space_df"]
    st.success(f"총 {len(df_space)} 개의 유효 분자를 생성했습니다.")
    st.dataframe(df_space)

    if not df_space.empty:
        top_n = st.slider(
            "이미지로 볼 개수",
            min_value=4,
            max_value=64,
            value=16,
            step=4,
            key="local_top_n",
        )
        top_df = df_space.head(top_n)
        smiles_for_image = top_df["mutated_smiles"].tolist()
        legends = [
            f"seed:{row.seed_index}, sim:{row.similarity_to_seed:.2f}"
            for _, row in top_df.iterrows()
        ]
        img = smiles_list_to_grid_image(
            smiles_for_image, legends=legends, mols_per_row=4
        )
        st.image(img, caption="Local chemical space (상위 분자들)", use_column_width=True)
else:
    st.info("생성 버튼을 눌러 결과를 확인하세요.")
