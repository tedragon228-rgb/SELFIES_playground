import streamlit as st

st.set_page_config(
    page_title="SELFIES Molecule Playground",
    layout="wide",
)

st.title("SELFIES 기반 분자 탐색 Playground")

st.markdown(
    """
이 앱은 SELFIES 표현을 이용해 분자 구조를 변형·탐색하는 웹 인터페이스입니다.

왼쪽 상단의 페이지 탭을 이용해 다음 기능을 사용할 수 있습니다.

1. **Local Chemical Space** – Seed 분자 주변의 변이 분자 생성
2. **Greedy Chemical Path** – 시작/목표 분자 사이를 SELFIES 변이로 잇는 경로 생성
3. **Median Molecule Search** – 여러 참조 분자에 동시에 유사한 median-like 분자 탐색

필요 패키지: `streamlit`, `streamlit-ketcher`, `selfies`, `rdkit-pypi`, `numpy`, `pandas`, `pillow`

터미널에서 `streamlit run app.py` 로 실행한 뒤 사이드바의 각 페이지를 선택하세요.
"""
)
