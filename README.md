# SELFIES Playground

Streamlit 기반의 SELFIES 분자 탐색 앱입니다. 다음과 같이 세 가지 페이지로 구성된 멀티 페이지 앱으로 동작합니다.

1. **Local Chemical Space** – Seed 분자를 기준으로 SELFIES 변이를 적용해 주변 분자 생성
2. **Greedy Chemical Path** – 시작/목표 분자 사이를 greedy하게 이동하는 SELFIES 경로 생성
3. **Median Molecule Search** – 여러 참조 분자에 동시에 유사한 median-like 분자 탐색

## 실행 방법

필요 패키지: `streamlit`, `streamlit-ketcher`, `selfies`, `rdkit-pypi`, `numpy`, `pandas`, `pillow`

```bash
pip install streamlit streamlit-ketcher selfies rdkit-pypi numpy pandas pillow
streamlit run app.py
```

좌측 사이드바의 페이지 목록에서 원하는 기능을 선택한 뒤, 설정을 변경하고 실행하세요. 생성된 결과는 세션 상태에 저장되므로 슬라이더나 숫자 입력값을 조절해도 결과 테이블이 초기화되지 않습니다.
