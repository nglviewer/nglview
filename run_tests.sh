echo "Ensure to install requirements in: "
echo "pip install -r pip-requirements-test.txt"
echo "sh conda-requirements-test.sh"
pytest -vs nglview/tests --cov=nglview --cov-report=html
