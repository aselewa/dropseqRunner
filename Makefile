ENV_VARIABLES ?= PYTHONPATH=${PWD}:$PYTHONPATH

run_test_workflow:
	env ${ENV_VARIABLES} python tests/test_workflow.py
	env ${ENV_VARIABLES} python tests/test_reproducibility.py
