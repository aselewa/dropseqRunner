ENV_VARIABLES ?= PYTHONPATH=${PWD}:${PYTHONPATH}

help:
	cat Makefile
	
run_test_workflow:
	env ${ENV_VARIABLES} python tests/test_workflow.py
	env ${ENV_VARIABLES} python tests/test_reproducibility.py
