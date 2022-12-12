test+cov:
	pytest --cov-report term --cov=pysfc tests/

test+htmlcov:
	pytest --cov-report html --cov=pysfc tests/