test+cov:
	pytest --cov-report term --cov=pysfc tests/

test+covhtml:
	pytest --cov-report html --cov=pysfc tests/