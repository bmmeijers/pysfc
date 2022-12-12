# nD-PC

## getting up and running

To get development started, clone the repository and switch to the dev branch:

```
git clone git@github.com:bmmeijers/pysfc.git
git checkout dev
```

Set up a python environment with all dependencies installed

```
python3 -m venv .env
source .env/bin/activate
python3 -m pip install -e .[test]
```

Now use make to run the test suite.
After the test suite has run the folder `htmlcov` should contain the coverage results (see `index.html` in that folder).

```
make test+htmlcov
```