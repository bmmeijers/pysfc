#!/bin/bash

# find . -name '*.py' | entr ./conttest.sh

echo '* Picked up change (waiting 1 sec)'
sleep 1
echo '* Running test suite'

#python -m unittest discover -s . -p 'test_*.py'
nosetests3 ./tests/test_*.py -v --with-coverage --cover-html
