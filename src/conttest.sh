#!/bin/bash

# find . -name '*.py' | entr ./conttest.sh

echo '* Picked up change (waiting for 2 secs)'
sleep 2
echo '* Running test suite'

#python -m unittest discover -s . -p 'test_*.py'
nosetests-2.7 ./tests/test_*.py -v --with-coverage --cover-html
