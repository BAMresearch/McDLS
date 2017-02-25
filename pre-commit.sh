#!/bin/sh
# some notes to be removed when done:
# FIXME: too aggressive atm, code has to be fixed first [TODO]
# TODO: quick&simple regression tests can be added here
#       The full reg.test (atm) needs too much time to run it at every commit
#       (due to the nature of the mc-algorithm)
#
# This pre-commit hook verifies the code base regarding naming conventions
# and formatting. Pylint performs also some static code analysis.
# http://www.pylint.org/
#
# Add this script as pre-commit hook to your local git repository:
# 
# ln -s ../../pre-commit.sh .git/hooks/pre-commit
# 
# Skip the pre-commit hook sometimes:
#
# git commit --no-verify
#
# Hints and recommendation from:
# http://codeinthehole.com/writing/tips-for-using-a-git-pre-commit-hook/

FILES="McSAS.py"

# put local not staged changes aside
git stash -q --keep-index

eval "./run_pylint.sh"
ret_code=$?

# restore local changes not staged yet
git stash pop -q

echo "pre-commit hook exits with code '$ret_code'."
exit $ret_code
