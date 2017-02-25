#!/bin/bash
# Perform some static code analysis with pylint
# Usage: ./run_script.sh [<filenames>]
# If no filenames were provided, it selects all *.py files.

FILES="$@"
if [ -z "$FILES" ];
then
    FILES=$(ls -1 *.py | grep -v '^ui_' | grep -v '^resources_')
fi;

for fname in $FILES;
do
	while true;
	do
                # do code analysis
		(echo "Pylint results for '$fname':";
                 eval "pylint --rcfile=.pylintrc $fname"
                 ret_code=$?) | less -S
		read -p "press ctrl-c for abort, enter anything for retry: "
		if [ -z "$REPLY" ];then
			break;
		fi;
	done
done

exit $ret_code
