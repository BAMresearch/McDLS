#!/bin/bash
# modifies system libraries to contain absolute paths for cxfreeze to find them
# to be run as root

LIBPATH="$@"
echo "processing '$LIBPATH'"
for oldname in $(otool -L "$LIBPATH" | grep -v ':' | awk '{print $1}');
do
	if [ -f "$oldname" ]; then
		continue;
	fi;
	echo -n "  $oldname -> "
	pattern="$oldname"
	if [ "${pattern##*.}" == "dylib" ]; then # test file ending
		pattern="${pattern##*/}"         # remove leading paths
	fi;
	if [[ "${pattern:0:1}" != "/" ]]; then   # prepend / if missing
		pattern="/$pattern"
	fi;
	absname="$(find -x / -path "*$pattern" 2> /dev/null | grep -v '/Applications' | head -n 1)"
	echo "$absname"
	install_name_tool -change "$oldname" "$absname" "$LIBPATH"
done;
echo -n "result: "
otool -L "$LIBPATH"

