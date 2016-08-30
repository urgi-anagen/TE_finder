#!/usr/bin/env bash 

if [ $# -lt 3 ]
	then
		echo "usage: $0 <old_string> <new_string> <file(s) or regex>"
	else
	    for i in $*
		do
		    if [ -f $i ]
			then
			    name=`echo $i | sed s/$1/$2/`
			    echo "change: " $i " to " $name
			    mv $i $name
		    fi
		done
fi

