#!/bin/bash

# This converts a string from the Gaussview 5 Atom Selection Editor
# into the python list format for use in TheoDORE analyses

if [[ $# -eq 0 ]] ; then
    echo 'This script converts strings from the Gaussview 5 Atom Selection Editor'
    echo "into the python list format for use in TheoDORE 'dens_ana' files"
    echo 'usage:'
    echo '      GEN_FRAG_LIST.sh string'
    exit 0
fi

# Python list format = [N1, N2, N3] etc.
fragstr="["
# split incoming string by commas
for part in $(echo $1 | sed "s/,/ /g"); do
  # ranges are separated by '-'
  BEGR=$( echo $part | sed 's/\-.*//' )
  ENDR=$( echo $part | cut -d "-" -f2 )
  # sequence with same start and end only prints once
  for i in $( seq ${BEGR} ${ENDR} ); do
    thispart=$( printf "%d, " $i )
    fragstr="${fragstr}${thispart}"
  done
done
# trim the final ', ' and add closing ]
fragstr="${fragstr%, }]"

echo $fragstr
