#!/bin/bash


if [[ $# -ne 1 ]]; then
    echo "invalid number of params!"
    echo "usage:"
    echo "    $ ./fix_permissions.sh path/to/a/folder"
    exit 1
fi

if [ ! -d $1 ]; then
    echo "error: ${1} is not a folder"
    exit 1
fi

chmod 775 $1
chgrp -R CM $1
find $1 -type f -exec chmod 664 -- {} +
find $1 -type d -exec chmod 775 -- {} +
