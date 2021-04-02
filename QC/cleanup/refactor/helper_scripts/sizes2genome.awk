#!/bin/awk -f 

BEGIN {
        FS = OFS = "\t"
}
{
print $1, 0, $2, $1
}
