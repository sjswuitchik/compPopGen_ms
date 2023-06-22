#!/bin/awk -f

BEGIN { 
        FS = "[\"\t ]"; OFS="\t" 
}

$(NF - 1) 
{
        print $1, $2, $3, $(NF-1)
}
