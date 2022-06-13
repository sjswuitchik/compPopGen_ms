#!/bin/awk -f

BEGIN {
        FS = OFS = "\t"
} 

{
if ($3 == "CDS")
        print $0
}
