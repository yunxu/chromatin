#!/bin/bash - 

# Remove head line of chr16.fa
# head -n 10001 chr16.fa

paste - - < alpha.fa | tr -d "\t" > 1.fa
# or
sed -e '1,1d' alpha.fa | paste - - | tr -d "\t" > 2.fa

# paste -s alpha.fa | tr -d "\t" | grep -i --color aagctt | vim -
# :%s/aagctt//gn
# 75 match
export GREP_COLORS="mt=01;31:ln=34"
gre p--color -i -n -T aagctt 1.fa 



#Attributes - attr
#    0 - Reset All Attributes (return to normal mode)
#    1 - Bright (Usually turns on BOLD)
#    2 - Dim
#    3 - Underline
#    5 - Blink
#    7 - Reverse
#    8 - Hidden
#
#Foreground color - fg
#    30 - Black
#    31 - Red
#    32 - Green
#    33 - Yellow
#    34 - Blue
#    35 - Magenta
#    36 - Cyan
#    37 - White
#
#Background color - bg
#    40 - Black
#    41 - Red
#    42 - Green
#    43 - Yellow
#    44 - Blue
#    45 - Magenta
#    46 - Cyan
#    47 - White

