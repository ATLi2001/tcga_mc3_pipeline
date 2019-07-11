#!/bin/bash
echo Your container args are: "$@"
echo $OUTPUT_DIR

#if random seed parameter is less than 0, don't pass it on to multievolve
if [ $r -lt 0 ]
then
    python /tmp/phylowgs/multievolve.py \
    -n $n \
    -I $I \
    --ssms $SSM \
    --cnvs $CNV \
    -B $B \
    -s $s \
    -i $i
else
    python /tmp/phylowgs/multievolve.py \
    -n $n \
    -r $r \
    -I $I \
    --ssms $SSM \
    --cnvs $CNV \
    -B $B \
    -s $s \
    -i $i
fi

cp -r ../workingdir/chains/* "$OUTPUT_DIR"
