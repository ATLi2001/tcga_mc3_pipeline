#!/bin/bash
echo Your container args are: "$@"
echo $INPUT
OUTPUT_DIR=$(dirname $OUTPUT_DIR)
echo $OUTPUT_DIR

python /tmp/phylowgs/tree_reader_multi_latest.py \
-t $INPUT \
-s $SSM \
-O $OUTPUT_DIR \
-k $k
