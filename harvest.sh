#!/usr/bin/env bash

for dir in res/*/; do
	base_dir=`basename $dir`
	if test -f $dir/annot_merged_filtered.vcf; then
		#echo $base_dir exists >/dev/null
		group_name=`echo $base_dir | sed 's/_.*$//'`
		mkdir -p /demo-mount/m2_results/$group_name
		mv $dir /demo-mount/m2_results/$group_name
	else
		echo $base_dir > /dev/null
	fi
done


