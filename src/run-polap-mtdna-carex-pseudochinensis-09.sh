#!/usr/bin/env bash

cat s1.fq s2.fq |
	awk 'NR % 4 == 2' | sort | tr NT TN | ropebwt2 -LR | tr NT TN |
	msbwt convert o/msbwt \
		>/dev/null 2>&1
