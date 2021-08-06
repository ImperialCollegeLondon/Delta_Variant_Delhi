#!/bin/bash

PLACE="Delhi"
ITER="2000"
WARM="1000"
JOB="base"
CHAIN="4"
SEED="12121"
DATE="2021-02-14"
UR="0.5"
ifrSD="0.02"
trunc_conv="51"
dir="R"
ifrMean="0.25"
WPAR="310"
RdiffPar1="5"
RdiffPar2="5"
crossPar1="2"
crossPar2="1"

Rscript $dir/base.R --vdate=$DATE \
                    --areaname=$PLACE \
                    --iter=$ITER \
                    --warmup=$WARM \
                    --seed=$SEED \
                    --jname=$JOB \
                    --wpar=$WPAR \
                    --nchains=$CHAIN \
                    --ur=$UR \
                    --ifrSD=$ifrSD \
                    --ifrMean=$ifrMean \
                    --trunc=$trunc_conv \
		    --RdiffPar1=$RdiffPar1 \
		    --RdiffPar2=$RdiffPar2 \
		    --crossPar1=$crossPar1 \
		    --crossPar2=$crossPar2 &
