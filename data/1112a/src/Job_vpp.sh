#!/bin/sh
# @$-g ish26b
# @$-q vppPrcsQueue
# @$-mi
VPP_MBX_SIZE=335544320; export VPP_MBX_SIZE
cd ${HOME}/zeus_rad/data/Targ
timex ./a.out.vpp -Wl,-g3000
