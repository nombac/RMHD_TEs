# --- system dependent variables ---
ifeq (${SYSTYPE}, sx)
CPUTYPE = vec
FC      = sxmpif90
LD      = sxmpif90
#FFLAGS += -C debug
FFLAGS += -C hopt
FFLAGS += -Wf"-L fmtlist"
FFLAGS += -Wf"-pvctl noloopchg"
#FFLAGS += -ftrace
LDFLAGS =
#LDFLAGS += -ftrace
QUEUE = P8
PRCS = 8
MEMORY = 5gb
FTRACE = NO
SCRATCH = /disk/sx8soto1/shirose
endif

ifeq (${SYSTYPE}, mstsx8rf)
CPUTYPE = vec
FC      = sxmpif90
LD      = sxmpif90
#FFLAGS += -C debug
FFLAGS += -C hopt
FFLAGS += -Wf"-L fmtlist"
FFLAGS += -Wf"-pvctl noloopchg"
#FFLAGS += -ftrace
LDFLAGS =
#LDFLAGS += -ftrace
QUEUE = vq16
PRCS = 8
MEMORY = 5gb
FTRACE = NO
SCRATCH = /sxwork/G10103/shirose
endif

ifeq (${SYSTYPE}, mac)
CPUTYPE = sca
FC      = mpif90
LD      = mpif90
#FC      = g95
#LD      = g95
FFLAGS  = -O0 -Wall
#FFLAGS += -fbounds-check
#FFLAGS += -ftrace=full
LDFLAGS =
QUEUE =
PRCS =
#SCRATCH = /Volumes/IODATA1/zeus_rad/data.rad
SCRATCH = /Volumes/Data/data
endif

ifeq (${SYSTYPE}, vpp)
CPUTYPE = vec
FC      = mpifrt
LD      = mpifrt
FFLAGS  += -O4 -Am -X9
#FFLAGS += -Wv,-m2 -Ps
LDFLAGS =
QUEUE = 
PRCS = 8
SCRATCH = /large4/hirosesg
WP = -Wp,
FFLAGS += ${WP}-DVPP
endif

ifeq (${SYSTYPE}, altix)
FC      = ifort
LD      = ifort
FFLAGS  = -O3 -mp -convert big_endian
LDFLAGS = -lmpi
QUEUE = P1
PRCS =
SCRATCH = /diskwk/store/shirose
endif

ifeq (${SYSTYPE}, mstatx00)
CPUTYPE = sca
FC      = ifort
LD      = ifort
FFLAGS  = -O3 -mp -convert big_endian
LDFLAGS = -lmpi -limf 
QUEUE = sq511
PRCS = 64
SCRATCH = /atwork/G10103/shirose
endif

ifeq (${SYSTYPE}, teragrid)
CPUTYPE = sca
FC      = mpif90
LD      = mpif90
FFLAGS  = -O2 -convert big_endian
LDFLAGS =
QUEUE = debug
WALLTIME = 00:30:00
PRCS = 4
NODES = 2
PPN = 2
RESOURCE = fastcpu
SCRATCH = ${TG_CLUSTER_PFS}
endif
