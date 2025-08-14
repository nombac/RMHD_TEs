SHELL = /bin/sh

# --- system dependent variables ---
include config.mk
ifeq ($(SYSTYPE), macosx_intel)
FFLAGS += -DDEBUGD
endif
ifeq ($(SYSTYPE), macosx)
FFLAGS += -DDEBUGD
endif

# --- common variables ---
#FFLAGS += -DSPEEDTEST
FFLAGS += -DMPI

# --- problem
#PROBLEM = pp
#PROBLEM = unstratified_rad
#PROBLEM = unstratified_wd
#PROBLEM = ws
#PROBLEM = test_selfgravshear
#PROBLEM = test_hdshearwave
#PROBLEM = test_stonemri
#PROBLEM = test_loop
#PROBLEM = self
#PROBLEM = sano
#PROBLEM = ws2
PROBLEM = self_rad
#PROBLEM = self_stratified
#PROBLEM = is
#PROBLEM = unstratified_iso
#PROBLEM = unstratified_ad
#PROBLEM = raddom
#PROBLEM = iso
#PROBLEM = okuzumi
#PROBLEM = shock
#PROBLEM = dummy
# --- problem dependent variables ---
include config.mk.$(PROBLEM)

# --- source files ---
SRCS = advection_module.F90 \
       bval_module.F90 \
       cfl_module.F90 \
       check_module.F90 \
       constant_module.F90 \
       convert_module.F90 \
       diff_module.F90 \
       digit_module.F90 \
       dump_module.F90 \
       field_module.F90 \
       floor_module.F90 \
       flux_module.F90 \
       gather_module.F90 \
       grid_module.F90 \
       hist_module.F90 \
       hydro_module.F90 \
       initial_module.F90 \
       interp_module.F90 \
       isothermal_module.F90 \
       magnetic_module.F90 \
       main.F90 \
       mpi_module.F90 \
       eos_module.F90 \
       opacity_module.F90 \
       param_module.F90 \
       poisson_fft_module.F90 \
       problem_module.F90 \
       radiation_module.F90 \
       remap_module.F90 \
       resist_module.F90 \
       root_module.F90 \
       selfgrav_module.F90 \
       shear_module.F90 \
       strtoi_module.F90 \
       test_module.F90 \
       viscous_module.F90

ifeq ($(PROBLEM), pp)
SRC_RES = ohmeq2/ohmeq2.c ohmeq2/xrioniz.c
OBJ_RES = $(SRC_RES:.c=.o)
endif

# --- object files ---
OBJS = $(SRCS:.F90=.o) $(OBJ_RES)
MODS = $(SRCS:.F90=.mod)

# --- execution file ---
EXEC = a.out.$(SYSTYPE)

# --- input file ---
Z3DINPUT = z3dinput

# --- execution directory ---
EXEDIR = ../data/$(TARG)

# --- target directory
TARG = default

# --- targets ---
all: $(EXEC)
	-mkdir $(EXEDIR)
	-mkdir $(EXEDIR)/init
	-\rm $(EXEDIR)/$(EXEC)
	\cp $(EXEC) $(EXEDIR)/$(EXEC)
	-mkdir $(EXEDIR)/src
	-mkdir $(EXEDIR)/data
	-mkdir $(EXEDIR)/h
	-mkdir $(EXEDIR)/hist
	-mkdir $(EXEDIR)/ps
	-mkdir $(EXEDIR)/tmp
	-mkdir $(EXEDIR)/vtk
	-mkdir $(EXEDIR)/gnuplot
	sed -e 's/Targ/$(TARG)/g'\
	    -e 's/Queue/$(QUEUE)/g' \
            -e 's/Prcs/$(PRCS)/g' \
            -e 's/Nodes/$(NODES)/g' \
            -e 's/Ppn/$(PPN)/g' \
            -e 's/Resource/$(RESOURCE)/g'\
            -e 's/Walltime/$(WALLTIME)/g'\
            -e 's/Memory/$(MEMORY)/g'\
            -e 's/Ftrace/$(FTRACE)/g'\
            -e 's/Systype/$(SYSTYPE)/g'\
            -e 's/Group/$(GROUP)/g'\
            Job_$(SYSTYPE).sh > $(EXEDIR)/$(TARG).sh
	\cp -pR . $(EXEDIR)/src/
ifdef SCRATCH
	-mkdir $(SCRATCH)/$(TARG).data
	-mkdir $(SCRATCH)/$(TARG).hist
	-\rm -r $(EXEDIR)/data
	-\rm -r $(EXEDIR)/hist
	ln -s $(SCRATCH)/$(TARG).data $(EXEDIR)/data
	ln -s $(SCRATCH)/$(TARG).hist $(EXEDIR)/hist
endif
	echo 0 > $(EXEDIR)/abort.data
	echo 0 > $(EXEDIR)/restart.data
	echo 0 > $(EXEDIR)/status.data

$(EXEC): $(OBJS)
	$(FC) -o $(EXEC) $(OBJS) $(LDFLAGS)

%.o: %.F90
	$(FC) -c $(FFLAGS) -o $@ $<

%.o: %.c
	$(CC) -c $(CFLAGS) -o $@ $<

# --- module dependencies ---
advection_module.o: bval_module.o convert_module.o field_module.o flux_module.o floor_module.o grid_module.o hydro_module.o interp_module.o param_module.o root_module.o
bval_module.o: field_module.o grid_module.o mpi_module.o root_module.o param_module.o
cfl_module.o: field_module.o flux_module.o grid_module.o mpi_module.o param_module.o root_module.o
check_module.o: mpi_module.o
constant_module.o : 
convert_module.o: bval_module.o field_module.o grid_module.o
digit_module.o :
diff_module.o: bval_module.o constant_module.o digit_module.o field_module.o grid_module.o mpi_module.o eos_module.o param_module.o problem_module.o
diffb_module.o: bval_module.o digit_module.o field_module.o grid_module.o mpi_module.o param_module.o
diffd_module.o: bval_module.o digit_module.o field_module.o grid_module.o mpi_module.o param_module.o
diffg_module.o: bval_module.o digit_module.o field_module.o grid_module.o mpi_module.o param_module.o
dump_module.o: constant_module.o field_module.o flux_module.o gather_module.o grid_module.o mpi_module.o param_module.o root_module.o strtoi_module.o problem_module.o
field_module.o : constant_module.o gather_module.o grid_module.o eos_module.o param_module.o root_module.o
floor_module.o: bval_module.o field_module.o flux_module.o grid_module.o param_module.o root_module.o problem_module.o
flux_module.o : grid_module.o
gather_module.o: grid_module.o mpi_module.o
grid_module.o: mpi_module.o param_module.o
hist_module.o: constant_module.o digit_module.o field_module.o flux_module.o grid_module.o mpi_module.o param_module.o root_module.o radiation_module.o strtoi_module.o shear_module.o
hydro_module.o: bval_module.o convert_module.o field_module.o flux_module.o floor_module.o grid_module.o interp_module.o param_module.o problem_module.o radiation_module.o root_module.o
interp_module.o: grid_module.o root_module.o
initial_module.o : constant_module.o dump_module.o field_module.o flux_module.o gather_module.o grid_module.o hist_module.o magnetic_module.o eos_module.o param_module.o problem_module.o radiation_module.o root_module.o selfgrav_module.o opacity_module.o shear_module.o
isothermal_module.o : constant_module.o bval_module.o param_module.o field_module.o grid_module.o
magnetic_module.o: bval_module.o convert_module.o field_module.o flux_module.o grid_module.o resist_module.o root_module.o
main.o: advection_module.o bval_module.o cfl_module.o constant_module.o dump_module.o floor_module.o hist_module.o hydro_module.o isothermal_module.o magnetic_module.o mpi_module.o root_module.o viscous_module.o shear_module.o selfgrav_module.o test_module.o
mpi_module.o: digit_module.o
eos_module.o: constant_module.o mpi_module.o
opacity_module.o: mpi_module.o
param_module.o: constant_module.o eos_module.o opacity_module.o
poisson_fft_module.o: constant_module.o mpi_module.o
problem_module.o: bval_module.o field_module.o grid_module.o eos_module.o param_module.o
radiation_module.o: bval_module.o constant_module.o diff_module.o dump_module.o field_module.o flux_module.o grid_module.o interp_module.o mpi_module.o param_module.o root_module.o check_module.o opacity_module.o shear_module.o
remap_module.o :
resist_module.o : bval_module.o constant_module.o convert_module.o field_module.o flux_module.o grid_module.o eos_module.o 
root_module.o : param_module.o
selfgrav_module.o: bval_module.o field_module.o flux_module.o grid_module.o mpi_module.o root_module.o gather_module.o poisson_fft_module.o remap_module.o
shear_module.o: bval_module.o field_module.o flux_module.o grid_module.o root_module.o 
strtoi_module.o:
test_module.o: field_module.o grid_module.o
viscous_module.o: bval_module.o convert_module.o field_module.o flux_module.o grid_module.o interp_module.o mpi_module.o param_module.o root_module.o

# --- clean your room ---
clean:
	-\rm *.o *.mod a.out.* *~ i.*.L *.lst *.s ah.r13/*.o ohmeq*/*.o ah.r13/*~ ohmeq*/*~
