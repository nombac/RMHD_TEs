
SHELL = /bin/sh

# --- problem
#PROBLEM = fuori
PROBLEM = ppd
#PROBLEM = raddom
#PROBLEM = iso
#PROBLEM = okuzumi


# --- system dependent variables ---
include config.mk


# --- problem dependent variables ---
include config.mk.$(PROBLEM)


# --- source files ---
SRCS = bval_module.F90 \
       cfl_module.F90 \
       check_module.F90 \
       constant_module.F90 \
       convert_module.F90 \
       diff_module$(CPUTYPE).F90 \
       digit_module.F90 \
       dump_module.F90 \
       dust_module.F90 \
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
       param_module.F90 \
       radiation_module.F90 \
       root_module.F90 \
       shear_module.F90 \
       strtoi_module.F90 \
       tracer_module.F90 \
       viscous_module.F90


# --- object files ---
OBJS = $(SRCS:.F90=.o)
MODS = $(SRCS:.F90=.mod)


# --- execution file ---
EXEC = a.out.$(SYSTYPE)


# --- misc ---
TARG = dummy
EXEDIR = ../data/$(TARG)
JOBSH = Job_$(SYSTYPE).sh


# --- targets ---
all: $(EXEC)
	-mkdir $(EXEDIR)
	mv $(EXEC) $(EXEDIR)
	-mkdir $(EXEDIR)/src
	-mkdir $(EXEDIR)/data
	-mkdir $(EXEDIR)/h
	-mkdir $(EXEDIR)/hist
	-mkdir $(EXEDIR)/ps
	-mkdir $(EXEDIR)/tmp
	sed -e 's/Targ/$(TARG)/g' -e 's/Queue/$(QUEUE)/g' \
            -e 's/Prcs/$(PRCS)/g' -e 's/Nodes/$(NODES)/g' \
            -e 's/Ppn/$(PPN)/g' -e 's/Resource/$(RESOURCE)/g'\
            -e 's/Walltime/$(WALLTIME)/g'\
            -e 's/Memory/$(MEMORY)/g'\
            -e 's/Ftrace/$(FTRACE)/g'\
            $(JOBSH) > $(EXEDIR)/$(TARG).sh
	cp -pR . $(EXEDIR)/src/
ifeq ($(SYSTYPE), ha8000)
else
ifeq ($(SYSTYPE), xt)
else
	-mkdir $(SCRATCH)/$(TARG)
	-mkdir $(SCRATCH)/$(TARG).h
	-mkdir $(SCRATCH)/$(TARG).hist
	-mkdir $(SCRATCH)/$(TARG).tmp
	-\rm -r $(EXEDIR)/data
	-\rm -r $(EXEDIR)/h
	-\rm -r $(EXEDIR)/hist
	-\rm -r $(EXEDIR)/tmp
	ln -s $(SCRATCH)/$(TARG) $(EXEDIR)/data
	ln -s $(SCRATCH)/$(TARG).h $(EXEDIR)/h
	ln -s $(SCRATCH)/$(TARG).hist $(EXEDIR)/hist
	ln -s $(SCRATCH)/$(TARG).tmp $(EXEDIR)/tmp
endif
endif

ifeq ($(PROBLEM), ppd)
RESIST_OBJ = 
ifeq ($(CHEMEQ), ohmeq)
RESIST_OBJ = ohmeq/ohmeq.o
endif
ifeq ($(CHEMEQ), ohmeq2)
RESIST_OBJ = ohmeq2/ohmeq2.o ohmeq2/xrioniz.o
endif
ifeq ($(CHEMEQ), resist1d)
RESIST_OBJ = ah.r13/adios.o ah.r13/lubksb.o ah.r13/nrutil.o ah.r13/pzextr.o ah.r13/simpr.o ah.r13/chemeq1d.o ah.r13/ludcmp.o ah.r13/odeint.o ah.r13/react.o ah.r13/stifbs.o
endif
else
RESIST_OBJ =
endif

$(EXEC): $(OBJS) $(RESIST_OBJ)
	$(LD) -o $(EXEC) $(OBJS) $(RESIST_OBJ) $(LDFLAGS)

%.o: %.F90
	$(FC) -c $(FFLAGS) -o $@ $<
#	$(FC) -c $(FFLAGS) $<

%.o: %.c
	$(CC) -c $(CFLAGS) -o $@ $<

# --- module dependencies ---
bval_module.o: field_module.o grid_module.o mpi_module.o root_module.o param_module.o
cfl_module.o: field_module.o grid_module.o mpi_module.o param_module.o root_module.o
check_module.o: mpi_module.o
constant_module.o : 
convert_module.o: bval_module.o field_module.o grid_module.o
digit_module.o :
diff_module$(CPUTYPE).o: bval_module.o digit_module.o field_module.o grid_module.o mpi_module.o param_module.o
dump_module.o: constant_module.o dust_module.o field_module.o flux_module.o gather_module.o grid_module.o mpi_module.o param_module.o root_module.o strtoi_module.o tracer_module.o
dust_module.o : bval_module.o field_module.o grid_module.o interp_module.o param_module.o root_module.o
field_module.o : constant_module.o gather_module.o grid_module.o param_module.o root_module.o
floor_module.o: bval_module.o field_module.o flux_module.o grid_module.o param_module.o root_module.o
flux_module.o : grid_module.o
gather_module.o: grid_module.o mpi_module.o
grid_module.o: mpi_module.o 
hist_module.o: constant_module.o digit_module.o field_module.o flux_module.o grid_module.o mpi_module.o param_module.o root_module.o radiation_module.o strtoi_module.o
hydro_module.o: bval_module.o convert_module.o field_module.o flux_module.o floor_module.o grid_module.o interp_module.o param_module.o root_module.o
interp_module.o: grid_module.o root_module.o
initial_module.o : constant_module.o dump_module.o field_module.o flux_module.o gather_module.o grid_module.o hist_module.o magnetic_module.o param_module.o radiation_module.o root_module.o tracer_module.o
isothermal_module.o : constant_module.o bval_module.o param_module.o field_module.o grid_module.o
magnetic_module.o: bval_module.o convert_module.o field_module.o flux_module.o grid_module.o root_module.o
main.o: bval_module.o cfl_module.o constant_module.o dump_module.o floor_module.o hist_module.o hydro_module.o isothermal_module.o magnetic_module.o mpi_module.o root_module.o viscous_module.o shear_module.o tracer_module.o check_module.o
mpi_module.o: digit_module.o
param_module.o: constant_module.o grid_module.o
radiation_module.o: bval_module.o constant_module.o diff_module$(CPUTYPE).o field_module.o flux_module.o grid_module.o interp_module.o mpi_module.o param_module.o root_module.o check_module.o
shear_module.o: bval_module.o field_module.o flux_module.o grid_module.o root_module.o 
strtoi_module.o:     
tracer_module.o: bval_module.o field_module.o grid_module.o
viscous_module.o: bval_module.o convert_module.o field_module.o flux_module.o grid_module.o interp_module.o mpi_module.o param_module.o root_module.o

# --- clean your room ---
clean:
	-\rm *.o *.mod a.out.* *~ i.*.L *.lst *.s ah.r13/*.o ohmeq*/*.o ah.r13/*~ ohmeq*/*~
