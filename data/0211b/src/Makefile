SHELL = /bin/sh

PROBLEM = raddom

# --- system dependent variables ---
include config.mk

# --- preprocessor macros
FFLAGS += $(WP)-DCOMPTON=compton$(CPUTYPE)

FFLAGS += $(WP)-DMPI

FFLAGS += $(WP)-DSHEAR

FFLAGS += $(WP)-DPARTETOT
FFLAGS += $(WP)-DPARTETOTK

FFLAGS += $(WP)-DFLOORV3
FFLAGS += $(WP)-DFLOORD
FFLAGS += $(WP)-DFLOORE

FFLAGS += $(WP)-DDAMPING

# --- source files ---
SRCS = bval_module.F90 \
       cfl_module.F90 \
       check_module.F90 \
       constant_module.F90 \
       convert_module.F90 \
       diff_module$(CPUTYPE).F90 \
       digit_module.F90 \
       dump_module.F90 \
       field_module.F90 \
       floor_module.F90 \
       flux_module.F90 \
       gather_module.F90 \
       grid_module.F90 \
       hist_module.F90 \
       initial_module.F90 \
       interp_module.F90 \
       magnetic_module.F90 \
       main.F90 \
       mpi_module.F90 \
       param_module.F90 \
       radiation_module.F90 \
       root_module.F90 \
       shear_module.F90 \
       source_module.F90 \
       strtoi_module.F90 \
       transport_module.F90 \
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
	cp -p * $(EXEDIR)/src/
ifneq ($(SYSTYPE), xt)
	-mkdir $(SCRATCH)/$(TARG)
	-mkdir $(SCRATCH)/$(TARG).h
	-mkdir $(SCRATCH)/$(TARG).tmp
	-\rm -r $(EXEDIR)/data
	-\rm -r $(EXEDIR)/h
	-\rm -r $(EXEDIR)/tmp
	ln -s $(SCRATCH)/$(TARG) $(EXEDIR)/data
	ln -s $(SCRATCH)/$(TARG).h $(EXEDIR)/h
	ln -s $(SCRATCH)/$(TARG).tmp $(EXEDIR)/tmp
endif

$(EXEC): $(OBJS)
	$(LD) -o $(EXEC) $(OBJS) $(LDFLAGS)

%.o: %.F90
	$(FC) -c $(FFLAGS) -o $@ $<


# --- module dependencies ---
bval_module.o: field_module.o grid_module.o mpi_module.o root_module.o param_module.o
cfl_module.o: field_module.o grid_module.o mpi_module.o param_module.o root_module.o
check_module.o: mpi_module.o
constant_module.o : 
convert_module.o: field_module.o floor_module.o grid_module.o
digit_module.o :
diff_module$(CPUTYPE).o: bval_module.o digit_module.o field_module.o grid_module.o mpi_module.o param_module.o
dump_module.o: constant_module.o field_module.o flux_module.o gather_module.o grid_module.o mpi_module.o param_module.o root_module.o strtoi_module.o
field_module.o : constant_module.o gather_module.o grid_module.o param_module.o root_module.o
floor_module.o: bval_module.o field_module.o flux_module.o grid_module.o param_module.o root_module.o
flux_module.o : grid_module.o
gather_module.o: grid_module.o mpi_module.o
grid_module.o: mpi_module.o 
hist_module.o: constant_module.o digit_module.o field_module.o flux_module.o grid_module.o mpi_module.o param_module.o radiation_module.o root_module.o strtoi_module.o
interp_module.o: grid_module.o root_module.o
initial_module.o : dump_module.o field_module.o flux_module.o gather_module.o grid_module.o hist_module.o param_module.o radiation_module.o root_module.o
magnetic_module.o: bval_module.o convert_module.o field_module.o flux_module.o grid_module.o root_module.o
main.o: bval_module.o cfl_module.o constant_module.o dump_module.o floor_module.o hist_module.o magnetic_module.o mpi_module.o root_module.o source_module.o transport_module.o viscous_module.o check_module.o
mpi_module.o: digit_module.o
param_module.o: constant_module.o grid_module.o
radiation_module.o: bval_module.o constant_module.o diff_module$(CPUTYPE).o field_module.o flux_module.o grid_module.o interp_module.o mpi_module.o param_module.o root_module.o check_module.o
shear_module.o: bval_module.o field_module.o flux_module.o grid_module.o root_module.o 
source_module.o: bval_module.o convert_module.o field_module.o flux_module.o grid_module.o interp_module.o param_module.o radiation_module.o root_module.o shear_module.o
strtoi_module.o:     
transport_module.o: bval_module.o convert_module.o field_module.o floor_module.o flux_module.o grid_module.o interp_module.o radiation_module.o root_module.o
viscous_module.o: bval_module.o convert_module.o field_module.o flux_module.o grid_module.o interp_module.o mpi_module.o param_module.o root_module.o


# --- clean your room ---
clean:
	-\rm *.o *.mod a.out.* *~ i.*.L *.lst
