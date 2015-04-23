#DEFAULT_SHELL := $(shell which sh)
#DEFAULT_CPP := $(shell which cpp)
#DEFAULT_FC := $(shell which gfortran)

DEFAULT_SHELL := sh
DEFAULT_CPP := cpp
DEFAULT_FC := mpifort

DEFAULT_HYPRE := /home-2/ncsa01/hyprelib_pgi153

DEFAULT_CPPFLAGS += -P -C -traditional

DEFAULT_FFLAGS := -O1
DEFAULT_FFLAGS_HOT := -fast
DEFAULT_FFLAGS_ACC := -acc -Minfo=acc -fast
DEFAULT_FFLAGS_HYPRE := -I$(DEFAULT_HYPRE)/include

DEFAULT_LDFLAGS := -acc -fast
DEFAULT_LIBS := -L$(DEFAULT_HYPRE)/lib -lHYPRE

DEFAULT_BIN_NAME := nhwave



CREATE_CONFIG ?=

ifneq ($(strip $(CREATE_CONFIG)),)



all:
	@mkdir -p config
	@(                                                   \
		echo                                        ;\
		echo 'SHELL := $$(DEFAULT_SHELL)'           ;\
		echo 'CPP := $$(DEFAULT_CPP)'               ;\
		echo 'FC := $$(DEFAULT_FC)'                 ;\
		echo                                        ;\
		echo 'HYPRE := $$(DEFAULT_HYPRE)'           ;\
		echo                                        ;\
		echo 'CPPFLAGS := $$(DEFAULT_CPPFLAGS)'     ;\
		echo                                        ;\
		echo 'FFLAGS := $$(DEFAULT_FFLAGS)'         ;\
		echo 'FFLAGS_HOT := $$(DEFAULT_FFLAGS_HOT)' ;\
		echo 'FFLAGS_ACC := $$(DEFAULT_FFLAGS_ACC)' ;\
		echo 'FFLAGS_HYPRE := $$(DEFAULT_FFLAGS_HYPRE)' ;\
		echo                                        ;\
		echo 'LDFLAGS := $$(DEFAULT_LDFLAGS)'       ;\
		echo 'LIBS := $$(DEFAULT_LIBS)'             ;\
		echo                                        ;\
		echo 'BIN_NAME := $$(DEFAULT_BIN_NAME)'     ;\
		echo                                        ;\
		echo 'OPTIONS += -DDOUBLE_PRECISION'        ;\
		echo 'OPTIONS += -DPARALLEL'                ;\
		echo 'OPTIONS += #-DLANDSLIDE'               ;\
		echo 'OPTIONS += #-DSALINITY'                ;\
		echo 'OPTIONS += #-DTEMPERATURE'             ;\
		echo 'OPTIONS += #-DBUBBLE'                  ;\
		echo 'OPTIONS += #-DSEDIMENT'                ;\
		echo 'OPTIONS += #-DVEGETATION'              ;\
		echo 'OPTIONS += #-DINTEL'                   ;\
		echo 'OPTIONS += #-DBALANCE2D'               ;\
		echo 'OPTIONS += #-DTWOLAYERSLIDE'           ;\
		echo 'OPTIONS += #-DCORALREEF'               ;\
		echo                                        ;\
	) > config/$(CREATE_CONFIG).mk

else

CONFIG ?= default
BUILD_DIR := build/$(CONFIG)

include config/$(CONFIG).mk





SRC_MODS := mod_global.F mod_util.F mod_timing.F
SRC_HOTSPOTS += hotspot_generate_coef_rhs.F
SRC_HOTSPOTS += hotspot_hypre_pres_solver.F
SRC_HOTSPOTS += hotspot_delxyzFun.F
SRC_ORIG := $(SRC_MODS) $(SRC_HOTSPOTS) nhwave.F

SRC := $(shell                                        \
         echo $(SRC_ORIG:.F=.f90) $(wildcard *.f90) | \
	 sed 's/ /\n/g'                             | \
         sed 's$$\(.*\)$$$(BUILD_DIR)/\1$$g'        | \
         sort -u                                      \
         )

OBJ_MODS := $(shell                                   \
		 echo $(SRC_MODS:.F=.o)              |\
                 sed 's/ /\n/g'                      |\
                 sed 's$$\(.*\)$$$(BUILD_DIR)/\1$$g' |\
                 sort -u                              \
             )

OBJ := $(SRC:.f90=.o)

BIN := $(BUILD_DIR)/$(BIN_NAME)

.INTERMEDIATE: $(SRC)






all: $(BIN)



$(BIN): $(OBJ_MODS) $(OBJ)
	$(FC) $(LDFLAGS) $^ -o $@ $(LIBS)

$(BUILD_DIR)/%.f90: %.F $(BUILD_DIR)
	$(CPP) $(CPPFLAGS) $(OPTIONS) $< -o $@

$(BUILD_DIR)/%.f90: %.f90 $(BUILD_DIR)
	cp $< $@

$(BUILD_DIR)/%.o: $(BUILD_DIR)/%.f90 $(BUILD_DIR)
	$(FC) -c $(FFLAGS) $< -o $@

target := $(BUILD_DIR)/hotspot_hypre_pres_solver
$(target).o: $(target).f90 $(BUILD_DIR)
	$(FC) -c $(FFLAGS) $(FFLAGS_HYPRE) $< -o $@

target := $(BUILD_DIR)/hotspot_delxyzFun
$(target).o: $(target).f90 $(BUILD_DIR)
	$(FC) -c $(FFLAGS_HOT) $< -o $@

target := $(BUILD_DIR)/hotspot_generate_coef_rhs
$(target).o: $(target).f90 $(BUILD_DIR)
	$(FC) -c $(FFLAGS_ACC) $< -o $@





$(BUILD_DIR):
	mkdir -p $@

clean:
	rm -rf $(BUILD_DIR)

clean-all:
	rm -rf build/*/

run:
	(                                            \
	cd ../examples/submerged_bar                ;\
	time mpirun -np 1 ../../source/$(BIN)       ;\
	cat reference_result.txt                     \
	)


test:
	@mkdir -p logs
	@(                           \
		set -e                  ;\
		source env/$(CONFIG).sh ;\
		make CONFIG=$(CONFIG)   ;\
		source run/$(CONFIG).sh  \
	) 2>&1 | tee logs/$(CONFIG).log

help:
	@(                                \
	echo 'make CREATE_CONFIG=xxx'    ;\
	echo 'make CONFIG=xxx'           ;\
	echo 'make <VAR_OVERRIDE>=xxx'   ;\
	)

help2:
	@echo $(FC)

	@echo $(SRC)
	@echo $(OBJ_MODS)
	@echo $(OBJ)
	@echo $(BIN)

endif

.SUFFIXES: .o .f90 .F
.PHONY: all clean clean-config test

