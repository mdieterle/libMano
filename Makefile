S             = src
I             = include
EXI           = include
O             = obj
L             = lib
AI            = -I./$(I) -I$(OSCAR)/$(I) -I$(PLUTOLIB)/$(S)

SRC           = $(wildcard $(S)/TM*.cxx)
INCD          = $(wildcard $(I)/TM*.h)
INC           = $(notdir $(INCD))
OBJD          = $(patsubst $(S)/%.cxx, $(O)/%.o, $(SRC))
OBJ           = $(notdir $(OBJD))
OBJD         += $(O)/G__MLIB.o

OSTYPE       := $(subst -,,$(shell uname))

ROOTGLIBS    := $(shell root-config --glibs) -lEG -lFoam
ROOTCFLAGS   := $(shell root-config --cflags)
ROOTLDFLAGS  := $(shell root-config --ldflags)

vpath %.cxx $(S)
vpath %.h  $(I)
vpath %.o  $(O)

DEP_LIB      := libCore.so libHist.so libGraf.so libGpad.so libGui.so libPhysics.so libFoam.so libEG.so libTree.so

LIB_MLIB = $(L)/libMano.so
SOFLAGS = -shared
POST_LIB_BUILD =


# -------------------------------- Compile options --------------------------------

CCCOMP      = g++
CXXFLAGS    = -g -O3 -Wall -fPIC $(ROOTCFLAGS) $(AI) $(foreach dir,$(EXI),-I$(dir))
LDFLAGS     = -g -O3 $(ROOTLDFLAGS)

# ------------------------------------ Targets ------------------------------------

all:	begin $(OBJ) $(LIB_MLIB) $(L)/libMano.rootmap end

begin:
	@echo
	@echo "Building Mano"
	@echo

end:
	@echo "Done."

%.o: %.cxx
	@echo "Compiling $(patsubst %.cxx,%,$(notdir $<))"
	@mkdir -p $(O)
	@$(CCCOMP) $(CXXFLAGS) -o $(O)/$@ -c $< 

$(O)/G__MLIB.o: $(S)/G__MLIB.cxx
	@echo "Compiling Dict"
	@$(CCCOMP) $(CXXFLAGS) -o $(O)/G__MLIB.o -c $(S)/G__MLIB.cxx

$(S)/G__MLIB.cxx: $(INCD) $(I)/LinkDef.h 
	@echo "Creating MLIB dictionary"
	@rootcint -v -f $@ -c $(AI) -p $(INC) LinkDef.h

$(LIB_MLIB): $(OBJD)
	@echo
	@echo "Building libMano"
	@mkdir -p $(L)
	@rm -f $(L)/libMano.*
	@$(CCCOMP) $(LDFLAGS) $(SOFLAGS) $(OBJD) -o $(LIB_MLIB)
	@$(POST_LIB_BUILD)


$(L)/libMano.rootmap: $(LIB_MLIB)
	@echo "Creating ROOT map"
	@rlibmap -o $(L)/libMano.rootmap -l $(LIB_MLIB) -d $(DEP_LIB) -c $(I)/LinkDef.h

clean:
	@echo "Cleaning Mano distribution"
	@rm -r -f $(S)/G__*
	@rm -f -r $(O)
	@rm -r -f $(L)
	@echo "Done."
