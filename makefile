#Compiler and compiler options
FC = gfortran
FCFLAGS = -fbounds-check -O0 -I $(MOD_DIR) 
LDFLAGS = -llapack 

#Definitions of directorys
SRC_DIR = src/
BLD_DIR = build/
MOD_DIR = mod/

#Source file decrlarations
SRC_FILES = atom_mod.f03 utils.f03 num_alg.f03 main.f03 
SRC = $(addprefix $(SRC_DIR),$(SRC_FILES))

#Object file declarations
O_FILES = $(patsubst $(SRC_DIR)%.f03, $(BLD_DIR)%.o, $(SRC))

#Executable name
EXECUTABLE = guass

$(BLD_DIR)$(EXECUTABLE) : $(O_FILES)
	$(FC) $^ $(LDFLAGS) $(FCFLAGS) -o $@ 
	mv $@ .

$(O_FILES) : $(BLD_DIR)%.o : $(SRC_DIR)%.f03
	$(FC)  -c $< -o $@ $(FCFLAGS) $(LDFLAGS) -J $(MOD_DIR) 

clean:
	rm -f $(BLD_DIR)*.o
	rm -f $(MOD_DIR)*.mod
	rm -f $(EXECUTABLE)
