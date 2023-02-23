#======================================================================
# Global Makefile 
#======================================================================

# directories
MAIN_DIR  = $(PWD)
BUILD_DIR = $(MAIN_DIR)/build
SRC_DIR   = $(MAIN_DIR)
OBJ_DIR   = $(BUILD_DIR)/obj
DEP_DIR   = $(BUILD_DIR)/dep
INC_DIR   = $(MCM_ROOT_DIR) $(EIGEN_ROOT_DIR) $(MKLROOT)/include
LIB_DIR   = $(MKLROOT)/lib/intel64

# extension
SRC_EXT = .cc
OBJ_EXT = .o
DEP_EXT = .d

# source, objective, dependency and executable files for unit tests
SRC_FILE  = $(wildcard $(SRC_DIR)/*$(SRC_EXT))
OBJ_FILE  = $(SRC_FILE:$(SRC_DIR)/%$(SRC_EXT)=$(OBJ_DIR)/%$(OBJ_EXT))
DEP_FILE  = $(SRC_FILE:$(SRC_DIR)/%$(SRC_EXT)=$(DEP_DIR)/%$(DEP_EXT))
EXEC_FILE = $(MAIN_DIR)/main

# compile and link settings
CXX      := g++
CXXFLAGS  = -std=c++17 -Wall -O3 -g
CXXFLAGS += -MT $@ -MMD -MP -MF $(DEP_DIR)/$*$(DEP_EXT)
CXXFLAGS += $(INC_DIR:%=-I %)
CXXFLAGS += -fopenmp
CXXFLAGS += -DWITHOUT_NUMPY -DMKL_ILP64
# CXXFLAGS += -DMCM_DEBUG
LDFLAGS   = $(LIB_DIR:%=-L %)
LDFLAGS  += -lmcm
LDFLAGS  += -m64 -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -liomp5 -lpthread -lm -ldl

all: $(EXEC_FILE) | build

.PHONY: build
build:
	@mkdir -p $(BUILD_DIR) $(OBJ_DIR) $(DEP_DIR)

-include $(DEP_FILE)

.PHONY: clean
clean:
	@rm -f $(EXEC_FILE)
	@rm -rf $(BUILD_DIR)

.PHONY: info
info:
	@echo "Source directory for test:\n$(DIR)\n"
	@echo "Build directory for test:\n$(BUILD_DIR)\n"
	@echo "Object directory for test:\n$(OBJ_DIR)\n"
	@echo "Dependency directory for test:\n$(DEP_DIR)\n"
	@echo "Source files for test:\n $(SRC_FILE:%=%\n)"
	@echo "Object files for test:\n $(OBJ_FILE:%=%\n)"
	@echo "Dependency files for test:\n $(DEP_FILE:%=%\n)"
	@echo "Exectuable file for test:\n $(EXEC_FILE:%=%\n)"

# linked object files to get an executable file
$(EXEC_FILE): $(OBJ_FILE) 
	@echo "[Link all objective files to executable file $@]"
	@$(CXX) -o $@ $^ $(LDFLAGS)

# compile all source files in test directory to get object and dependency files
$(OBJ_DIR)/%$(OBJ_EXT): $(SRC_DIR)/%$(SRC_EXT) | build
	@echo "[Compile $<]"
	@$(CXX) $(CXXFLAGS) -c $< -o $@ 

#======================================================================
# Fast numerical tests
#======================================================================
DATA_DIR = $(MAIN_DIR)/data

.PHONY: run
run:
	@make && $(MAIN_DIR)/main -Nx 4 -Ny 4 -T 1.00 -CFL 0.66 -O 1 
	# @gnuplot -p -c $(DATA_DIR)/plot_u_2D.gp

#======================================================================
# Git macros
#======================================================================
REPOS_NAME = linear_advection_2d_sldg
REPOS_URL = git@github.com:escapetiger/$(REPOS_NAME).git
REPOS_BRANCH = main
REPOS_REMOTE = origin

COMMIT_MSG = Updated at $(shell date +"%Y.%m.%d")

.PHONY: git_init
git_init:
	@echo "# $(REPOS_NAME)" >> README.md
	@git init
	@git add README.md
	@git commit -m "first commit"
	@git branch -M $(REPOS_BRANCH)
	@git remote add $(REPOS_REMOTE) $(REPOS_URL)
	@git push -u $(REPOS_REMOTE) $(REPOS_BRANCH)

.PHONY: git_reset
git_reset:
	@git reset HEAD~ 

.PHONY: git_add_commit
git_add_commit:
	@git add .
	@git commit -m "$(COMMIT_MSG)"

.PHONY: git_push
git_push:
	@git push -u $(REPOS_REMOTE) $(REPOS_BRANCH)

.PHONY: git_l2w
git_l2w:
	@git add .
	@git commit -m "$(COMMIT_MSG)"
	@git push $(REPOS_REMOTE)

.PHONY: git_w2l
git_w2l:
	@git fetch $(REPOS_REMOTE) $(REPOS_BRANCH)
	@git log $(REPOS_BRANCH).. $(REPOS_REMOTE)/$(REPOS_BRANCH)
	@git merge $(REPOS_REMOTE)/$(REPOS_BRANCH)