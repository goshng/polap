# Variables
CC = gcc
CFLAGS = -Wall -g
SRC_DIR = .
BUILD_DIR = ../build/src
POLAP_DATA = polap polap-data-aflye polap-data-cflye polap-data-dflye polap-data-taxon bolap
SRCS = $(wildcard $(SRC_DIR)/*.c)
# OBJS = $(patsubst $(SRC_DIR)/%.c,$(BUILD_DIR)/%.o,$(SRCS),$(POLAP_DATA))
OBJS = $($(POLAP_DATA))

# Default target
# all: $(BUILD_DIR) $(OBJS)
all: $(BUILD_DIR) $(POLAP_DATA)
	chmod +x $(POLAP_DATA)

# Create the build directory if it doesn't exist
$(BUILD_DIR):
	mkdir -p $(BUILD_DIR)

# Rule to compile source files into object files in the build directory
$(BUILD_DIR)/%.o: $(SRC_DIR)/%.c
	$(CC) $(CFLAGS) -c $< -o $@

polap: polap.sh
	sed 's/_POLAP_RELEASE:=0/_POLAP_RELEASE:=1/' $< > $@

polap-data-aflye: polap-data-v1.sh polaplib/polap-data-v1.csv polaplib/polap-bash-figure-latex.sh
	sed 's/_POLAP_RELEASE:=0/_POLAP_RELEASE:=1/' $< > $@
	cp -p polaplib/polap-data-v1.csv polaplib/polap-data-aflye.csv

polap-data-cflye: polap-data-v2.sh polaplib/polap-data-v2.csv polaplib/polap-bash-figure-latex.sh
	sed 's/_POLAP_RELEASE:=0/_POLAP_RELEASE:=1/' $< > $@
	cp -p polaplib/polap-data-v2.csv polaplib/polap-data-cflye.csv

polap-data-taxon: polap-data-v3.sh polaplib/polap-data-v3.csv polaplib/polap-bash-figure-latex.sh
	sed 's/_POLAP_RELEASE:=0/_POLAP_RELEASE:=1/' $< > $@
	cp -p polaplib/polap-data-v3.csv polaplib/polap-data-taxon.csv

polap-data-dflye: polap-data-v4.sh polaplib/polap-data-v4.csv polaplib/polap-bash-figure-latex.sh
	sed 's/_POLAP_RELEASE:=0/_POLAP_RELEASE:=1/' $< > $@
	cp -p polaplib/polap-data-v4.csv polaplib/polap-data-dflye.csv

bolap: bolap.sh
	sed 's/_POLAP_RELEASE:=0/_POLAP_RELEASE:=1/' $< > $@

# Clean build files
clean:
	rm -rf $(BUILD_DIR) $(POLAP_DATA)
	rm -f polap*~
	@echo "Cleaned $(BUILD_DIR) in subdir1."
	@echo "Cleaned $(POLAP_DATA) in subdir1."

# Phony targets
.PHONY: all clean

