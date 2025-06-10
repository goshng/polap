CC = gcc
CFLAGS = -Wall -g
SRC_DIR = src
BUILD_DIR = build
SUBDIRS = src
# SRCS = $(wildcard $(SRC_DIR)/*.c)
# OBJS = $(patsubst $(SRC_DIR)/%.c,$(BUILD_DIR)/%.o,$(SRCS))
TARGET = polap

# Default target
# all: $(BUILD_DIR) $(SUBDIRS) $(TARGET)
all: $(BUILD_DIR) $(SUBDIRS)

# Create the build directory if it doesn't exist
$(BUILD_DIR):
	mkdir -p $(BUILD_DIR)

# Build the target program
$(TARGET): $(OBJS)
	$(CC) $(CFLAGS) -o $@ $^

# Rule to compile source files into object files in the build directory
$(BUILD_DIR)/%.o: $(SRC_DIR)/%.c
	$(CC) $(CFLAGS) -c $< -o $@

# Recursively call Makefiles in subdirectories
$(SUBDIRS):
	@for dir in $(SUBDIRS); do \
		$(MAKE) -C $$dir; \
	done

# Clean all build files
clean:
	rm -rf $(BUILD_DIR) $(TARGET)
	@for dir in $(SUBDIRS); do \
		$(MAKE) -C $$dir clean; \
	done
	@echo "Cleaned all build artifacts."

# Phony targets
.PHONY: all clean $(SUBDIRS)

