FC = gfortran
FFLAGS = -J$(MOD_DIR) -Wall -O2
LDFLAGS = -lfftw3 -lfftw3_threads

# Directories
SRC_DIR = src
TEST_DIR = test
BUILD_DIR = build
OBJECTS_DIR = $(BUILD_DIR)/obj
OBJECTS_TEST_DIR = $(BUILD_DIR)/obj_test
MOD_DIR = $(BUILD_DIR)/mod
TEST_BINS = $(BUILD_DIR)/test
DIRS = $(BUILD_DIR) $(OBJECTS_DIR) $(MOD_DIR) $(OBJECTS_TEST_DIR) $(TEST_BINS)


# List of all source files
SRCS := $(SRC_DIR)/NAFPack_constant.f90 \
		$(SRC_DIR)/Array_creation/NAFPack_meshgrid.f90 \
		$(SRC_DIR)/algebra/NAFPack_matricielle.f90 \
		$(SRC_DIR)/algebra/NAFPack_matrix_decomposition.f90 \
		$(SRC_DIR)/algebra/NAFPack_linear_algebra.f90 \
		$(SRC_DIR)/FFT/FFTW3.f90 \
		$(SRC_DIR)/FFT/NAFPack_fft.f90

# List of all object files
OBJS := $(OBJECTS_DIR)/NAFPack_constant.o \
		$(OBJECTS_DIR)/NAFPack_meshgrid.o \
		$(OBJECTS_DIR)/NAFPack_matricielle.o \
		$(OBJECTS_DIR)/NAFPack_matrix_decomposition.o \
		$(OBJECTS_DIR)/NAFPack_linear_algebra.o \
		$(OBJECTS_DIR)/FFTW3.o \
		$(OBJECTS_DIR)/NAFPack_fft.o


# List of all test files
TEST_SRCS = $(TEST_DIR)/test_NAFPack_linear_algebra.f90 \
			$(TEST_DIR)/test_NAFPack_fft.f90 \
			$(TEST_DIR)/main_test.f90 \

# List of all test object files
TEST_OBJS = $(OBJECTS_TEST_DIR)/test_NAFPack_linear_algebra.o \
			$(OBJECTS_TEST_DIR)/test_NAFPack_fft.o \
			$(OBJECTS_TEST_DIR)/main_test.o \

# Build the main executable
modules: check_fftw $(OBJS)

# Create the directories if they do not exist
$(DIRS):
	mkdir -p $@

# Check if FFTW library is installed
check_fftw:
	@echo "Checking for FFTW library..."
ifeq ($(OS),Windows_NT)
	@./check_fftw/check_fftw.bat
else
	@sh ./check_fftw/check_fftw.sh
endif

$(OBJECTS_DIR)/%.o: $(SRC_DIR)/%.f90 | $(OBJECTS_DIR) $(MOD_DIR)
	$(FC) $(FFLAGS) -c $< -o $@

$(OBJECTS_DIR)/%.o: $(SRC_DIR)/*/%.f90 | $(OBJECTS_DIR) $(MOD_DIR)
	$(FC) $(FFLAGS) -c $< -o $@


# Build the test object files
test: modules $(TEST_OBJS) $(TEST_BINS)/main_test.exe

$(OBJECTS_TEST_DIR)/%.o: $(TEST_DIR)/%.f90 | $(OBJECTS_TEST_DIR)
	$(FC) $(FFLAGS) -c $< -o $@

$(TEST_BINS)/main_test.exe: $(TEST_OBJS) $(OBJS) | $(TEST_BINS)
	$(FC) $^ -o $@ $(LDFLAGS)

# Run the tests
run_test: test
	$(TEST_BINS)/main_test.exe

# Clean up the build directories
clean:
	rm -rf $(BUILD_DIR)

.PHONY: modules test run_test clean check_fftw