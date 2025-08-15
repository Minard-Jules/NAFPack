.PHONY: all win linux clean

ifeq ($(OS),Windows_NT)
CP = cp
RM = rm --r --fo
CONFIG_FILE = fpm_configs\fpm_windows.toml
C_FLAGS = 
F_FLAGS = 
else
CP = cp
RM = rm -rf
CONFIG_FILE = fpm_configs/fpm_linux.toml
C_FLAGS = -fPIE
F_FLAGS = 
endif

all: build

build:
	$(CP) $(CONFIG_FILE) fpm.toml
	fpm build --flag "$(F_FLAGS)" --c-flag "$(C_FLAGS)"

test: build
	fpm test --flag "$(F_FLAGS)" --c-flag "$(C_FLAGS)"


clean: clean_build clean_fpm

clean_build:
	$(RM) build

clean_fpm:
	$(RM) fpm.toml