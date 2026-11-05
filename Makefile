CONFIG ?= Debug

# Forward dependency overrides only when given on the make command line;
# an ambient shell value is ignored.
ifeq ($(origin FORTIO_REF),command line)
  _FORTIO_REF_ARG := -DFORTIO_REF=$(FORTIO_REF)
endif
ifeq ($(origin LIBNEO_REF),command line)
  _LIBNEO_REF_ARG := -DLIBNEO_REF=$(LIBNEO_REF)
endif
ifeq ($(origin LIBNEO_PATH),command line)
  _LIBNEO_PATH_ARG := -DLIBNEO_PATH=$(LIBNEO_PATH)
endif

.PHONY: all build test test_failed install clean plot golden golden_run dist internal python
all: build

build/CMakeCache.txt:
	cmake -S . -B build -DCMAKE_BUILD_TYPE=$(CONFIG) $(_FORTIO_REF_ARG) $(_LIBNEO_REF_ARG) $(_LIBNEO_PATH_ARG)

build: build/CMakeCache.txt
	cmake --build build

install: build
	cmake --install build

test: build
	ctest --test-dir build/test --output-on-failure -L quick

test_failed: build
	ctest --test-dir build/test --output-on-failure -V --rerun-failed

test_slow: build
	ctest --test-dir build/test --output-on-failure -L slow

test_all: build
	ctest --test-dir build/test --output-on-failure -L "quick|slow"

plot: build
	ctest --test-dir build/test --output-on-failure -L plot

current: build
	ctest --test-dir build/test --output-on-failure -V -L current

golden: build
	ctest --test-dir build/test --output-on-failure -V -L golden

golden_run: build
	ctest --test-dir build/test --output-on-failure -V -R GoldenRecordRun

golden_update: build/test/golden/rabe.nc build/test/golden/electric/electric_rabe.nc
	cp $^ test/golden/expected/.

clean:
	rm -rf build

dist:
	git archive --format=tar.gz --prefix=rabe/ HEAD -o rabe.tar.gz

python:
	cmake -S . -B build -DCMAKE_BUILD_TYPE=Release -DBUILD_PYTHON_BINDINGS=ON
	cmake --build build

# Private test data and external comparisons (not available in public repo)
internal: build
	ctest --test-dir build/test -V -L internal
