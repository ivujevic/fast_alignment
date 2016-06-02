.PHONY: \
    all \
    build \
    build-coverage \
    clean \
    compile \
    rebuild \
    test \
    update \
    packages

all: compile
	@:

clean:
	rm -rf build build_debug build_cov

build_cov/Makefile:
	mkdir -p build_cov && cd build_cov && cmake -DCODE_COVERAGE=ON ..

build_cov: build_cov/Makefile
	@:

build/Makefile:
	mkdir -p build && cd build && cmake .. -DCMAKE_CXX_COMPILER=g++-4.9 -DCMAKE_C_COMPILER=gcc-4.9

compile: update build/Makefile
	@make -C build -j$$NBUILDJOBS

build_debug/Makefile:
	mkdir -p build_debug && cd build_debug && cmake -DCMAKE_BUILD_TYPE=Debug ..

debug: update build_debug/Makefile
	@make -C build_debug -j$$NBUILDJOBS

rebuild: clean compile
	@:

test: build
	make -C build test

packages:
	bash ./deps/install_deps.sh
