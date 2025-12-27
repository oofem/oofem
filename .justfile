#
# justfile for various local build/test tasks; works on Debian stable
#
FAIRLY_COMPLETE_FLAGS := "-DUSE_XML=1 -DUSE_SHARED_LIB=0 -DUSE_OOFEM_EXE=1 -DUSE_SM=1 -DUSE_TM=1 -DUSE_MPM=1 -DUSE_TRACE_FIELDS=1 -DCMAKE_BUILD_TYPE=RelWithDebInfo -GNinja"

wheel:
	python -m build
wheel-test:
	uv --version || $( curl -LsSf https://astral.sh/uv/install.sh | sh )
	rm -rf .cache/venv
	uv venv --python 3.12 .cache/venv
	cd .cache/venv; uv pip install `ls ../../dist/oofem-*.whl |tail -n1`; . bin/activate; cd ../../bindings/python/tests; python -m pytest
pyodide:
	#!/bin/bash
	set -e -x
	# remove old wheels
	rm -f dist/oofem-*-pyodide_*_wasm32.whl 
	uv --version || $( curl -LsSf https://astral.sh/uv/install.sh | sh )
	# create venv (use uv since we can easily specify python version which is not installed locally)
	UV_VENV_CLEAR=1 uv venv --python 3.13 .cache/pyodide-venv
	source .cache/pyodide-venv/bin/activate
	# install pyodide-build and pip
	uv pip install pyodide-build pip
	# install emscripten SDK
	[ -d .cache/emsdk ] || git clone https://github.com/emscripten-core/emsdk.git .cache/emsdk
	cd .cache/emsdk/
	export PYODIDE_EMSCRIPTEN_VERSION=$(pyodide config get emscripten_version)
	./emsdk install ${PYODIDE_EMSCRIPTEN_VERSION}
	# activate emscripten SDK
	./emsdk activate ${PYODIDE_EMSCRIPTEN_VERSION}
	source emsdk_env.sh
	cd ../..
	pyodide build
	# create & activate runtime venv (≠ pyodide-build venv)
	# this must be done inside the .cache/pyodide-venv venv, since that's where pyodide command lives
	[ -d build/venv-pyodide ] || pyodide venv build/venv-pyodide
	# exit .cache/pyodide-venv before proceeding
	deactivate
	# enter the wasm32 venv
	source build/venv-pyodide/bin/activate
	python -m pip install --force-reinstall dist/oofem-*-pyodide_*_wasm32.whl
	cd bindings/python/tests; python -m pytest
	deactivate # exit the wasm32 venv
shared:
	mkdir -p build-shared
	cmake -Bbuild-shared -H. -GNinja  -DCMAKE_EXPORT_COMPILE_COMMANDS=1 -DCMAKE_BUILD_TYPE=Release -DUSE_PYBIND_BINDINGS=1 -DUSE_PYTHON_EXTENSION=1 -DUSE_SHARED_LIB=1 -DUSE_OOFEM_EXE=1 -DUSE_SM=1 -DUSE_TM=1 -DUSE_MPM=1 -DUSE_XML=1
	ninja -C build-shared/
	ctest --test-dir build-shared/ --parallel 16 --output-on-failure
eigen:
	mkdir -p build-eigen
	cmake -Bbuild-eigen -H. -DUSE_LAPACK=0 -DUSE_EIGEN=1 -DUSE_OPENMP_PARALLEL=0 {{FAIRLY_COMPLETE_FLAGS}}
	ninja -C build-eigen/
	ctest --test-dir build-eigen/ --parallel 16 --output-on-failure
openmp:
	mkdir -p build-openmp
	cmake -Bbuild-openmp -H. -DUSE_LAPACK=0 -DUSE_EIGEN=1 -DUSE_OPENMP_PARALLEL=1 {{FAIRLY_COMPLETE_FLAGS}}
	ninja -C build-openmp/
	ctest --test-dir build-openmp/ --parallel 16 --output-on-failure
lapack:
	# those don't work with BLAS (yet): -DUSE_PYBIND_BINDINGS=0 -DUSE_PYTHON_EXTENSION=0 -DUSE_LAPACK=1
	mkdir -p build-lapack
	cmake -Bbuild-lapack -H. -DUSE_LAPACK=1 -DUSE_EIGEN=0 -DUSE_OPENMP_PARALLEL=0 {{FAIRLY_COMPLETE_FLAGS}}
	ninja -C build-lapack/
	ctest --test-dir build-lapack/ --parallel 16 # --verbose --output-on-failure
callgrind:
	valgrind --tool=callgrind --callgrind-out-file=callgrind.lapack build-lapack/oofem -f tests/sm/spring01.in
	valgrind --tool=callgrind --callgrind-out-file=callgrind.eigen build-eigen/oofem -f tests/sm/spring01.in
	valgrind --tool=callgrind --callgrind-out-file=callgrind.openmp build-openmp/oofem -t1 -f tests/sm/spring01.in
	kcachegrind callgrind.lapack &
	kcachegrind callgrind.eigen &
	kcachegrind callgrind.openmp &
static:
	# those don't work with BLAS (yet): -DUSE_PYBIND_BINDINGS=0 -DUSE_PYTHON_EXTENSION=0 -DUSE_LAPACK=1
	mkdir -p build-static
	cmake -Bbuild-static -H. -DUSE_LAPACK=1 -DUSE_EIGEN=0 -DUSE_OPENMP_PARALLEL=0 {{FAIRLY_COMPLETE_FLAGS}}
	ninja -C build-static/
	ctest --test-dir build-static/ --parallel 16 # --verbose --output-on-failure
nanobind:
	cmake -Bbuild-nanobind -H. -GNinja -DCMAKE_BUILD_TYPE=Release -DUSE_PYBIND_BINDINGS=1 -DUSE_NANOBIND=1 -DUSE_SHARED_LIB=0 -DUSE_OOFEM_EXE=1 -DUSE_SM=1 -DUSE_TM=1 -DUSE_MPM=1
	ninja -C build-nanobind
	ctest --test-dir build-nanobind/ --output-on-failure --parallel=16
act-install:
	#!/bin/bash
	[ -f .cache/nektos-act ] || mkdir -p .cache && wget https://github.com/nektos/act/releases/download/v0.2.76/act_Linux_x86_64.tar.gz && tar xvfz act_Linux_x86_64.tar.gz act && rm act_Linux_x86_64.tar.gz && mv ./act .cache/nektos-act && chmod a+x .cache/nektos-act
act-linux:
	.cache/nektos-act -P ubuntu-latest=ghcr.io/catthehacker/ubuntu:act-24.04 --workflows ./.github/workflows/linux.yml
act-full:
	.cache/nektos-act --reuse -P ubuntu-latest=ghcr.io/catthehacker/ubuntu:act-24.04 --workflows ./.github/workflows/full.yml
act-python:
	.cache/nektos-act --reuse -P ubuntu-latest=ghcr.io/catthehacker/ubuntu:act-24.04 --workflows ./.github/workflows/python.yml
act-cibuildwheel:
	.cache/nektos-act --reuse -P ubuntu-latest=ghcr.io/catthehacker/ubuntu:act-24.04 --workflows ./.github/workflows/cibuildhweel.yml
msvc:
	#!/bin/bash
	[ -d /opt/msvc ] || ( git clone https://github.com/mstorsjo/msvc-wine.git && cd msvc-wine && ./vsdownload.py --accept-license --dest /opt/msvc && ./install.sh /opt/msvc )
	export PATH=/opt/msvc/bin/x64:$PATH
	rm -rf build-msvc
	cmake -Bbuild-msvc -GNinja -DCMAKE_BUILD_TYPE=Release -DCMAKE_SYSTEM_NAME=Windows -DCMAKE_C_COMPILER=cl -DCMAKE_CXX_COMPILER=cl
	cmake --build ./build-msvc --parallel
pytest:
	#!/bin/bash
	PYTHONPATH=build:bindings/python python -m pytest bindings/python/tests
test-ext:
	# ctest --test-dir build/ -VV --output-on-failure -R test_sm_python_usrdefboundaryload01.in
	gdb -ex=run -args build-shared/oofem -f tests/sm/python/usrdefboundaryload01.in
xml:
	ninja -C build-static
	# gdb -ex=run -args
	# gdb -ex=run -args
	build-static/oofem -f tests/sm/spring01.in
	build-static/oofem -f tests/sm/spring01.xml
mpm:
	ninja -C build-eigen
	build-eigen/oofem -f tests/mpm/cook2_u1p0_2.in
	echo '------------------------------------- XML -------------------------------'
	build-eigen/oofem -f tests/mpm/cook2_u1p0_2.xml
	diff -I '^User time consumed by solution step .*' ./cook2_u1p0_2.out ./cook2_u1p0_2.xml.out && echo "NO DIFFERENCE :)"
mpm-gdb:
	ninja -C build-eigen
	DEBUGINFOD_URLS= gdb -ex=run -args build-eigen/oofem -f tests/mpm/cook2_u1p0_2.xml
