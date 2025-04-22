#
# justfile for various local build/test tasks; works on Debian stable
#
wheel:
	python -m build
wheel-test:
	uv --version || $( curl -LsSf https://astral.sh/uv/install.sh | sh )
	rm -rf .cache/venv
	uv venv --python 3.12 .cache/venv
	cd .cache/venv; uv pip install `ls ../../dist/oofem-*.whl |tail -n1`; . bin/activate; cd ../../bindings/python/tests; python -m pytest
pyodide:
	#!/bin/bash
	pip3 install pyodide-build
	[ -d .cache/emsdk ] || git clone https://github.com/emscripten-core/emsdk.git build/emsdk
	cd .cache/emsdk/
	export PYODIDE_EMSCRIPTEN_VERSION=$(pyodide config get emscripten_version)
	./emsdk install ${PYODIDE_EMSCRIPTEN_VERSION}
	./emsdk activate ${PYODIDE_EMSCRIPTEN_VERSION}
	source emsdk_env.sh
	cd ../..
	pyodide build
	[ -d build/venv-pyodide ] || pyodide venv build/venv-pyodide
	source build/venv-pyodide/bin/activate
	pip install --force-reinstall dist/oofem-*-pyodide_*_wasm32.whl
	pip install pytest numpy
	cd bindings/python/tests; python -m pytest
shared:
	mkdir -p build
	cmake -Bbuild -H. -GNinja -DCMAKE_BUILD_TYPE=Release -DUSE_PYBIND_BINDINGS=1 -DUSE_PYTHON_EXTENSION=1 -DUSE_SHARED_LIB=1 -DUSE_OOFEM_EXE=1 -DUSE_SM=1 -DUSE_TM=1 -DUSE_MPM=1
	ninja -C build/
	ctest --test-dir build/ --parallel=16 --output-on-failure
notshared:
	mkdir -p build
	cmake -Bbuild -H. -GNinja -DCMAKE_BUILD_TYPE=Release -DUSE_PYBIND_BINDINGS=1 -DUSE_SHARED_LIB=0 -DUSE_OOFEM_EXE=1 -DUSE_SM=1 -DUSE_TM=1 -DUSE_MPM=1
	ninja -C build/
	ctest --test-dir build/ --parallel=16
nanobind:
	cmake -Bbuild-nanobind -H. -GNinja -DCMAKE_BUILD_TYPE=Release -DUSE_PYBIND_BINDINGS=1 -DUSE_NANOBIND=1 -DUSE_SHARED_LIB=0 -DUSE_OOFEM_EXE=1 -DUSE_SM=1 -DUSE_TM=1 -DUSE_MPM=1
	ninja -C build-nanobind
	ctest --test-dir build/ --output-on-failure --parallel=16
act-install:
	#!/bin/bash
	[ -f .cache/nektos-act ] || mkdir -p .cache && wget https://github.com/nektos/act/releases/download/v0.2.76/act_Linux_x86_64.tar.gz && tar xvfz act_Linux_x86_64.tar.gz act && rm act_Linux_x86_64.tar.gz && mv ./act .cache/nektos-act && chmod a+x .cache/nektos-act
act-linux:
	.cache/nektos-act -P ubuntu-latest=ghcr.io/catthehacker/ubuntu:act-24.04 --workflows ./.github/workflows/linux.yml
act-python:
	.cache/nektos-act --reuse -P ubuntu-latest=ghcr.io/catthehacker/ubuntu:act-24.04 --workflows ./.github/workflows/python.yml
act-cibuildwheel:
	.cache/nektos-act --reuse -P ubuntu-latest=ghcr.io/catthehacker/ubuntu:act-24.04 --workflows ./.github/workflows/cibuildhweel.yml
msvc:
	#!/bin/bash
	[ -d /opt/msvc ] || ( git clone https://github.com/mstorsjo/msvc-wine.git && cd msvc-wine && ./vsdownload.py --accept-license --dest /opt/msvc && ./install.sh )
	export PATH=/opt/msvc/bin/x64:$PATH
	rm -rf build-msvc
	cmake -Bbuild-msvc -GNinja -DCMAKE_BUILD_TYPE=Release -DCMAKE_SYSTEM_NAME=Windows -DCMAKE_C_COMPILER=cl -DCMAKE_CXX_COMPILER=cl
	cmake --build ./build-msvc --parallel
pytest:
	#!/bin/bash
	PYTHONPATH=build:bindings/python python -m pytest bindings/python/tests
test-ext:
	# ctest --test-dir build/ -VV --output-on-failure -R test_sm_python_usrdefboundaryload01.in
	gdb -ex=run -args build/oofem -f tests/sm/python/usrdefboundaryload01.in
