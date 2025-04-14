#
# justfile for various local build/test tasks; works on Debian stanble
#
wheel:
	python -m build
	pip3 install --force-reinstall dist/oofem-2.5.0.dev1-cp312-cp312-linux_x86_64.whl
	cd bindings/python/tests; python -m pytest
pyodide:
	#!/bin/bash
	pip3 install pyodide-build
	[ -d build/emsdk ] || git clone https://github.com/emscripten-core/emsdk.git build/emsdk
	cd build/emsdk/
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
	cmake -Bbuild -H. -GNinja -DCMAKE_BUILD_TYPE=Release -DUSE_PYBIND_BINDINGS=1 -DUSE_SHARED_LIB=1 -DUSE_OOFEM_EXE=1 -DUSE_SM=1 -DUSE_TM=1 -DUSE_MPM=1
	ninja -C build/
	ctest --test-dir build/ --parallel=16
notshared:
	mkdir -p build
	cmake -Bbuild -H. -GNinja -DCMAKE_BUILD_TYPE=Release -DUSE_PYBIND_BINDINGS=1 -DUSE_SHARED_LIB=0 -DUSE_OOFEM_EXE=1 -DUSE_SM=1 -DUSE_TM=1 -DUSE_MPM=1
	ninja -C build/
	ctest --test-dir build/ --parallel=16
act-install:
	#!/bin/bash
	[ -f nektos-act ] || wget https://github.com/nektos/act/releases/download/v0.2.76/act_Linux_x86_64.tar.gz && tar xvfz act_Linux_x86_64.tar.gz act && mv ./act nektos-act && chmod a+x nektos-act
act-linux:
	./nektos-act -P ubuntu-latest=ghcr.io/catthehacker/ubuntu:act-24.04 --workflows ./.github/workflows/linux.yml
act-python:
	./nektos-act --reuse -P ubuntu-latest=ghcr.io/catthehacker/ubuntu:act-24.04 --workflows ./.github/workflows/python.yml
act-cibuildwheel:
	./nektos-act --reuse -P ubuntu-latest=ghcr.io/catthehacker/ubuntu:act-24.04 --workflows ./.github/workflows/cibuildhweel.yml
msvc:
	#!/bin/bash
	[ -d /opt/msvc ] || ( git clone https://github.com/mstorsjo/msvc-wine.git && cd msvc-wine && ./vsdownload.py --accept-license --dest /opt/msvc && ./install.sh )
	export PATH=/opt/msvc/bin/x64:$PATH
	rm -rf build-msvc
	cmake -Bbuild-msvc -GNinja -DCMAKE_BUILD_TYPE=Release -DCMAKE_SYSTEM_NAME=Windows -DCMAKE_C_COMPILER=cl -DCMAKE_CXX_COMPILER=cl -DCMAKE_CXX_FLAGS="/wd4275 /wd4267 /wd4458 /wd4456 /wd5205 /wd4244 /wd4101 /EHsc"
	cmake --build ./build-msvc --parallel
