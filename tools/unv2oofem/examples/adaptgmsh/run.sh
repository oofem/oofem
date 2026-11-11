#!/bin/sh
# Adaptive analyses of a notched three-point bending beam with Gmsh remeshing.
# Usage: ./run.sh [2d|tet|hex|nl2d|nl3d|ls2d|lstet|lshex]   (default 2d)
#   2d, tet, hex       : notched beam, adaptive linear static analysis (adaptlinearstatic, ZZ error estimator)
#   nl2d, nl3d         : notched beam, adaptive nonlinear static analysis (adaptnlinearstatic, nonlocal damage)
#   ls2d, lstet, lshex : L-shaped body, adaptive linear static analysis (adaptlinearstatic, ZZ error estimator)
# Requires: gmsh executable and python3 on PATH, oofem executable (OOFEM env. variable or on PATH).
set -e
HERE=$(cd "$(dirname "$0")" && pwd)
GMSH2OOFEM="$HERE/../../gmsh2oofem.py"
OOFEM=${OOFEM:-oofem}

case "${1:-2d}" in
  2d)   GEO=beam.geo;   CTRL=beam.ctrl;        OPTS="--dim 2" ;;
  tet)  GEO=beam3d.geo; CTRL=beam3d_tet.ctrl;  OPTS="--dim 3" ;;
  hex)  GEO=beam3d.geo; CTRL=beam3d_hex.ctrl;  OPTS="--dim 3 --elemtype hex" ;;
  nl2d) GEO=beam.geo;   CTRL=beam_nl.ctrl;     OPTS="--dim 2" ;;
  nl3d) GEO=beam3d.geo; CTRL=beam3d_nl.ctrl;   OPTS="--dim 3" ;;
  ls2d)  GEO=lshape.geo;   CTRL=lshape.ctrl;       OPTS="--dim 2" ;;
  lstet) GEO=lshape3d.geo; CTRL=lshape3d_tet.ctrl; OPTS="--dim 3" ;;
  lshex) GEO=lshape3d.geo; CTRL=lshape3d_hex.ctrl; OPTS="--dim 3 --elemtype hex" ;;
  *)    echo "usage: $0 [2d|tet|hex|nl2d|nl3d|ls2d|lstet|lshex]"; exit 1 ;;
esac
IN="${CTRL%.ctrl}.in"

cd "$HERE"
# initial mesh generated from geometry
python3 "$GMSH2OOFEM" $OPTS "$GEO" "$CTRL" "$IN"
# adaptive analysis: analysis -> error estimation -> remeshing -> analysis ...
"$OOFEM" -f "$IN"
