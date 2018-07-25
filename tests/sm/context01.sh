#
# this test check basic save/restore context functionality
#
OOFEM=$1
echo "target executable: $OOFEM"
pwd

echo "Command: $OOFEM -f context01.in.0 -c"
# run target on input and store context file
$OOFEM -f context01.in.0 -c
echo "Command: $OOFEM -f context01.in.0 -r 1"
# run target on the same file, but restarting from step 2
$OOFEM -f context01.in.0 -r 1

