# Lxplus Batch Job Script

export WMASS_SO="XXX"
export PY_WRAPPER="YYY"
export EOS_DIR="ZZZ"
export BASEDIR="AAA"
export FILENAME="BBB"
export OUTDIR="CCC"

echo "BATCH: i am here"
$PWD
echo "BATCH: this is in the directory"
$LS
echo "BATCH: getting the shared object..."
cp $WMASS_SO .
echo "ensuring the output directory exists"
mkdir -p $BASEDIR/$OUTDIR
echo "BATCH: copying the executable here"
cp $BASEDIR/$PY_WRAPPER .
echo "BATCH: starting the job"
python $PY_WRAPPER $EOS_DIR  -o $BASEDIR/$OUTDIR/ -f $FILENAME
