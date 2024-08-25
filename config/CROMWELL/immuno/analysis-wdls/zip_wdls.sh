SRC_DIR="$(dirname $0)"
OLD_DIR=$PWD

cd $SRC_DIR/definitions/
zip -r ../workflows.zip . 1> /dev/null
cd $OLD_DIR
