echo "Make sure cmake is installed and is available in PATH"
wget http://lemon.cs.elte.hu/pub/sources/lemon-1.3.1.tar.gz
LEMON_SRC="lemon-1.3.1"
LEMON_DIR="lemon"
WORK_DIR=`pwd`
mkdir -p ${LEMON_DIR}
tar -zvxf ${LEMON_SRC}.tar.gz

cd ${LEMON_SRC}
cmake . -DCMAKE_INSTALL_PREFIX="../${LEMON_DIR}"
make 
make install
cd ${WORK_DIR}

echo "Compiling now"
g++ solve_mcf.cc -I ${LEMON_DIR}/include/ -L ${LEMON_DIR}/lib/ -lemon -o mcf
