# Install GDL (in gdl-1.2/GDL) 
cd ext
tar -zxf gdl-1.1.tar.gz
cd gdl-1.1
sh autogen.sh
./configure --prefix=$PWD/GDL/
make
make install
cd ../..

# Install tools 
make
