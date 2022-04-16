rm -r build
mkdir build/
cd build
cmake ..
make

echo "DA-utils build done!\n"

cd ../../LETKF
rm -r build
mkdir build/
cd build
cmake ..
make

echo "Letkf build done!\n"
