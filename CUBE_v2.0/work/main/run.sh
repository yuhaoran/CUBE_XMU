source ../utilities/module_load_brew.sh 

cd ../utilities/
make clean
make ic.x cicpower.x
./ic.x

cd ../main/
make clean
make
./main.x

cd ../utilities/
./cicpower.x
source run_fof.sh

cd ../main/
