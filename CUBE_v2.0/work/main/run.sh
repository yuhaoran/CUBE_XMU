source ../utilities/module_load_brew.sh 

cd ../utilities/
make clean
make ic.x cicpower.x fof.x lpt1.x
./ic.x

cd ../main/
make clean
make
./main.x

cd ../utilities/
./cicpower.x
./fof.x
./lpt1.x

cd ../main/
