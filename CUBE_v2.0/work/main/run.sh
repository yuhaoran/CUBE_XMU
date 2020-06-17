source ../utilities/module_load_brew.sh 

cd ../utilities/
make clean
make ic.x dsp.x cicpower.x fof.x lpt1.x cicrsd.x
./ic.x

cd ../main/
make clean
make
./main.x

cd ../utilities/
./dsp.x
./cicpower.x
#./fof.x
#./lpt1.x
#./cicrsd.x

cd ../main/
