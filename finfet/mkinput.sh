#/usr/bin/bash -w

#make link file from the input subdirectory.
#input parameter as the pathname .

grid_target=$1"/grid.txt"
config_target=$1"/dg.txt"
plot_target=$1"/plot3d"

if [ -f ${grid_target} ] && [ -f ${config_target} ]
then
       `ln -sf ${grid_target} lgrid.txt`
       `ln -sf ${config_target} ldg.txt`
       `cp -r ${plot_target} plot3d`
fi
