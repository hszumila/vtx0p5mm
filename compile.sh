#!/bin/sh

# last argument is always sourcecode's filename to compile
# preceeding arguments will be passed to compiler
# output executable has same name as source w/o suffix

narg=${#@}
source=${@:$narg:$narg}
exe=${source%.*}
clopts=${@:1:$narg-1}

if [ ! -e "$source" ];
then
    echo "Usage:  compile.sh [compiler clopts] c-file"
    exit
fi

#OSX needs c++:
compiler=c++
#compiler=gcc
#compiler=g++
#compiler=gfortran

opts='-W -Wall '
#opts+='-Wshadow '
opts+='-Wstrict-aliasing '

opts+='-Wno-gnu-static-float-init '
opts+='-Wno-unused-private-field '

rootflags=`root-config --cflags --libs`

cmd="$compiler
     $opts
     $clopts
     $rootflags
     -I/Applications/root/build
     -o $exe $source"

echo $cmd
`$cmd`

