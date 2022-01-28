#!/bin/bash

#
#
#  This file is part of the Sky Background software package.
#  Copyright (C) 2009-2018 European Southern Observatory
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA
#
#   Author:      Wolfgang Kausch & ESO in-kind team Innsbruck
#   Since:       01 Jun 2010
#   Last update: 06 Oct 2015
#

# defining console output text properties
txtbld=$(tput bold)       # Bold
txtrst=$(tput sgr0)       # Text reset
txtred=$(tput setaf 1)    # Red
txtgrn=$(tput setaf 2)    # Green
txtylw=$(tput setaf 3)    # Yellow
sourcedir=`pwd`
# preparation
echo "${txtbld} "
echo "*************************************************************************"
echo "*  Welcome to the SM-01 'Advanced Sky Background Model' setup routine   *"
echo "*                                                                       *"
echo "*                                                                       *"
echo "*                                                                       *"
echo "*  Austrian ESO in-kind Team Innsbruck                                  *"
echo "*  M. Barden, W. Kausch, S.Kimeswenger (head), S. Noll                  *"
echo "*  Institute of Astro- and Particlephysics                              *"
echo "*  University of Innsbruck / Austria                                    *"
echo "*                                                                       *"
echo "*************************************************************************"
echo "${txtrst} "

while [ ! $ANSWER != "n" ]
do
    echo -n "Choose installation directory of SM-01 software: "
    read -e INST_DIR
    LEN=$(echo ${#INST_DIR})
    echo -n "Is this OK [Y/n]? "
    read -e ANSWER
    if [ "$ANSWER" == "" -o "$ANSWER" == "y" -o "$ANSWER" == "Y" ]; then
         ANSWER="y"
    else
        ANSWER="n"
    fi
done

# checking write permissions
BASEDIR=`dirname $INST_DIR`

if  [ -w $BASEDIR ]; then
    echo " " >/dev/null
else
    echo "${txtbld}${txtred}ERROR: Cannot access '$INST_DIR' (no permissions) - Stop.${txtrst}"
    exit 1
fi

# creating directories
if [ -d $INST_DIR ]; then
    if [ -d $INST_DIR/sm-01_mod1 -o -d $INST_DIR/sm-01_mod2 ]; then
        echo "${txtbld}${txtred}Oops, seems that the installation already exists...."
        echo -n "Proceed (will overwrite existing files without further warning) (y/N):${txtrst}"
        read -e ANSWER
        if [ "$ANSWER" == "y" -o "$ANSWER" == "Y" ]; then
            echo "Ok, proceeding...."
        else
            exit 1
        fi
    else
        mkdir $INST_DIR/sm-01_mod1
        mkdir $INST_DIR/sm-01_mod1/bin
        mkdir $INST_DIR/sm-01_mod1/data
        mkdir $INST_DIR/sm-01_mod1/doc
        mkdir $INST_DIR/sm-01_mod1/doc/log
        mkdir $INST_DIR/sm-01_mod1/config

        mkdir $INST_DIR/sm-01_mod2
        mkdir $INST_DIR/sm-01_mod2/config
        mkdir $INST_DIR/sm-01_mod2/data
        mkdir $INST_DIR/sm-01_mod2/data/lib
        mkdir $INST_DIR/sm-01_mod2/doc
        mkdir $INST_DIR/sm-01_mod2/src
        mkdir $INST_DIR/sm-01_mod2/src/test
    fi
else
    echo -n "${txtbld}WARNING:${txtrst} '$INST_DIR' does not exist! Should it be created (Y/n)? "
    read -e ANSWER
    if [ "$ANSWER" == "" -o "$ANSWER" == "y" -o "$ANSWER" == "Y" ]; then
        mkdir $INST_DIR
        mkdir $INST_DIR/sm-01_mod1
        mkdir $INST_DIR/sm-01_mod1/bin
        mkdir $INST_DIR/sm-01_mod1/data
        mkdir $INST_DIR/sm-01_mod1/doc
        mkdir $INST_DIR/sm-01_mod1/doc/log
        mkdir $INST_DIR/sm-01_mod1/config

        mkdir $INST_DIR/sm-01_mod2
        mkdir $INST_DIR/sm-01_mod2/config
        mkdir $INST_DIR/sm-01_mod2/data
        mkdir $INST_DIR/sm-01_mod2/data/lib
        mkdir $INST_DIR/sm-01_mod2/doc
        mkdir $INST_DIR/sm-01_mod2/src
        mkdir $INST_DIR/sm-01_mod2/test
    else
        echo "Cannot continue without installation directory - stopped."
        exit 1
    fi
    echo "  "
fi

echo " "
echo "${txtbld}INFO:${txtrst} SM-01 Module 1 will be installed to $INST_DIR/sm-01_mod1"
echo "      SM-01 Module 2 will be installed to $INST_DIR/sm-01_mod2"
echo " "


# checking CPL
echo -n "Input CPL installation directory : "
read -e CPL_DIR

if [ -d ${CPL_DIR} ]; then
    echo " "
else
    echo "No CPL found! Cannot continue."
    exit 0
fi

# checking for cpl files
if [ ! -e ${CPL_DIR}/lib/libcplcore.a -o \
     ! -e ${CPL_DIR}/lib/libcpldfs.a -o \
     ! -e  ${CPL_DIR}/lib/libcext.a -o \
     ! -e  ${CPL_DIR}/lib/libcpldrs.a -o \
     ! -e  ${CPL_DIR}/lib/libcplui.a -o \
     ! -e  ${CPL_DIR}/include/cpl.h ]; then
    echo "Some CPL files not found! Please restart installation script and provide the correct CPL path."
    exit 0
fi
LD_LIBRARY_PATH=${CPL_DIR}/lib/
# +++++++++++++++++++++++++++++++++  Module 1 ++++++++++++++++++++++++++++++++++
cd sm-01_mod1
# start copying files
echo "${txtbld}Starting installation of Module 1......${txtrst}"
echo -n "Copying config  files to <MOD1-DIR>/config"
cp -r config/* $INST_DIR/sm-01_mod1/config
echo " - ${txtbld}${txtgrn}done${txtrst}"

echo -n "Copying data    files to <MOD1-DIR>/data (may take a while) "
cp data/* $INST_DIR/sm-01_mod1/data 2>/dev/null
echo " - ${txtbld}${txtgrn}done${txtrst}"

echo -n "Copying documentation to <MOD1-DIR>/doc  "
cp doc/*.pdf $INST_DIR/sm-01_mod1/doc 2>/dev/null
echo " - ${txtbld}${txtgrn}done${txtrst}"

# start compilation
file="./sm-01_mod1_gcc.log"
if [ -e $file ]; then
    mv sm-01_mod1_gcc.log sm-01_mod1_gcc.log.bak
    echo "Copying sm-01_mod1_gcc.log to sm-01_mod1_gcc.log.bak - ${txtbld}${txtgrn}done${txtrst}"
fi

echo -n "Compiling sources (compiler output --> sm-01_mod1_gcc.log) "

touch sm-01_mod1_gcc.log
./bootstrap >> sm-01_mod1_gcc.log
./configure --prefix=$INST_DIR/sm-01_mod1/ --with-cpl=${CPL_DIR} CFLAGS="-std=c99 -Wno-error" >> sm-01_mod1_gcc.log
make  >> sm-01_mod1_gcc.log
make install  >> sm-01_mod1_gcc.log
make clean  >> sm-01_mod1_gcc.log
cp sm-01_mod1_gcc.log $INST_DIR/sm-01_mod1/doc/log/
# checking success
file="$INST_DIR/sm-01_mod1/bin/create_speclib"
if [ -e $file ]; then
    echo "    - ${txtbld}${txtgrn}done${txtrst}"
    echo "Installation of Module 1 (to '$INST_DIR/sm-01_mod1/') - ${txtbld}${txtgrn}SUCCESSFUL${txtrst}"
else
    echo "    ${txtbld}${txtred}- NOT successfull!${txtrst}"
    echo "${txtbld}${txtred}Check $PWD/sm-01_mod1_gcc.log for more details!${txtrst}"
    echo "Installation of Module 1 (to '$INST_DIR/sm-01_mod1/') - ${txtbld}${txtred}NOT SUCCESSFUL${txtrst}"
fi
echo " "
cd ..

#===========================================================================================
#
#
#
# Installing LNFL / LBLRTM / AER package
#
# This script is optimised for LNFL V2.6 / LBLRTM 12.2, and AER_V3.2.
# There might be changes in the way of compilation of the external software.
# Modify this section if the LNFL / LBLRTM package is updated.
# Read the README files in the tree of the sources.
#
#
#
#===========================================================================================
echo " "
echo "${txtbld} Installing LNFL/LBLRTM and AER line list........${txtrst}"
echo -n "Unpacking......"
cd $sourcedir/sm-01_mod1/third_party_code/
tar zxf lnfl_lblrtm_aer.tar.gz
echo " - ${txtbld}${txtgrn}SUCCESSFUL${txtrst}"

echo -n "Compiling and installing lnfl package"
file="$INST_DIR/sm-01_mod1/doc/log/lnfl_compilation.log"
touch $file
cd $sourcedir/sm-01_mod1/third_party_code/lnfl/build
make -f make_lnfl linuxGNUsgl 1>$file 2>$file
cd ..
if [ ! -e ./lnfl_v2.6_linux_gnu_sgl ];
then
    echo " - ${txtbld}${txtred}- NOT successful!${txtrst}"
    echo " Please check the file '$file' for more details and try to install"
    echo " LNFL manually (see SM-01 User Manual)."
else
    echo "   - ${txtbld}${txtgrn}SUCCESSFUL${txtrst}"
    cp ./lnfl_v2.6_linux_gnu_sgl $INST_DIR/sm-01_mod1/bin/
    cd $INST_DIR/sm-01_mod1/bin/
    if [ -e lnfl ];then
        rm lnfl
        ln -s ./lnfl_v2.6_linux_gnu_sgl lnfl
    else
        ln -s ./lnfl_v2.6_linux_gnu_sgl lnfl
    fi
fi


cd $sourcedir/sm-01_mod1/third_party_code/lblrtm/build/
echo -n "Compiling and installing lblrtm package"
file="$INST_DIR/sm-01_mod1/doc/log/lblrtm_compilation.log"
touch $file
make -f make_lblrtm linuxGNUsgl 1>$file 2>$file
cd ..

if [ ! -e ./lblrtm_v12.2_linux_gnu_sgl ];
then
    echo " - ${txtbld}${txtred}- NOT successful!${txtrst}"
    echo " Please check the file '$file' for more details and try to install"
    echo " LNFL manually (see SM-01 User Manual)."
else
    echo " - ${txtbld}${txtgrn}SUCCESSFUL${txtrst}"
    cp ./lblrtm_v12.2_linux_gnu_sgl $INST_DIR/sm-01_mod1/bin/
    cd $INST_DIR/sm-01_mod1/bin/
    if [ -e lblrtm ];then
        rm lblrtm
        ln -s ./lblrtm_v12.2_linux_gnu_sgl lblrtm
    else
        ln -s ./lblrtm_v12.2_linux_gnu_sgl lblrtm
    fi
fi

echo -n "Preparing and installing aer line list"
cd $sourcedir/sm-01_mod1/third_party_code/aer_v_3.2/line_file
if [ ! -e ./aer_v_3.2 ];
then
    echo " - ${txtbld}${txtred}- NOT successful!${txtrst}"
    echo " Please check the file '$sourcedir/third_party_code/aer_v_3.2/line_file/aer_v_3.2' and try to copy it to"
    echo "  $INST_DIR/sm-01_mod1/data/hitran/ manually (see SM-01 User Manual)."
else
    cp ./aer_v_3.2 $INST_DIR/sm-01_mod1/data/
    echo "  - ${txtbld}${txtgrn}SUCCESSFUL${txtrst}"
fi
cd $sourcedir/

#===========================================================================================
#
#
#
#     END OF SECTION TO BE MODIFIED
#
#
#
#===========================================================================================

# +++++++++++++++++++++++++++++++++  Module 2 ++++++++++++++++++++++++++++++++++
cd sm-01_mod2
# start copying files
echo "${txtbld}Starting installation of Module 1a and Module 2......${txtrst}"
echo -n "Copying config files  to $INST_DIR/sm-01_mod2/config"
cp -r config/* $INST_DIR/sm-01_mod2/config
echo " - ${txtbld}${txtgrn}done${txtrst}"

echo -n "Copying data   files  to $INST_DIR/sm-01_mod2/data "
cp data/* $INST_DIR/sm-01_mod2/data 2>/dev/null
echo " - ${txtbld}${txtgrn}done${txtrst}"

echo -n "Copying documentation to $INST_DIR/sm-01_mod2/doc"
cp -r doc/*.pdf $INST_DIR/sm-01_mod2/doc
echo " - ${txtbld}${txtgrn}done${txtrst}"

echo -n "Copying test files    to $INST_DIR/sm-01_mod2/test"
cp -r test $INST_DIR/sm-01_mod2/
echo " - ${txtbld}${txtgrn}done${txtrst}"

echo -n "Creating dir             $INST_DIR/sm-01_mod2/output"
mkdir $INST_DIR/sm-01_mod2/output
echo " - ${txtbld}${txtgrn}done${txtrst}"

# start compilation
file="./sm-01_mod2_gcc.log"
if [ -e $file ]; then
    mv sm-01_mod2_gcc.log sm-01_mod2_gcc.log.bak
    echo "Copying sm-01_mod2_gcc.log to sm-01_mod2_gcc.log.bak            - ${txtbld}${txtgrn}done${txtrst}"
fi

echo -n "Compiling sources (compiler output --> sm-01_mod2_gcc.log) "

touch sm-01_mod2_gcc.log
./bootstrap >> sm-01_mod2_gcc.log
./configure --prefix=$INST_DIR/sm-01_mod2/ --with-cpl=${CPL_DIR} CFLAGS="-std=c99 -Wno-error -Wno-format -Wno-format-y2k" >> sm-01_mod2_gcc.log
make  >> sm-01_mod2_gcc.log
make install  >> sm-01_mod2_gcc.log
#cp test/test_skyemcomp $INST_DIR/sm-01_mod2/test/ 2>/dev/null
make clean  >> sm-01_mod2_gcc.log

# checking success
file="$INST_DIR/sm-01_mod2/bin/preplinetrans"
if [ -e $file ]; then
    echo "    - ${txtbld}${txtgrn}done${txtrst}"
    echo "Installation of Module 2 (to '$INST_DIR/sm-01_mod2/')  - ${txtbld}${txtgrn}SUCCESSFUL${txtrst}"
else
    echo "    ${txtbld}${txtred}- NOT successfull!${txtrst}"
    echo "${txtbld}${txtred}Check $PWD/sm-01_mod2_gcc.log for more details!${txtrst}"
    echo "Installation of Module 2 (to '$INST_DIR/sm-01_mod2/')  - ${txtbld}${txtred}NOT SUCCESSFUL${txtrst}"
fi

cd ..
# +++++++++++++++++++++++++++++++++  sky library ++++++++++++++++++++++++++++++++++
#echo " "
#cd lib/
#echo -n "Moving sky library to <MOD2-DIR>/data/lib/"
#mv LBL.R60k/*.fits $INST_DIR/sm-01_mod2/data/lib/ 2>/dev/null
#mv LBL.R300k/*.fits $INST_DIR/sm-01_mod2/data/lib/ 2>/dev/null
#mv LBL.R1000k/*.fits $INST_DIR/sm-01_mod2/data/lib/ 2>/dev/null
#mv LBL.PWV/*.fits $INST_DIR/sm-01_mod2/data/lib/ 2>/dev/null
#mv L-files/*.fits $INST_DIR/sm-01_mod2/data/lib/ 2>/dev/null
#echo " - ${txtbld}${txtgrn}done${txtrst}"
#
#cd ..

echo " "
echo "${txtbld}NOTE:${txtrst} Don't forget to check the directory $INST_DIR/sm-01_mod1/bin/"
echo "      for the correct binaries of lnfl and lblrtm (see SM-01 User Manual for "
echo "      more information on compilation and required naming)"
echo " "
echo "      Please also make sure that the path ${CPL_DIR}/lib is added to "
echo "      the LD_LIBRARY_PATH environment variable."
echo " "
echo "Installation finished. "

exit 0

