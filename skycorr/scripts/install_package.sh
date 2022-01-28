#!/bin/bash
#
#   This file is part of the SKYCORR software package.
#   Copyright (C) 2009-2013 European Southern Observatory
#
#   This programme is free software; you can redistribute it and/or modify
#   it under the terms of the GNU General Public License as published by
#   the Free Software Foundation; either version 2 of the License, or
#   (at your option) any later version.
#
#   This programme is distributed in the hope that it will be useful,
#   but WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
#   GNU General Public License for more details.
#
#   You should have received a copy of the GNU General Public License
#   along with this programme. If not, see <http://www.gnu.org/licenses/>.

set -e
sourcedir=$PWD
tmpdir="$sourcedir/$(mktemp -d skycorr_XXXXXX || exit 1)"
trap "rm -rf $tmpdir" EXIT
set +e

# defining console output text properties
txtbld=$(tput bold)       # Bold
txtrst=$(tput sgr0)       # Text reset
txtred=$(tput setaf 1)    # Red
txtgrn=$(tput setaf 2)    # Green
txtylw=$(tput setaf 5)    # magenta

BINTAR=$(ls skycorr_binaries_*.tar.gz)
SRCTAR=$(ls skycorr-*.tar.gz)
# contains directory the installer was run from, not PWD when run from self
# extracting tarball
BASEDIR=${BASEDIR:-"$PWD"}

if [ $# -eq 1 ]; then
    INST_DIR="$1"
    ANSWER="y"
elif [ $# -eq 3 ]; then
    INST_DIR="$1"
    BINTAR="$2"
    SRCTAR="$3"
    ANSWER="y"
fi
echo "Installing $BINTAR and $SRCTAR into $INST_DIR"

if [ ! -r "$BINTAR" ]; then
  echo "Binary tarball $BINTAR not found, please execute script in same " \
       "folder the package was extracted to."
  exit 1
fi
if [ ! -r "$SRCTAR" ]; then
  echo "Source tarball $SRCTAR not found, please execute script in same " \
       "folder the package was extracted to."
  exit 1
fi

set -e

# preparation
echo "${txtbld} "
echo "*************************************************************************"
echo "*  Welcome to the SKYCORR setup routine                                 *"
echo "*                                                                       *"
echo "*                                                                       *"
echo "*                                                                       *"
echo "*  Austrian ESO in-kind Team Innsbruck                                  *"
echo "*  W. Kausch, S. Noll, S. Kimeswenger (head),  M. Barden, C. Szyszka,   *"
echo "*  Amy M. Jones                                                         *"
echo "*  Institute of Astro- and Particle Physics                             *"
echo "*  University of Innsbruck                                              *"
echo "*  Technikerstr. 25                                                     *"
echo "*  A-6020 Innsbruck / Austria                                           *"
echo "*                                                                       *"
echo "*************************************************************************"
echo "${txtrst} "

while [ ! $ANSWER != "n" ] || [ -z "$INST_DIR" ]
do
    echo "The software package 'skycorr' is aimed at correcting the sky emission "
    echo "signature from scientific spectra. It is an enhanced development of the"
    echo "sky correction algorithm by R. Davies (2007, MNRAS, 375, 1099)."
    echo " "
    echo "This script installs the skycorr package to a definable installation"
    echo "directory <INST_DIR> creating the following directory structure:"
    echo " "
    echo "    <INST_DIR>"
    echo "            ├── bin"
    echo "            ├── config"
    echo "            ├── doc"
    echo "            ├── examples"
    echo "            │   ├── config"
    echo "            │   ├── data"
    echo "            │   └── sysdata -> <INST_DIR>/sysdata/"
    echo "            ├── gui"
    echo "            ├── output"
    echo "            ├── reflex"
    echo "            └── sysdata"
    echo " "
    echo " "

    echo -n "${txtbld}Step 1:${txtrst} Choose root installation directory <INST_DIR>: "
    read -e INST_DIR


    LEN=$(echo ${#INST_DIR})
    export nogo="~"

    # checking write permissions
    if [  "${INST_DIR}" = ""  ];then
        echo "Cannot install in directory with empty name: '${INST_DIR}'"
        echo "Exiting."
        exit
    fi

    if [ "${INST_DIR:0:1}" = "${nogo}" ]; then
        if [  "${HOME}" = ""  ];then
            echo "'~' not defined. Please specify path without '~'."
            echo "Exiting."
            exit
        fi
        echo "${txtbld}${txtylw}WARNING: ${txtrst} No '~' in directory path allowed."
        echo "          Should it be replaced by your home directory '$HOME'?"
        echo "          New installation path would be: $HOME/${INST_DIR:1:$LEN}"
        export INST_DIR=$HOME/${INST_DIR:1:$LEN}
    fi
    echo -n "Is this OK [Y/n]? "

        read -e ANSWER
    if [ "$ANSWER" == "" -o "$ANSWER" == "y" -o "$ANSWER" == "Y" ]; then
         ANSWER="y"
        if [ $LEN -gt 40 ]; then
            echo "${txtbld}${txtylw}WARNING: ${txtrst} Using '$INST_DIR' "
            echo "          probably exceeds the 80 character limit of path+filenames "
            echo "          introduced by the Fortran codes LBLRTM."
            echo "          May cause problems using the radiative transfer codes!"
            echo -n "          Proceed with this installation directory (y/N)? "
            read -e ANSWER
            if [ "$ANSWER" == "y" -o "$ANSWER" == "Y" ]; then
                ANSWER="y"
                echo " "
            else
                ANSWER="n"
                echo " "
            fi
        fi
    else
        ANSWER="n"
    fi
done
if [ "${INST_DIR:0:1}" != "/" ]; then
  INST_DIR="$BASEDIR/$INST_DIR"
fi

# creating directories
if [ -d "$INST_DIR" ]; then
    echo "${txtbld}${txtred}Oops, it seems that the installation already exists...."
    echo -n "Proceed (will overwrite existing files without further warning) (y/N):${txtrst}"
    read -e ANSWER
    if [ "$ANSWER" == "y" -o "$ANSWER" == "Y" ]; then
        echo "Ok, proceeding...."
    else
        exit 1
    fi
fi

mkdir -p "$INST_DIR"

# start copying files
echo " "
echo "${txtbld}Step 2: Starting installation of skycorr dependencies......${txtrst}"
tar -C "$INST_DIR" -x -z -f "$BINTAR"

# revert opensuse weirdness
mv $INST_DIR/lib32/* $INST_DIR/lib 2>/dev/null || true
mv $INST_DIR/lib64/* $INST_DIR/lib 2>/dev/null || true
rmdir $INST_DIR/lib32 $INST_DIR/lib64 2>/dev/null || true
export LD_LIBRARY_PATH="$INST_DIR/lib:$LD_LIBRARY_PATH"
export DYLD_LIBRARY_PATH="$INST_DIR/lib:$DYLD_LIBRARY_PATH"

# start copying files
echo " "
echo "${txtbld}Step 3: Starting installation of skycorr${txtrst}"

tar -C "$tmpdir" --strip-components=1 -x -z -f $SRCTAR
cd "$tmpdir"

# start compilation
file="$INST_DIR/doc/log/skycorr_gcc.log"
if [ -e $file ]; then
    cp $file $file.bak
    echo "Copying $file to $file.bak - ${txtbld}${txtgrn}done${txtrst}"
fi

echo -n "Compiling sources (see $file for compiler output) "
mkdir -p $INST_DIR/doc/log

set +e
(
set -ex
./configure --prefix=$INST_DIR --with-cpl=${INST_DIR} --libdir=$INST_DIR/lib
make
make install
) >> $file 2>>$file
if [ $? -ne 0 ]; then
    echo "    ${txtbld}${txtred}- NOT successful!${txtrst}"
    echo
    echo "===================="
    echo "Last lines of $file:"
    echo "===================="
    echo
    tail -n 20 "$file"
    echo "==========" >> $file
    echo "config.log" >> $file
    echo "==========" >> $file
    # add some extra info
    cat config.log >> $file
    echo "Installed binary: $BINTAR" >> $file
    echo "Installed source: $SRCTAR" >> $file
    echo
    echo "Aborting."
    exit 1;
fi
echo "- ${txtbld}${txtgrn}done${txtrst}"

# Finishing
if [ "$(uname -a | grep -i darwin)" ]; then
  LIBVAR="DYLD_LIBRARY_PATH"
else
  LIBVAR="LD_LIBRARY_PATH"
fi

# add environment setup wrappers for user executables
for exe in skycorr; do
mv "$INST_DIR/bin/${exe}" "$INST_DIR/bin/${exe}.real"
cat << EOF >> "$INST_DIR/bin/${exe}"
#!/bin/sh
export $LIBVAR="$INST_DIR/lib\${$LIBVAR:+":\$$LIBVAR"}"
exec $INST_DIR/bin/${exe}.real "\$@"
EOF
chmod u+x "$INST_DIR/bin/${exe}"
done

echo "Skycorr has been successfully installed into $INST_DIR"
echo
echo "Run an example file using this command:"
echo
echo "$INST_DIR/bin/skycorr $INST_DIR/examples/config/sctest_sinfo_H.par"
echo
