#!/bin/sh

touch "$1"

echo "INST_DIR=$2" >> "$1"
echo "INPUT_OBJECT_SPECTRUM=$3" >> "$1"
echo "INPUT_SKY_SPECTRUM=$4" >> "$1"
echo "OUTPUT_DIR=$5" >> "$1"
echo "OUTPUT_NAME=$6" >> "$1"
echo "COL_NAMES=$7" >> "$1"
echo "DEFAULT_ERROR=$8" >> "$1"
echo "WLG_TO_MICRON=$9" >> "$1"
echo "VAC_AIR=${10}" >> "$1"
echo "DATE_KEY=${11}" >> "$1"
echo "TIME_KEY=${12}" >> "$1"
echo "TELALT_KEY=${13}" >> "$1"
echo "LINETABNAME=${14}" >> "$1"
echo "VARDATNAME=${15}" >> "$1"
echo "SOLDATURL=${16}" >> "$1"
echo "SOLDATNAME=${17}" >> "$1"
echo "SOLFLUX=${18}" >> "$1"
echo "FWHM=${19}" >> "$1"
echo "VARFWHM=${20}" >> "$1"
echo "LTOL=${21}" >> "$1"
echo "MIN_LINE_DIST=${22}" >> "$1"
echo "FLUXLIM=${23}" >> "$1"
echo "FTOL=${24}" >> "$1"
echo "XTOL=${25}" >> "$1"
echo "WTOL=${26}" >> "$1"
echo "CHEBY_MAX=${27}" >> "$1"
echo "CHEBY_MIN=${28}" >> "$1"
echo "CHEBY_CONST=${29}" >> "$1"
echo "REBINTYPE=${30}" >> "$1"
echo "WEIGHTLIM=${31}" >> "$1"
echo "SIGLIM=${32}" >> "$1"
echo "FITLIM=${33}" >> "$1"
echo "PLOT_TYPE=N" >> "$1"

#$INST_DIR/reflex/create_sc_parfile.sh $parfile $INST_DIR $INPUT_OBJECT_SPECTRUM $INPUT_SKY_SPECTRUM $OUTPUT_DIR $OUTPUT_NAME "$COL_NAMES" $DEFAULT_ERROR $WLG_TO_MICRON $VAC_AIR "$DATE_KEY" "$TIME_KEY" "$TELALT_KEY" $LINETABNAME $VARDATNAME $SOLDATURL $SOLDATNAME $SOLFLUX $FWHM $VARFWHM $LTOL $MIN_LINE_DIST $FLUXLIM $FTOL $XTOL $WTOL $CHEBY_MAX $CHEBY_MIN $CHEBY_CONST $REBINTYPE $WEIGHTLIM $SIGLIM $FITLIM
