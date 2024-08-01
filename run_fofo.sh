#!/usr/bin/env bash

###############################################################
### run_fofo                                version 240724   ##
### (C) 2024 Adams Vallejos                                  ##
###############################################################

BASEDIR="/PATH/TO/resampling"
METHOD="B"

OUTDIR="$BASEDIR/output/${LIGHT}/${METHOD}"

DARK='dark-reference'
LIGHT="${@: -1}"

FOFODIR="analysis_fofo_0"

if [ ! $(which phenix.fobs_minus_fobs_map) ]
then
        echo "Source phenix and relaunch"
        exit 1
else
        echo "Calculating FoFo maps in ${OUTDIR}"
fi


for MTZIN in `find ${OUTDIR} -name '*_uniq.mtz'|sort`
do
        MTZDIR=$(dirname $MTZIN)
        FOFODIR=${MTZDIR}/${FOFODIR}
        MTZFILENAME=$(basename $MTZIN .mtz)
        # MTZFILENAME=${MTZIN:0:-4}
        echo Working on ${MTZIN}
        cd ${MTZDIR}
        mkdir -p ${FOFODIR}
        PREFIX=${FOFODIR}/${MTZFILENAME}
        NEGATIVEMTZ="${MTZIN//$LIGHT/$DARK}"
        if ! phenix.fobs_minus_fobs_map \
                f_obs_1_file=${MTZIN} \
                f_obs_1_label='F,SIGF' \
                f_obs_2_file=${NEGATIVEMTZ} \
                f_obs_2_label='F,SIGF' \
                high_resolution=2.0 \
                low_resolution=38.9 \
                phase_source="${NEGATIVEMTZ:0:-4}"_001.pdb \
                output_file=${PREFIX}_fofo.mtz \
                output_dir=$(dirname $MTZIN) \
                file_name_prefix="${PREFIX}"_fofo \
                >> ${PREFIX}_fofo_map_coeffs.log
        then
                echo "Failed to calculate FoFo"
                exit 1
        fi
        
        pdb_file="${PREFIX}_cad_uniq_001.pdb"
        mtz_file="${PREFIX}_fofo.mtz"
        cat << EOF > "${PREFIX}_fofo_map_coeffs.py"
COORDS = handle_read_draw_molecule("$pdb_file")
MAP = make_and_draw_map("$mtz_file","FoFo","PHFc","",0,1)
export_map(MAP,"${PREFIX}_fofo_map_coeffs.ccp4")
coot_no_state_real_exit(0)
EOF
        # if
        #         ! coot --no-graphics --script "${PREFIX}_fofo_map_coeffs.py" > "${PREFIX}_fofo_map_coeffs.ccp4.log"
        # then
        if ! phenix.mtz2map\
            mtz_file=${PREFIX}_fofo.mtz \
            labels='FoFo,PHFc' \
            scale=sigma \
            output.prefix=${PREFIX}_fofo_map_coeffs \
            output.format=ccp4 \
            > ${PREFIX}_fofo_map_coeffs.log

        then
                "Failed to generate map"
                exit 1
        fi
        mv "${PREFIX}_fofo_map_coeffs_1.ccp4" "${PREFIX}_fofo_map_coeffs.ccp4" 

    # # #   --- modified from mapTool: processMaps.sh ---
    # 1. EXTEND INITIAL MAP TO COVER FULL CELL (necessary for mapprot to run)
    #----------------------------------------
    # Run MAPMASK from CCP4 suite
        echo "Extending maps"
        if ! mapmask \
                MAPIN ${PREFIX}_fofo_map_coeffs.ccp4 \
                MAPOUT ${PREFIX}_fofo_map_coeffs_fullCell.ccp4 \
                <<-EOF > ${PREFIX}_fofo_map_coeffs_fullCell.log
                XYZLIM CELL
EOF
        then
                echo "Failed to extend map!"
                exit 1
        fi

        # 2. TRANSLATE MAP TO CARTESIAN COORDINATES
        #    Set CELL WORK / GRID WORK / XYZLIM to cover the volume of interest
        #-----------------------------------------
        # CELL WORK: A B C 90 90 90
        #   ->  A B C sides of cell in Ångström (= size of pdb + 5Å border)
        #   angles 90 degree 
        # GRID WORK: 4*A 4*B 4*C 
        #   -> 0.25 Å distance between grid points
        # XYZLIM: 4*xmin 4*xmax 4*ymin 4*ymax 4*zmin 4*zmax
        #-----------------------------------------
        # The numbers below are obtained by pdb_size in xshuffle.py
        echo "Translating map to cartesian cordinates"
        INMAP=${OUTMAP}
        OUTMAP+="_cartesian"
        # 0.016 
        # CELL WORK 43 55 78 90 90 90
        # GRID WORK 172 220 312
        # XYZLIM -144 28 -180 40 -136 176
        # 0.760 & 1725
        # CELL WORK 44 55 78 90 90 90
        # GRID WORK 176 220 312
        # XYZLIM -148 28 -180 40 -136 176
        if ! maprot \
                MAPIN ${PREFIX}_fofo_map_coeffs_fullCell.ccp4 \
                WRKOUT ${PREFIX}_fofo_map_coeffs_fullCell_cartesian.map \
                <<-EOF > ${PREFIX}_fofo_map_coeffs_fullCell_cartesian.log
                        MODE FROM
                        CELL WORK 98 72 74 90 90 90
                        GRID WORK 392 288 296
                        XYZLIM -192 200 -132 156 -36 260
                        SYMM WORK 1
                        AVER
                        ROTA POLAR 0 0 0
                        TRANS 0 0 0
EOF
        then
                echo "Failed to translate map!"
                exit 1
        fi

        # 3. CONVERT MAPS TO HDF5 (for matlab)
        #----------------------------------------------------------------------------
        echo "Converting map to h5"
        if ! maptool-resampling/map_to_h5.py \
                ${PREFIX}_fofo_map_coeffs_fullCell_cartesian.map\
                ${PREFIX}_fofo_map_coeffs_fullCell_cartesian.h5
        then
                echo "Failed to convert map!"
                exit 1
        fi

        # 4. EXTRACT INFORMATION ON MAP DIMESIONS
        #    cell dimensions, grid sampling, start and stop, axis order
        #    and sigma (=RMSD) from new and original map
        #    (Sigma of the original map is used in the matlab calculations since the 
        #    sigma value for the cartesian map is slightly changed in the conversion.)
        #----------------------------------------------------------------------------
        # Example of data in mapdump output:
        # Cell dimensions .................................   42.0000    54.0000    77.0000    90.0000    90.0000    90.0000
        # Grid sampling on x, y, z ........................  168  216  308
        # Start and stop points on columns, rows, sections  -136  172 -144   24 -176 
        # Fast, medium, slow axes  3  1  2
        # Rms deviation from mean density .................     0.01195
        #----------------------------------------------------------------------------
        echo "Extracting map dimensions"
        if ! mapdump \
                MAPIN ${PREFIX}_fofo_map_coeffs_fullCell_cartesian.map \
                <<-EOF > ${PREFIX}_fofo_map_coeffs_fullCell_cartesian_header.txt
EOF
        then
                echo "Failed to extract map dimensions!"
                exit 1
        fi

        IN=${IN}
        OUT="${IN}_original_header"
        if ! mapdump \
                MAPIN ${PREFIX}_fofo_map_coeffs.ccp4 \
                <<-EOF > ${PREFIX}_fofo_map_coeffs_original_header.txt
EOF
        then
                echo "Failed to extract original dimensions!"
                exit 1
        fi

        grep 'Cell dimensions' ${PREFIX}_fofo_map_coeffs_fullCell_cartesian_header.txt | awk '{ print $4,$5,$6,$7,$8,$9}' > ${PREFIX}_fofo_map_coeffs_XYZinfo.dat
        grep 'Grid sampling' ${PREFIX}_fofo_map_coeffs_fullCell_cartesian_header.txt | awk '{ print $8,$9,$10}' >> ${PREFIX}_fofo_map_coeffs_XYZinfo.dat
        grep 'Start and stop' ${PREFIX}_fofo_map_coeffs_fullCell_cartesian_header.txt | awk '{ print $9,$10,$11,$12,$13,$14}' >> ${PREFIX}_fofo_map_coeffs_XYZinfo.dat
        grep '^ Fast' ${PREFIX}_fofo_map_coeffs_fullCell_cartesian_header.txt | awk '{ print $5,$6,$7}' >> ${PREFIX}_fofo_map_coeffs_XYZinfo.dat
        grep 'Rms deviation' ${PREFIX}_fofo_map_coeffs_fullCell_cartesian_header.txt | awk '{ print $7}'  >> ${PREFIX}_fofo_map_coeffs_XYZinfo.dat
        grep 'Rms deviation' ${PREFIX}_fofo_map_coeffs_original_header.txt | awk '{ print $7}'  >> ${PREFIX}_fofo_map_coeffs_XYZinfo.dat

done

echo Task Finished
