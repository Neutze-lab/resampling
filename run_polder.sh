#!/usr/bin/env bash

###############################################################
### run_polder                              version 240724   ##
### (C) 2024 Adams Vallejos                                  ##
###############################################################

BASEDIR="/PATH/TO/resampling"
METHOD="B"

SELECTION='chain a and (resid 400 or resid 401)'

DATASET="${@: -1}"

SOURCEDIR="$BASEDIR/output/${DATASET}/${METHOD}"
cd $SOURCEDIR

if [ ! $(which phenix.polder) ]
then
        printf "Source phenix and relaunch\n"
        exit 1
else
        printf "Calculating Polder-maps in:\n\t ${SOURCEDIR}\n"
fi


for MTZIN in `find * -name '*_001.mtz'|sort`
do
        printf "Currently at ${MTZIN}\n"
        DIR_MTZ=$(dirname $MTZIN)
        DIR_POLDER=$DIR_MTZ/analysis_polder_2
        mkdir -p ${DIR_POLDER}
        MTZ_NAME=$(basename $MTZIN)
        OUT_PREFIX=$DIR_POLDER/${MTZ_NAME:0:-4}
        MAPLABELS='F-obs-filtered,SIGF-obs-filtered'
        # All the atoms
        SELECTION=''
        
        # With only chlorides and sulphurs
        printf "\r\033[K\t1/6 Calculating Polder-map\r"
	# echo File is "${MTZIN:0:-4}.pdb"
	# exit
        if ! phenix.polder \
                model_file_name=${MTZIN:0:-4}.pdb \
                solvent_exclusion_mask_selection="${SELECTION}" \
                reflection_file_name=${MTZIN} \
		data_labels=${MAPLABELS} \
                r_free_flags_labels=R-free-flags \
                high_resolution=2.0 \
                low_resolution=38.9 \
                output_file_name_prefix=${OUT_PREFIX} \
		output.overwrite=True \
                > ${OUT_PREFIX}_polder_map_coeffs.log
        then
                printf "\tFailed to calculate Polder map.\n"
                exit 1
        fi

	        printf "\r\033[K\t2/6 Converting to CCP4 format\r"
        if ! phenix.mtz2map \
                mtz_file=${OUT_PREFIX}_polder_map_coeffs.mtz \
                labels='mFo-DFc_polder,PHImFo-DFc_polder' \
                scale=sigma \
                output.prefix=${OUT_PREFIX}_polder_map_coeffs \
                output.format=ccp4 \
                > ${OUT_PREFIX}_polder_map_coeffs.ccp4.log

        then
                printf "\tFailed to generate map.\n"
                exit 1
        fi
        mv "${OUT_PREFIX}_polder_map_coeffs_1.ccp4" "${OUT_PREFIX}_polder_map_coeffs.ccp4"

    # # #   --- modified from mapTool: processMaps.sh ---
    # 1. EXTEND INITIAL MAP TO COVER FULL CELL (necessary for mapprot to run)
    #----------------------------------------
    # Run MAPMASK from CCP4 suite
        printf "\r\033[K\t3/6 Extending map\r"
        if ! mapmask \
                MAPIN ${OUT_PREFIX}_polder_map_coeffs.ccp4 \
                MAPOUT ${OUT_PREFIX}_polder_map_coeffs_fullCell.ccp4 \
                <<-EOF > ${OUT_PREFIX}_polder_map_coeffs_fullCell.log
                XYZLIM CELL
EOF
        then
                printf "\tFailed to extend map!\n"
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
        printf "\r\033[K\t4/6 Translating map to cartesian cordinates\r"
        INMAP=${OUTMAP}
        OUTMAP+="_cartesian"
        if ! maprot \
                MAPIN ${OUT_PREFIX}_polder_map_coeffs_fullCell.ccp4 \
                WRKOUT ${OUT_PREFIX}_polder_map_coeffs_fullCell_cartesian.map \
                <<-EOF > ${OUT_PREFIX}_polder_map_coeffs_fullCell_cartesian.log
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
                printf "\tFailed to translate map!\n"
                exit 1
        fi

        # 3. CONVERT MAPS TO HDF5 (for matlab)
        #----------------------------------------------------------------------------
        printf "\r\033[K\t5/6 Converting map to h5\r"
        if ! ${BASEDIR}/maptool-cco/map_to_h5.py \
                ${OUT_PREFIX}_polder_map_coeffs_fullCell_cartesian.map\
                ${OUT_PREFIX}_polder_map_coeffs_fullCell_cartesian.h5
        then
                printf "\tFailed to convert map!\n"
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
        printf "\r\033[K\t6/6Extracting map dimensions\r"
        if ! mapdump \
                MAPIN ${OUT_PREFIX}_polder_map_coeffs_fullCell_cartesian.map \
                <<-EOF > ${OUT_PREFIX}_polder_map_coeffs_fullCell_cartesian_header.txt
EOF
        then
                printf "\tFailed to extract map dimensions!\n"
                exit 1
        fi

        IN=${IN}
        OUT="${IN}_original_header"
        if ! mapdump \
                MAPIN ${OUT_PREFIX}_polder_map_coeffs.ccp4 \
                <<-EOF > ${OUT_PREFIX}_polder_map_coeffs_original_header.txt
EOF
        then
                printf "\tFailed to extract original dimensions!\n"
                exit 1
        fi

        grep 'Cell dimensions' ${OUT_PREFIX}_polder_map_coeffs_fullCell_cartesian_header.txt | awk '{ print $4,$5,$6,$7,$8,$9}' > ${OUT_PREFIX}_polder_map_coeffs_XYZinfo.dat
        grep 'Grid sampling' ${OUT_PREFIX}_polder_map_coeffs_fullCell_cartesian_header.txt | awk '{ print $8,$9,$10}' >> ${OUT_PREFIX}_polder_map_coeffs_XYZinfo.dat
        grep 'Start and stop' ${OUT_PREFIX}_polder_map_coeffs_fullCell_cartesian_header.txt | awk '{ print $9,$10,$11,$12,$13,$14}' >> ${OUT_PREFIX}_polder_map_coeffs_XYZinfo.dat
        grep '^ Fast' ${OUT_PREFIX}_polder_map_coeffs_fullCell_cartesian_header.txt | awk '{ print $5,$6,$7}' >> ${OUT_PREFIX}_polder_map_coeffs_XYZinfo.dat
        grep 'Rms deviation' ${OUT_PREFIX}_polder_map_coeffs_fullCell_cartesian_header.txt | awk '{ print $7}'  >> ${OUT_PREFIX}_polder_map_coeffs_XYZinfo.dat
        grep 'Rms deviation' ${OUT_PREFIX}_polder_map_coeffs_original_header.txt | awk '{ print $7}'  >> ${OUT_PREFIX}_polder_map_coeffs_XYZinfo.dat

done

printf "\r\033[K\tTask Finished\n"
