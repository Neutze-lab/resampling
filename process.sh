#!/usr/bin/env bash

###############################################################
### XPROCESS: process                       version 240705   ##
### (C) 2024 Adams Vallejos                                  ##
###############################################################

TITLE=""
PROJECT=""
PROTEIN=""
CELL="${1}"
SPACEGROUP="${2}"
POINTGROUP="${3}"
SAMPLENAME="${4}"
COORDINATES="${5}"
CIFLIB="${6}"
CELLWORK="${7}"
GRIDWORK="${8}"
XYZLIM="${9}"
OUT="$(dirname $SAMPLENAME)"
RFREE=$(dirname $(dirname $(dirname $OUT)))/rfree.mtz
NPROC=10
NCYCLE=20
PARAMETERS="phenix_por.params"
W="\033[0m"
R="\033[31m"
G="\033[32m"
# Y="\033[33m"
# #   --- TEST SECTION

# #   --- CRYSTFEL SECTION ###################################################
# #   --- MERGING REFLECTIONS ---
# # INFO: Merge reflections by scaling and post refinement.

PARTIALATOR=1

# If partialator is true, run process_hkl with split
SPLIT=$((PARTIALATOR^1))
FOM=1

OUT+="/merge"
if [ ${PARTIALATOR} = 1 ]
then
    # #   --- Run PARTIALATOR ---
    # # example: https://www.desy.de/~twhite/crystfel/manual-partialator.html
    printf "%s Merging intensities with partialator\n" "$(date +"%y-%m-%d %H:%M")"
    if ! partialator \
            --input="${SAMPLENAME}" \
            --output="${OUT}".hkl \
            --symmetry="${POINTGROUP}" \
            --iterations=3 \
            --model=unity \
            -j ${NPROC} \
            > "${OUT}".hkl.log 2>&1
    then
        printf "%bPartialator failed!%b\n" "${R}" "${W}"
        exit 1
    fi
else
    # #   --- Run PROCESS_HKL ---
    printf "%s Merging intensities with process_hkl\n" "$(date +"%y-%m-%d %H:%M")"
    if ! process_hkl \
            --input="${SAMPLENAME}" \
            --output="${OUT}".hkl \
            --scale \
            --min-res=3 \
            --symmetry="${POINTGROUP}" \
            > "${OUT}".hkl.log 2>&1
    then
        printf "%bprocess_hkl failed!%b\n" "${R}" "${W}"
        exit 1
    fi
fi


if [ ${SPLIT} = 1 ]
then
    # #   --- Split half datasets for figures of merit_
    # #   ---  Odd only
    # printf "%s Spliting reflections: Odd only\n" "$(date +"%y-%m-%d %H:%M")"
    printf "%s Spliting reflections\n" "$(date +"%y-%m-%d %H:%M")"
    if ! process_hkl \
        --input="${SAMPLENAME}" \
        --output="${OUT}".hkl1 \
        --odd-only \
        --scale \
        --min-res=3 \
        --max-adu=18500 \
        --symmetry="${POINTGROUP}" \
        > "${OUT}".hkl1.log 2>&1
    then
        printf "%bOdd splitting failed!%b\n" "${R}" "${W}"
        exit 1
    fi
    # #   ---  Even only
    # printf "%s Spliting reflections: Even only\n" "$(date +"%y-%m-%d %H:%M")"
    if ! process_hkl \
        --input="${SAMPLENAME}" \
        --output="${OUT}".hkl2 \
        --even-only \
        --scale \
        --min-res=3 \
        --max-adu=18500 \
        --symmetry="${POINTGROUP}" \
        > "${OUT}".hkl2.log 2>&1
    then
        printf "%bEven slitting failed!%b\n" "${R}" "${W}"
        exit 1
    fi
else
    printf "%s Skipping reflection splitting.\n" "$(date +"%y-%m-%d %H:%M")"
fi



if [ ${FOM} = 1 ]
then
    # #   --- Compute FOM5n 
    printf "%s Computing figures of merit\n" "$(date +"%y-%m-%d %H:%M")"
    if ! check_hkl "${OUT}".hkl \
            -y "${POINTGROUP}" \
            -p "${COORDINATES}" \
            --shell-file="${OUT}"_check.dat \
            --nshells=20 \
            > "${OUT}"_check.log 2>&1
    then
        printf "%bcheck_hkl failed!%b\n" "${R}" "${W}"
        exit 1
    fi

    # printf "%s Computing Rsplit\n" "$(date +"%y-%m-%d %H:%M")"
    # Options: R1I, R1F, R2, Rsplit, CC, CCstar, CCano, CRDano, Rano, 
    #           Rano/Rsplit, d1sig, d2sig
    for fom in 'ccstar' 'rsplit';
    # for fom  in {'CChalf', 'ccstar', 'rsplit'};
    do
        # printf "%s Computing CC*\n" "$(date +"%y-%m-%d %H:%M")"
        if ! compare_hkl "${OUT}".hkl1 "${OUT}".hkl2 \
                -y "${POINTGROUP}" \
                -p "${COORDINATES}" \
                --fom="$fom" \
                --shell-file="${OUT}"_"$fom".dat \
                --nshells=20 \
                > "${OUT}"_"$fom".log 2>&1
        then
            printf "%bcompare_hkl failed!%b\n" "${R}" "${W}"
            exit 1
        fi
        printf "%s\n" "$(grep 'Overall' "${OUT}"_"$fom".log)"
    done
    # printf "%s\n" "$(grep 'Overall' "${OUT}"_"$fom".log)"
fi


# #   --- Remove stream file, resampled crystals logged in yaml file

if ! rm "${SAMPLENAME}"
then
    printf "%bSample stream could not be removed%b\n" "${R}" "${W}"
    exit 1
fi


# #   --- CCP4/PHENIX SECTION ################################################
# #   --- Convert HKL to MTZ ------------------------------------------------
# # INFO: Convert a formatted reflection file to MTZ format.
# # Based on CrystFEL's 'create-mtz' , deprecated from version 0.10.1
# # example: ccp4/8.0/examples/unix/runnable/f2mtz.exam 
printf "%s Writing mtz\n" "$(date +"%y-%m-%d %H:%M")"
sed -n '/End\ of\ reflections/q;p' "${OUT}".hkl > "${OUT}".temp.hkl
# get_hkl -i file.hkl -p my.cell -o out.mtz --output-format=mtz
if ! f2mtz \
    HKLIN "${OUT}".temp.hkl \
    HKLOUT "${OUT}".mtz \
    <<-ENDf2mtz > "${OUT}".log
        TITLE ${TITLE}
        NAME PROJECT ${PROJECT} CRYSTAL ${PROTEIN} DATASET $(basename "${SAMPLENAME}")
        CELL ${CELL}
        SYMM ${SPACEGROUP}
        SKIP 3
        LABOUT H K L I SIGI
        CTYPE  H H H J     Q
        FORMAT '(3(F4.0,1X),F10.2,10X,F10.2)'
        END
ENDf2mtz
then
    printf "%bFailed to write mtz!%b\n" "${R}" "${W}"
    exit 1
fi
rm "${OUT}".temp.hkl


# # FIRST TIME DO uniqueify {-p fraction} sample_0.mtz sample_0_rFree.mtz
# # ^ default fraction is 0.05
if [ ! -f "$RFREE" ]
then
    printf "%s Generating FreeR_flag column\n" "$(date +"%y-%m-%d %H:%M")"
    if ! uniqueify -p 0.05 "${OUT}".mtz "${RFREE}"
    then
        printf "%bFailed to generate FreeR_flag column!%b\n" "${R}" "${W}"
        exit 1
    fi
    if ! mv -f "$(basename "$RFREE" .mtz)".log "$(dirname "${RFREE}")"/.
    then
        printf "%bFailed to move log file!%b\n" "${R}" "${W}"
    fi
fi

# #   --- run CTRUNCATE ---
# # INFO: Intensity to amplitude conversion and data statistics.
# # example: ccp4/8.0/examples/unix/runnable/ctruncate.exam
printf "%s Converting intensities\n" "$(date +"%y-%m-%d %H:%M")"

IN=${OUT}
OUT+="_trunc"
if ! ctruncate -stdin \
    <<-ENDctruncate > ${OUT}.log
        hklin ${IN}.mtz
        hklout ${OUT}.mtz
        colin /*/*/[I,SIGI]
ENDctruncate
then
    printf "%bFailed to convert intensities!%b\n" "${R}" "${W}"
    exit 1
fi


# #   --- Run CAD ---
# # INFO: Collect and sort crystallographic reflection data from several files,
# #       ^ to generate a single set.
# # example: ccp4/8.0/examples/unix/runnable/cad.exam
printf "%s Combining amplitudes\n" "$(date +"%y-%m-%d %H:%M")"
IN=${OUT}
OUT+="_cad"
if ! cad \
    hklin1 ${IN}.mtz \
    hklin2 "${RFREE}" \
    hklout ${OUT}.mtz \
    <<-ENDcad > ${OUT}.log
        TITLE ${TITLE}
        LABIN FILE 1 ALLIn
        LABIN FILE 2 E1=FreeR_flag
        END
ENDcad
then
    printf "%bFailed to combine amplitudes!%b\n" "${R}" "${W}"
    exit 1
fi


# #   --- Run UNIQUEIFY ---
# # INFO: Generate a unique list of reflections
# # example: ccp4/8.0/examples/unix/runnable/unique-free-R
printf "%s Listing unique reflections\n" "$(date +"%y-%m-%d %H:%M")"
IN=${OUT}
OUT+="_uniq"
if ! uniqueify -f FreeR_flag ${IN}.mtz ${OUT}.mtz
then
    printf "%bFailed to list unique reflections!%b\n" "${R}" "${W}"
    exit 1
fi
mv "$(basename ${OUT}.mtz .mtz)".log "$(dirname "${OUT}")"/.


# #   --- a. partial occupancy refinement ---
# #  INFO: A 
printf "%s Refinement partial occupancy\n" "$(date +"%y-%m-%d %H:%M")"
IN=${OUT}
# OUT+="_por"
HKLIN=${IN}.mtz
# echo Refining light with "$OCCUPANCY" partial occupancy
phenix.refine \
    pdb.file_name="${COORDINATES}" \
    xray_data.file_name="${HKLIN}" \
    xray_data.labels="F,SIGF" \
    output.prefix="${OUT}" \
    map.file_name="${OUT}".ccp4 \
    "$PARAMETERS" \
    --quiet --overwrite


OUT+="_001"

R_FACTORS=$(grep 'R-work' ${OUT}.log|tail -n 1 |xargs)
if [[ "${R_FACTORS:17:4}" -lt 2500 ]]
then
    printf "%s    %b%s%b\n" "$(date +"%y-%m-%d %H:%M")" "${G}" "${R_FACTORS}" "${W}"
else
    printf "%s    %b%s%b\n" "$(date +"%y-%m-%d %H:%M")" "${R}" "${R_FACTORS}" "${W}"
fi

# # CONVERT HETATM entries to ATOM entries
sed 's/HETATM/ATOM  /' ${OUT}.pdb > ${OUT}_hetatm2atom.pdb
