#!/bin/bash

# Single switch for system selection.
# Allowed examples: O-O, O-O-DCAz05, O-O-TPCCrossedRows90, O-O-TPCClusters70, O-O-Chi^2_TPCClusters4, O-O-Nch, Ne-Ne, Ne-Ne-Nch, p-O, p-O-Nch, pp, pp-Nch
ACTIVE_SYSTEM="${1:-O-O-TPCClusters70}"

echo "[run.sh] Active system: ${ACTIVE_SYSTEM}"

root -l -b -q "Process_dPhidEta.cxx(\"${ACTIVE_SYSTEM}\")"
root -l -b -q "Process_CreateBootstrapSample.cxx(\"${ACTIVE_SYSTEM}\")"
root -l -b -q "Process_TemplateFit.cxx(\"${ACTIVE_SYSTEM}\")"
root -l -b -q "Process_3times2PC.cxx(\"${ACTIVE_SYSTEM}\")"

case "${ACTIVE_SYSTEM}" in
	"O-O-Nch")
		root -l -b -q Plot3times2PC_EtaDiff_OO_Nch10_50_Rings.cxx
		;;
	"Ne-Ne-Nch")
		root -l -b -q Plot3times2PC_EtaDiff_NeNe_Nch10_50_Rings.cxx
		;;
	"p-O-Nch")
		root -l -b -q Plot3times2PC_EtaDiff_pO_Nch10_50_Rings.cxx
		;;
	"pp-Nch")
		root -l -b -q Plot3times2PC_EtaDiff_pp_Nch10_50_Rings.cxx
		;;
esac

if [ "${ACTIVE_SYSTEM}" = "O-O-DCAz05" ]; then
	root -l -b -q Compare3times2PC_EtaDiff_TwoSystems_OO_DCAz05_Rings.cxx
	root -l -b -q Compare_Barlow_OO_DCAz05.cxx
fi

if [ "${ACTIVE_SYSTEM}" = "O-O-TPCCrossedRows90" ]; then
	root -l -b -q Compare3times2PC_EtaDiff_TwoSystems_OO_TPCCrossedRows90_Rings.cxx
	root -l -b -q Compare_Barlow_OO_TPCCrossedRows90_FinalRings.cxx
fi

if [ "${ACTIVE_SYSTEM}" = "O-O-TPCClusters70" ]; then
	root -l -b -q Compare3times2PC_EtaDiff_TwoSystems_OO_TPCClusters70_Rings.cxx
	root -l -b -q Compare_Barlow_OO_TPCClusters70_FinalRings.cxx
fi

if [ "${ACTIVE_SYSTEM}" = "O-O-Chi^2_TPCClusters4" ]; then
	root -l -b -q Compare3times2PC_EtaDiff_TwoSystems_OO_Chi2_TPCClusters4_Rings.cxx
	root -l -b -q Compare_Barlow_OO_Chi2_TPCClusters4_FinalRings.cxx
fi

# Run this after all four systems are generated to build the shared comparison plot.
# root -l -b -q Compare3times2PC_EtaDiff_FourSystems_Nch.cxx
# root -l -b -q Compare3times2PC_EtaDiff_TwoSystems_OO_NeNe_Rings.cxx
# root -l -b -q Compare3times2PC_EtaDiff_TwoSystems_OO_NeNe.cxx
# root -l -b -q Process_FourierFit.cxx
# root -l -b -q Compare3times2PC_EtaDiff_FourSystems.cxx
# root -l -b -q Compare3times2PC_EtaDiff_ThreeSystems.cxx
# root -l -b -q Plot3times2PC_EtaDiff_OO_Nch10_50_Rings.cxx
# root -l -b -q Plot3times2PC_EtaDiff_NeNe_Nch10_50_Rings.cxx
# root -l -b -q Plot3times2PC_EtaDiff_pO_Nch10_50_Rings.cxx
# root -l -b -q Plot3times2PC_EtaDiff_pp_Nch10_50_Rings.cxx