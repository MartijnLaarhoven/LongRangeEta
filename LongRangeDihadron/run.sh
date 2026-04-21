#!/bin/bash

# Single switch for system selection.
# Allowed examples: O-O, O-O-DCAz05, O-O-Nch, Ne-Ne, Ne-Ne-Nch, p-O, p-O-Nch, pp, pp-Nch
ACTIVE_SYSTEM="${1:-O-O}"

echo "[run.sh] Active system: ${ACTIVE_SYSTEM}"

root -l -b -q "Process_dPhidEta.cxx(\"${ACTIVE_SYSTEM}\")"
root -l -b -q "Process_CreateBootstrapSample.cxx(\"${ACTIVE_SYSTEM}\")"
root -l -b -q "Process_TemplateFit.cxx(\"${ACTIVE_SYSTEM}\")"
root -l -b -q "Process_3times2PC.cxx(\"${ACTIVE_SYSTEM}\")"

case "${ACTIVE_SYSTEM}" in
	"O-O-DCAz05")
		root -l -b -q Compare_Barlow_OO_DCAz05.cxx
		;;
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

# Run this after all four systems are generated to build the shared comparison plot.
# root -l -b -q Compare3times2PC_EtaDiff_FourSystems_Nch.cxx
root -l -b -q Compare3times2PC_EtaDiff_TwoSystems_OO_NeNe_Rings.cxx
# root -l -b -q Compare3times2PC_EtaDiff_TwoSystems_OO_NeNe.cxx
# root -l -b -q Process_FourierFit.cxx
# root -l -b -q Compare3times2PC_EtaDiff_FourSystems.cxx
# root -l -b -q Compare3times2PC_EtaDiff_ThreeSystems.cxx
# root -l -b -q Plot3times2PC_EtaDiff_OO_Nch10_50_Rings.cxx
# root -l -b -q Plot3times2PC_EtaDiff_NeNe_Nch10_50_Rings.cxx
# root -l -b -q Plot3times2PC_EtaDiff_pO_Nch10_50_Rings.cxx
# root -l -b -q Plot3times2PC_EtaDiff_pp_Nch10_50_Rings.cxx