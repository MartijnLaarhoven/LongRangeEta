#!/bin/bash

root -l -b -q Process_dPhidEta.cxx
root -l -b -q Process_CreateBootstrapSample.cxx
root -l -b -q Process_TemplateFit.cxx
root -l -b -q Process_FourierFit.cxx
root -l -b -q Process_3times2PC.cxx
root -l -b -q Compare3times2PC_EtaDiff_FourSystems.cxx