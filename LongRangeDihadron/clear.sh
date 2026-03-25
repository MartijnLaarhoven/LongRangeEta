#!/bin/bash

# Clear output directories and log files
echo "Clearing output directories and log files..."

rm -rf ProcessOutput/
rm -rf 3times2PC/
rm -rf FourierFit/
rm -rf TemplateFit/

rm -f core_dump_*

echo "Done! All output files cleared."
