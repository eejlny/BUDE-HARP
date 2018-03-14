#!/bin/bash
HOSTNAME=`hostname`
echo "Executing OpenCL Code on node $HOSTNAME"
source bdw_fpga_pilot_example_design/bude/init.sh
cd bdw_fpga_pilot_example_design/bude
bin/bude
