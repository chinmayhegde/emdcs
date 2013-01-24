In order to use the Matlab interface to the EMD flow algorith, do the
following:

1) Build the lemon library (used for min-cost max flow)
1.1) cd to emd_flow/lemon
1.2) execute the prepare script in emd_flow/lemon: ./prepare.sh

2) Build the emd_flow code and the matlab interface
2.1) cd back to emd_flow
2.2) make sure mex (the matlab compiler) is on your $PATH
2.3) build the mex file: make mexfile
