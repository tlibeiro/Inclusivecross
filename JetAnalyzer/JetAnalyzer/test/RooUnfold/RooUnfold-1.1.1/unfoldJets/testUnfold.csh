#!/bin/csh
#root -l 
gSystem->Load("../libRooUnfold.so");
.include ../src
.L UnfoldJets4.cxx+
UnfoldJets3()
Getline("Type <return> to exit: ");
##EOF

