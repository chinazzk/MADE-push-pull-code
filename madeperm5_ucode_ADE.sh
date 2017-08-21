#!/bin/bash
#My first script

matlab -nodisplay -r "parameter"

cd ~/Madesite_Ben/3_upscaling/madeperm_5slim/SLIM_Facies/ADE

tclsh madesite_slim_inject.tcl

cd ~/Madesite_Ben/3_upscaling/madeperm_5slim/SLIM_Facies/ADE

tclsh madesite_slim_wait.tcl


cd ~/Madesite_Ben/3_upscaling/madeperm_5slim/SLIM_Facies/ADE

tclsh madesite_slim_extract.tcl

cd madesite_slim_5

matlab -nodisplay -r "smoother"

exit