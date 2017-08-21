# MADE-push-pull-code

The following will contain code for:
1. Creating high resolution K field using Fractional Brownian motion.
- BridgeFast3d_cube_facies_NoFigs.m
- KrigsFast3d.m
- conditioning3d.txt
2. Averaging high resolution K field, and adding vadose zone and bedrock
- PerFK.m
- scaling_discretization.m
- upscaledK.m
3. Parflow input scripts to create pressure files for the duration of the test
- MADE_Facies_K_v1.tcl
- sres_homog.pfb
- vga_homog.pfb
- vgn_homog.pfb
4. Slim-FAST input files for the injection, wait, and extract periods for an ADE run. 
- madesite_slim_inject.tcl
- madesite_slim_wait.tcl
- madesite_slim_extract.tcl
-facies_inject_press2.txt
- facies_wait_press2.txt
- facies_extract_press2.txt
- inject_time.txt
- wait_time.txt
- extract_time.txt
5. PEST files for ADE and t-fADE model optimization. Control file, instruction file, template file
- pest5.pst
- pest.ins
- pest_ADE.tpl
- pest_ADE.par
- parameter.m
6. Smoother function to smooth simulated BTC
- smoother.m
