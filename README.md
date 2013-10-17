Stacking
========
<dl>
<dt>
Repository for the Stacking Method and the Caustic Technique
</dt>
</dl>
______


**caustic_mass_stack2D.py:**


  This is the "main program". It writes out the final data to a specified location in pickle files, which are ###.pkl. 

  

**caustic_class_stack2D.py:**
  
  This script contains classes and functions used specifically by the self-stacking method. 
  

**caustic_universal_stack2D.py:**
  
  This script contains functions used universally by a self-stacking, bin-stacking or any other kind of stacking method.
  

**flux_stack_pbs.sh:**
  
  This is the .sh file that is submitted to flux to run caustic_mass_stack2D.py over a multiple job array.
  

**flux_stack_recovery.py:**
  
  This script contains ss_recover() and bs_recover() functions, which correspond to self stack recover and bin stack recover respectively.
  These functions take the ####.pkl files that were output from caustic_mass_stack2D.py, and 're-loads' them into a working python shell to interactively work with the data.
  
