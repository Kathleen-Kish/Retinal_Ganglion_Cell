# Retinal_Ganglion_Cell

Kathleen Kish
July 29, 2022

Source files to model extracellular stimulation of a retinal ganglion cell. Associated with the following paper: Kathleen Kish, Scott Lempka, and James Weiland. "Modelling extracellular stimulation of retinal ganglion cells: theoretical and practical aspects." 

Included is the Python code (jupyter notebook) to generate a morphologically-realistic mammalian retinal ganglion cell, as well as .mod files for biophysical modeling in the NEURON software environment. 

The file "RGC_Model.ipynb" is the base script that generates the neuron model, and allows for 
  1. Intracellular Current Injection
  2. Extracellular Stimulation from an Ideal Point Source
  3. Extracellular Stimulation from a Disc Electrode (where extracellular voltage is imported from a .txt file)
  
The files "Axon_Ellipse.ipynb" and "Extracellular_Stimulation.ipynb" contain functions called by the base script that generate the axon trajectory and apply an extracellular stimulus pulse, respectively. 

The files "capump.mod", "mammalian_spike.mod", and "xtra.mod" contain the biophysical equations used by the NEURON simulation environment (https://neuron.yale.edu/neuron/). These .mod files must be compiled prior to using the model. 

The files "interpxyz.hoc", "setpointers.hoc", and "RGC_morph_4.hoc" contain functions called by the base script that relate to the 3-D morphology of the cell. 

The file "Voltages_from_COMSOL.txt" contains the extracellular voltage at the center of each neural compartment, and was generated by our finite element model of an epiretinal disc electrode near the retinal surface.
