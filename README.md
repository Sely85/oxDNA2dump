oxDNA2dump
==========

Welcome to oxDNA2dump code!
This code converts oxDNA (dna.physics.ox.ac.uk) configuration into a dump file that can be read by OVITO visualization software (ovito.org). 


BUILD (Linux)
-------------
g++ -lm oxDNA2dump.cpp -o oxDNA2dump


USAGE
-----
./oxDNA2dump oxDNAtopology oxDNAconfiguration


EXAMPLE
-------
As an example a double stranded sequence called S1S2 (GCGTCATACAGTGC)
has been generated using the generate-sa.py utils of oxDNA software.
The generated topology (S1S2.top) and configuration (S1S2.dat) can be
converted into a dump file as follows:

./oxDNA2dump S1S2.top S1S2.dat 

The output files are a dump file (oxDNAconv.dump) which can be
uploaded in OVITO Visualization Tool. If the oxDNA configuration file
is a trajectory a dump containing multiple frame will be created. The
column mapping is: id, type, position (x, y, z), quaternion (qw, qx,
qy, qz), aspherical shape (x, y, z), linear velocities (x, y, z),
angular velocity (x, y, z), moleculeID.
An energy file (oxDNAconv.ene) is also written, to collect energy data
contained in the oxDNA configurations file. The order of column is
mantained (timestep, total energy, potential energy and kinetic energy). 
