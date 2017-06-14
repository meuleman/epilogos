This is the repository for the software that computes Epilogos components
from observations of states at sites across the genome.

More info will be added to this README file soon.

In the meantime, here's a very brief tutorial:

* A single script, `computeEpilogos.sh`, does all the processing.
Before you run the script, you need to make the executable programs
that the script calls.  This is done by running `./Makefile` from
this directory.  This will make 3 programs, and place them into a
parallel directory named `bin`.
* After you have made these programs, and before you run the script,
you need to do one more thing:  Enable the script to find the `bin`
directory so it can access and run those programs. To do this,
add the full path to the `bin` directory to your PATH environment
variable.  (If/when all of this becomes available as a module,
this will be taken care of for you when you load the module.)
* Run the script `computeEpilogos.sh` with no arguments to see
the list of arguments it requires and descriptions of them.
* Then you can run the script and supply your arguments and get your
results.  All processing will be done on the cluster, managed by SLURM.
For the time being, the cluster/queue name is hardcoded in the script.


