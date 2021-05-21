import sys
import click
from epilogos.epilogos import main as runEpilogos


print("""\n
                  d8b 888                                     
                  Y8P 888                                     
                      888                                     
 .d88b.  88888b.  888 888  .d88b.   .d88b.   .d88b.  .d8888b  
d8P  Y8b 888 "88b 888 888 d88""88b d88P"88b d88""88b 88K      
88888888 888  888 888 888 888  888 888  888 888  888 "Y8888b. 
Y8b.     888 d88P 888 888 Y88..88P Y88b 888 Y88..88P      X88 
 "Y8888  88888P"  888 888  "Y88P"   "Y88888  "Y88P"   88888P' 
         888                            888                   
         888                       Y8b d88P                   
         888                        "Y88P"                    
""", flush=True)

if len(sys.argv) == 1:
    print("Run 'epilogos -h' for help")
    sys.exit()

@click.command(context_settings=dict(help_option_names=['-h', '--help']))
@click.option("-m", "--mode", "mode", type=click.Choice(["single", "paired"]), default=["single"], show_default=True,
              multiple=True, help="single for single group epilogos and paired for 2 group epilogos")
@click.option("-l", "--local", "commandLineBool", is_flag=True, multiple=True,
              help="If enabled, Epilogos will run locally in your terminal rather than on a SLURM cluster")
@click.option("-i", "--input-directory", "inputDirectory", type=str, multiple=True,
              help="Path to directory that contains files to read from (ALL files in this directory will be read in)")
@click.option("-a", "--directory-one", "inputDirectory1", type=str, multiple=True,
              help="Path to first directory that contains files to read from (ALL files in this directory will be read in)")
@click.option("-b", "--directory-two", "inputDirectory2", type=str, multiple=True,
              help="Path to second directory that contains files to read from (ALL files in this directory will be read in)")
@click.option("-o", "--output-directory", "outputDirectory", type=str, multiple=True,
              help="Output Directory (CANNOT be the same as input directory)\n")
@click.option("-n", "--state-info", "stateInfo", type=str, multiple=True, help="State model info file")
@click.option("-s", "--saliency", "saliency", type=int, default=[1], show_default=True, multiple=True,
              help="Desired saliency level (1, 2, or 3)")
@click.option("-c", "--num-cores", "numProcesses", type=int, default=[0], multiple=True,
              help="The number of cores to run on [default: 0 = Uses all cores]")
@click.option("-x", "--exit", "exitBool", is_flag=True, multiple=True,
              help="If flag is enabled program will exit upon submission of all SLURM processes rather than completion of " +
                   "all processes")
@click.option("-d", "--diagnostic-figures", "diagnosticBool", is_flag=True, multiple=True,
              help="If flag is enabled, program will output some diagnostic figures of the gennorm fit on the null data and " +
                   "comparisons between the null and real data")
@click.option("-t", "--num-trials", "numTrials", type=int, default=[101], show_default=True, multiple=True,
              help="The number of times subsamples of the scores are fit")
@click.option("-z", "--sampling-size", "samplingSize", type=int, default=[100000], show_default=True, multiple=True,
              help="The size of the subsamples on which the scores are fit")
@click.option("-q", "--quiescent-state", "quiescentState", type=int, multiple=True, 
              help="If a bin contains only states of this value, it is treated as quiescent and not factored into fitting." + 
                   "If set to 0, filtering is not done. [default: last state]")
@click.option("-g", "--group-size", "groupSize", type=int, default=[-1], show_default=True, multiple=True,
              help="In pairwise epilogos controls the sizes of the shuffled arrays. Default is sizes of the input groups")
def main(mode, commandLineBool, inputDirectory, inputDirectory1, inputDirectory2, outputDirectory, stateInfo, saliency,
    numProcesses, exitBool, diagnosticBool, numTrials, samplingSize, quiescentState, groupSize):
    """
    Information-theoretic navigation of multi-tissue functional genomic annotations

    Written by Jacob Quon and Wouter Meuleman
    """
    runEpilogos(mode, commandLineBool, inputDirectory, inputDirectory1, inputDirectory2, outputDirectory, stateInfo, saliency,
    numProcesses, exitBool, diagnosticBool, numTrials, samplingSize, quiescentState, groupSize)


if __name__ == "__main__":
    prog_name="epilogos"
    main()