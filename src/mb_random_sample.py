##############################################################################
#
#   A small script that samples one number from a standard
#   Gaussian distribution
#
#   AUTHOR: Maciej_Bak
#   AFFILIATION: Swiss_Institute_of_Bioinformatics
#   CONTACT: maciej.bak@unibas.ch
#   CREATED: 22-11-2019
#   LICENSE: GPL v3.0
#   USAGE: python mb_random_sample.py --outfile [outfile]
#
##############################################################################

# imports
import sys
import os
import time
import logging
from argparse import ArgumentParser, RawTextHelpFormatter
import numpy as np

def parse_arguments():
    '''Parser of the command-line arguments.'''
    parser = ArgumentParser(description=__doc__,
                            formatter_class=RawTextHelpFormatter)
    parser.add_argument("-v",
                        "--verbosity",
                        dest="verbosity",
                        choices=\
                        ('DEBUG', 'INFO', 'WARN', 'ERROR', 'CRITICAL'),
                        default='ERROR',
                        help="Verbosity/Log level. Defaults to ERROR")
    parser.add_argument("-l",
                        "--logfile",
                        dest="logfile",
                        help="Store log to this file.")
    parser.add_argument("--outfile",
                        dest="outfile",
                        required=True,
                        help="Path for the output file.")
    return parser

##############################################################################

def main():
    '''Main body of the script.'''
    with open(options.outfile,"w") as f:
        sample = np.random.normal() # by default: mu=0, sigma=1
        f.write(str(sample))

##############################################################################

if __name__ == '__main__':

    try:
        # parse the command-line arguments
        options = parse_arguments().parse_args()

        # set up logging during the execution
        formatter = logging.Formatter(fmt="[%(asctime)s] %(levelname)s\
                                      - %(message)s",
                                      datefmt="%d-%b-%Y %H:%M:%S")
        console_handler = logging.StreamHandler()
        console_handler.setFormatter(formatter)
        logger = logging.getLogger('uniprot_to_json')
        logger.setLevel(logging.getLevelName(options.verbosity))
        logger.addHandler(console_handler)
        if options.logfile is not None:
            logfile_handler = logging.handlers.RotatingFileHandler(\
                options.logfile, maxBytes=50000, backupCount=2)
            logfile_handler.setFormatter(formatter)
            logger.addHandler(logfile_handler)

        # execute the body of the script
        start_time = time.time()
        logger.info("Starting script")
        main()
        seconds = time.time() - start_time

        # log the execution time
        minutes, seconds = divmod(seconds, 60)
        hours, minutes = divmod(minutes, 60)
        logger.info("Successfully finished in {hours} hour(s) \
{minutes} minute(s) and {seconds} second(s)".format(
        hours=int(hours),
        minutes=int(minutes),
        seconds=int(seconds) if seconds > 1.0 else 1
        ))
    # log the execution time in case of interruption
    except KeyboardInterrupt as e:
        logger.warn("Interrupted by user after {hours} hour(s) \
{minutes} minute(s) and {seconds} second(s)".format(
        hours=int(hours),
        minutes=int(minutes),
        seconds=int(seconds) if seconds > 1.0 else 1
        ))
        sys.exit(-1)
    # log the exception in case it happens
    except Exception as e:
        logger.exception(str(e))
        raise e
