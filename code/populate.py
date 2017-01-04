"""Code for creating a population of FRBs"""

import sys

from argparse import ArgumentParser

from log import Log

assert sys.version_info >= (3,0), 'Please run with Python3'

if __name__ == '__main__':
    
    parser = ArgumentParser(description='Generate a population of FRBs')

    parser.add_argument('-n',
                        type=int, 
                        required=False,
                        help='number of FRBs to generate/detect')

    parser.add_argument('-nl',
                        '--nolog', 
                        action='store_true',
                        help="don't save log to file")

    parser.add_argument('-ll',
                        '--logloc',
                        default=None,
                        help='location of log file')
                      
    parser.add_argument('-v',
                        '--verbose',
                        action='store_true',
                        help='get more output in the console')
                        
    parser.add_argument('-q',
                        '--quiet',
                        action='store_true',
                        help='turn off output in the console')
    
    args = parser.parse_args()

    logger = Log(save=not args.nolog, 
                 verbose=args.verbose,
                 quiet=args.quiet,
                 loc=args.logloc).logger()
            
    logger.debug('Executing ' + ' '.join(sys.argv))
    logger.info('Hmm')
