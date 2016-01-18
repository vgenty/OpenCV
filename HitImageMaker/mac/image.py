import sys

if len(sys.argv) < 2:
    msg  = '\n'
    msg += "Usage 1: %s $INPUT_ROOT_FILE(s)\n" % sys.argv[0]
    msg += '\n'
    sys.stderr.write(msg)
    sys.exit(1)

from larlite import larlite as fmwk

my_proc = fmwk.ana_processor()

for x in xrange(len(sys.argv)-1):    
    my_proc.add_input_file(sys.argv[x+1])

my_proc.set_io_mode(fmwk.storage_manager.kREAD)

my_unit = fmwk.HitImageManager()

# Set the image maker
image_maker = fmwk.AlgoImageMaker()
my_unit.SetImageMaker(image_maker)

# Add to ana_processor
my_proc.add_process(my_unit)

my_proc.run()

sys.exit(0)
