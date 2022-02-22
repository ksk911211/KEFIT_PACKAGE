import time, sys

# update_progress() : Displays or updates a console progress bar
## Accepts a float between 0 and 1. Any int will be converted to a float.
## A value under 0 represents a 'halt'.
## A value at 1 or bigger represents 100%
def update_progress(progress):
    barLength = 50 # Modify this to change the length of the progress bar
    status = ""
    if isinstance(progress, int):
        progress = float(progress)
    if not isinstance(progress, float):
        progress = 0
        status = "error: progress var must be float\r\n"
    if progress < 0:
        progress = 0
        status = "Halt...\r\n"
    if progress >= 1:
        progress = 1
        status = "Done...\r\n"
    block = int(round(barLength*progress))
    text = "\rProgress: [{0}] {1}% {2}".format( "#"*block + "_"*(barLength-block), round(100*progress,1), status)
    sys.stdout.write(text)
    sys.stdout.flush()


# update_progress test script
#print("progress : 'hello'")
#update_progress("hello")
#time.sleep(1)

#for i in range(10):
#	update_progress(0.033*i)
#	time.sleep(1)
#print()
