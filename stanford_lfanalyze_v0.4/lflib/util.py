# This file contains functionl of general utility used throughout LFAnalyze.

# Run an individual command.  Command output is returned as a
# string, along with a bool indicated whether the command ran
# successfully.
def run_command_wrapper(cmd):
    try:
        import subprocess
        result = subprocess.check_output(cmd, stderr=subprocess.STDOUT)  # Grab stderr as well...
        return (True, '\n'+result, ''.join(["%s " % el for el in cmd]))
    except subprocess.CalledProcessError,e:                              # Handle error conditions
        return (False, str(e) + '\nProgram output was:\n' + e.output, ''.join(["%s " % el for el in cmd]))

def ensure_path(path):
    """ Check if a path exists. If not, create the necessary directories, 
    but if the path includes a file, don't create the file"""
    import os, errno
    dir_path = os.path.dirname(path)
    if len(dir_path) > 0:  # Only run makedirs if there is a directory to create!
        try:
            os.makedirs(dir_path)
        except OSError as exception:
            if exception.errno != errno.EEXIST:
                raise
    return path

if __name__ is "__main__":
    pass
