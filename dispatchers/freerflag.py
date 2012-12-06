#!/usr/bin/env python
#
#     Copyright (C) 2012 David Waterman
#
#     This code is distributed under the terms and conditions of the
#     CCP4 Program Suite Licence Agreement as a CCP4 Application.
#     A copy of the CCP4 licence can be obtained by writing to the
#     CCP4 Secretary, Daresbury Laboratory, Warrington WA4 4AD, UK.
#
"""CCP4 dispatcher for the binary 'freerflag'"""

import os
import sys
import shlex
import subprocess
from threading import Thread
from Queue import Queue, Empty
from time import sleep

class Dispatcher:
    """A class to handle dispatch to 'freerflag'"""

    def __init__(self, communicate = True):
        
        self.prog = "freerflag"
        
        # A string giving the command required for dispatch 
        self.command = os.path.join(os.environ["CBIN"], "freerflag")
        
        # Work out the path of this dispatcher and its parent package
        self.dispatcher_path = os.path.realpath(__file__)
        self.dispatcher_package_path = os.path.dirname(self.dispatcher_path)
        
        # This dispatcher points to the specified binary located in:
        self.prog_dir_path = os.environ["CBIN"]
        
        # The default call has empty command line arguments
        self.cmd_args = []
        
        # The default call passes no keywords
        self.keywords = None
        
        # handle to the subprocess
        self.process = None
        
        # is the process running?
        self.isRunning = False

        # Interact with the process by controlling stdin and stdout?
        self.communicate = communicate
        self.stdin = subprocess.PIPE if communicate else None
        self.stdout = subprocess.PIPE if communicate else None
        self.stderr = subprocess.PIPE if communicate else None
        
        # stdout and stderr data from the call (only used if !communicate)
        self.stdout_data = None
        self.stderr_data = None
        
        # The exit code of the executable
        self.call_val = None
        
        # The exception from a failed call
        self.call_err = None
        
        # Suppression of console windows on Windows
        try:
            self.startupinfo = subprocess.STARTUPINFO()
            try:
                self.startupinfo.dwFlags |= subprocess.STARTF_USESHOWWINDOW
            except: #Python > 2.7
                self.startupinfo.dwFlags |= subprocess._subprocess.STARTF_USESHOWWINDOW
        except:
            self.startupinfo = None
        
    def set_env(self):
        """Method to set the environment variables"""
        
        # Set specified environment variables
        os.environ['CCP4I_TCLTK']='/home/rmk65/opt/tcltk/8.4.19/bin'
        os.environ['CCP4_HELPDIR']='/home/rmk65/opt/ccp4/6.2.992/ccp4-6.2.992/help/'
        os.environ['MANPATH']='/home/rmk65/opt/ccp4/6.2.992/ccp4-6.2.992/man:/home/rmk65/opt/ccp4/6.2.992/ccp4-6.2.992/man:/home/rmk65/opt/ccp4/6.2.992/ccp4-6.2.992/man:'
        os.environ['CHTML']='/home/rmk65/opt/ccp4/6.2.992/ccp4-6.2.992/html'
        os.environ['XIA2CORE_ROOT']='/home/rmk65/opt/ccp4/6.2.992/ccp4-6.2.992/share/xia2/xia2core/'
        os.environ['XIA2_ROOT']='/home/rmk65/opt/ccp4/6.2.992/ccp4-6.2.992/share/xia2/xia2/'
        os.environ['PUBLIC_FONT84']='/home/rmk65/opt/ccp4/6.2.992/ccp4-6.2.992/lib/font84.dat'
        os.environ['CCP4']='/home/rmk65/opt/ccp4/6.2.992/ccp4-6.2.992'
        os.environ['CLIB']='/home/rmk65/opt/ccp4/6.2.992/ccp4-6.2.992/lib'
        os.environ['CEXAM']='/home/rmk65/opt/ccp4/6.2.992/ccp4-6.2.992/examples'
        os.environ['IMOSFLM_VERSION']='1.0.5'
        os.environ['CCP4_LIB']='/home/rmk65/opt/ccp4/6.2.992/ccp4-6.2.992/lib'
        os.environ['LD_LIBRARY_PATH']='/home/rmk65/opt/ccp4/6.2.992/ccp4-6.2.992/lib:/home/rmk65/opt/ccp4/6.2.992/ccp4-6.2.992/lib'
        os.environ['BALBES_ROOT']='/home/rmk65/opt/ccp4/6.2.992/ccp4-6.2.992/share/BALBES/Package'
        os.environ['CCP4_BIN']='/home/rmk65/opt/ccp4/6.2.992/ccp4-6.2.992/bin'
        os.environ['HARVESTHOME']='/home/rmk65'
        os.environ['RASMOLPATH']='/home/rmk65/opt/ccp4/6.2.992/ccp4-6.2.992/x-windows/RasMol'
        os.environ['CCP4I_HELP']='/home/rmk65/opt/ccp4/6.2.992/ccp4-6.2.992/share/ccp4i/help'
        os.environ['CCP4I_TOP']='/home/rmk65/opt/ccp4/6.2.992/ccp4-6.2.992/share/ccp4i'
        os.environ['CDOC']='/home/rmk65/opt/ccp4/6.2.992/ccp4-6.2.992/doc'
        os.environ['PISA_CONF_FILE']='/home/rmk65/opt/ccp4/6.2.992/ccp4-6.2.992/share/pisa/pisa.cfg'
        os.environ['PATH']='/home/rmk65/opt/arpwarp/arp_warp_7.2/bin/bin-x86_64-Linux:/home/rmk65/opt/arpwarp/arp_warp_7.2/bin/bin-x86_64-Linux:/home/rmk65/opt/arpwarp/arp_warp_7.2/bin/bin-x86_64-Linux:/home/rmk65/opt/ccp4/6.2.992/ccp4-6.2.992/share/dbccp4i/bin:/home/rmk65/opt/ccp4/6.2.992/ccp4-6.2.992/share/ccp4i/bin:/home/rmk65/opt/ccp4/6.2.992/ccp4-6.2.992/bin:/home/rmk65/opt/ccp4/6.2.992/ccp4-6.2.992/etc:/usr/local/bin:/usr/bin:/bin:/usr/games:/usr/local/sbin:/usr/sbin:/sbin:/home/rmk65/bin:/home/rmk65/opt/qtmg/ccp4mg-2.5.3/bin:/home/rmk65/opt/ccp4/6.2.992/ccp4-6.2.992/share/xia2/xia2core//Test:/home/rmk65/opt/ccp4/6.2.992/ccp4-6.2.992/share/xia2/xia2//Applications:/home/rmk65/.local/bin:/home/rmk65/bin:/home/rmk65/bin:/home/rmk65/opt/qtmg/ccp4mg-2.5.3/bin:/home/rmk65/opt/ccp4/6.2.992/ccp4-6.2.992/share/xia2/xia2core//Test:/home/rmk65/opt/ccp4/6.2.992/ccp4-6.2.992/share/xia2/xia2//Applications:/home/rmk65/opt/ccp4/6.2.992/ccp4-6.2.992/share/xia2/xia2core//Test:/home/rmk65/opt/ccp4/6.2.992/ccp4-6.2.992/share/xia2/xia2//Applications'
        os.environ['XIA2_HOME']='/home/rmk65/opt/ccp4/6.2.992/ccp4-6.2.992/share/xia2'
        os.environ['CCP4_CRANK']='1'
        os.environ['MOSFLM_EXEC']='/home/rmk65/opt/ccp4/6.2.992/ccp4-6.2.992/bin/ipmosflm'
        os.environ['MOSFLM_WISH']='/home/rmk65/opt/tcltk/8.4.19/bin/wish'
        os.environ['DBCCP4I_TOP']='/home/rmk65/opt/ccp4/6.2.992/ccp4-6.2.992/share/dbccp4i'
        os.environ['CCP4_OPEN']='UNKNOWN'
        os.environ['CLIBS']='/home/rmk65/opt/ccp4/6.2.992/ccp4-6.2.992/lib/src'
        os.environ['CETC']='/home/rmk65/opt/ccp4/6.2.992/ccp4-6.2.992/etc'
        os.environ['PLOT_COMMAND']='lp -s'
        os.environ['VERSION']='0.3.4.0'
        os.environ['CLIBD_MON']='/home/rmk65/opt/ccp4/6.2.992/ccp4-6.2.992/lib/data/monomers/'
        os.environ['CCP4_SCR']='/tmp/rmk65'
        os.environ['MOLREPLIB']='/home/rmk65/opt/ccp4/6.2.992/ccp4-6.2.992/lib/data/monomers/'
        os.environ['SHLVL']='1'
        os.environ['GFORTRAN_UNBUFFERED_ALL']='1'
        os.environ['CRANK']='/home/rmk65/opt/ccp4/6.2.992/ccp4-6.2.992/ccp4i/crank'
        os.environ['DYLD_LIBRARY_PATH']='/home/rmk65/opt/ccp4/6.2.992/ccp4-6.2.992/lib:/home/rmk65/opt/ccp4/6.2.992/ccp4-6.2.992/lib:/home/rmk65/opt/ccp4/6.2.992/ccp4-6.2.992/lib'
        os.environ['PYTHONPATH']='/home/rmk65/opt/ccp4/6.2.992/ccp4-6.2.992/share/python:/home/rmk65/opt/ccp4/6.2.992/ccp4-6.2.992/share/python:/home/rmk65/opt/ccp4/6.2.992/ccp4-6.2.992/share/python:'
        os.environ['CLIBD']='/home/rmk65/opt/ccp4/6.2.992/ccp4-6.2.992/lib/data'
        os.environ['CLASSPATH']='/home/rmk65/opt/ccp4/6.2.992/ccp4-6.2.992/bin:/home/rmk65/opt/ccp4/6.2.992/ccp4-6.2.992/bin:/home/rmk65/opt/ccp4/6.2.992/ccp4-6.2.992/bin:'
        os.environ['CBIN']='/home/rmk65/opt/ccp4/6.2.992/ccp4-6.2.992/bin'
        os.environ['CINCL']='/home/rmk65/opt/ccp4/6.2.992/ccp4-6.2.992/include'
        os.environ['CCP4_BROWSER']='firefox'
        os.environ['CCP4_MASTER']='/home/rmk65/opt/ccp4/6.2.992'
        os.environ['MMCIFDIC']='/home/rmk65/opt/ccp4/6.2.992/ccp4-6.2.992/share/ccif/cif_mmdic.lib'
        os.environ['CPROG']='/home/rmk65/opt/ccp4/6.2.992/ccp4-6.2.992/src'
        os.environ['PRINT_COMMAND']='lp -s'
        os.environ['BINSORT_SCR']='/tmp/rmk65'
        os.environ['_']='/usr/bin/printenv'
        
        # Prepend dispatcher package to the PATH, if not included already
        try:
            curr_path = os.environ["PATH"]
        except KeyError:
            curr_path = ""
        if not self.dispatcher_package_path in curr_path:
            os.environ["PATH"] = (self.dispatcher_package_path +
                                  os.pathsep + curr_path)
            
        # Prepend dispatcher package to PYTHONPATH, if not included already
        try:
            curr_pythonpath = os.environ["PYTHONPATH"]
        except KeyError:
            curr_pythonpath = ""
        if not self.dispatcher_package_path in curr_pythonpath:
            os.environ["PYTHONPATH"] = (self.dispatcher_package_path +
                                        os.pathsep + curr_pythonpath)

    def get_bin(self):
        """Method to get the path to 'freerflag'"""
        
        return os.path.join(self.prog_dir_path, self.prog)
    
    def set_cmd_args(self, args):
        """Method to append arguments to the program call"""
        
        # If the arguments are supplied as a string, split it suitably.
        # If it is a list, pass it on untouched
        if isinstance(args, str):
            # Preserve backspace separated filenames for Windows
            self.cmd_args = shlex.split(args, posix=(os.name != "nt"))
        elif isinstance(args, list):
            self.cmd_args = args
    
    def set_keywords(self, kwstring):
        """Method to set keywords to pass to 'freerflag'"""

        if isinstance(kwstring, list):
            kwstring = "\n".join(kwstring)

        if not kwstring.endswith("\n"):
            kwstring = kwstring + "\n"

        self.keywords = kwstring
        
    def call(self, wait=True):
        """
        Method to execute freerflag.
        
        When communicating with the subprocess, wait=True avoids the need
        to repeatedly call monitor, but blocks until the subprocess completes.
        """
        
        cmd = shlex.split(self.command, posix=(os.name != "nt"))
        popenargs = cmd + self.cmd_args
        
        try:
            self.process = subprocess.Popen(popenargs, stdin=self.stdin,
                                 stdout=self.stdout,
                                 stderr=self.stderr,
                                 startupinfo=self.startupinfo)
        except OSError, e: 
            self.call_err = sys.exc_info()
            return self.call_err
        
        self.call_val = None # Reset this in case this was not the first call
        self.isRunning = True
        
        if self.communicate:

            if wait:
                (self.stdout_data,
                 self.stderr_data) = self.process.communicate(self.keywords)
                self.isRunning = False
                self.stdout_data = self.stdout_data.splitlines()
                self.stderr_data = self.stderr_data.splitlines()
                 
            else:
                if self.keywords:
                    self.process.stdin.write(self.keywords)
                self.process.stdin.close()
                
                # create lists to be populated by the monitor method
                self.stdout_data = []
                self.stderr_data = []
                
                # queues to buffer stdout and stderr data
                self.stdout_queue = Queue()
                self.stderr_queue = Queue()
                
                # start threads to read output into the queues
                self._stop_reader = False
                self._stdout_reader = Thread(target = self.__enqueue,
                                       args = (self.process.stdout, self.stdout_queue))
                self._stderr_reader = Thread(target = self.__enqueue,
                                       args = (self.process.stderr, self.stdout_queue))
                self._stdout_reader.daemon = True
                self._stderr_reader.daemon = True
                self._stdout_reader.start()
                self._stderr_reader.start()
                
                # set returncode if already finished
                self.process.poll()
                    
        else:
            self.process.wait()
            self.isRunning = False
            
        self.call_val = self.process.returncode
        return self.call_val
    
    def __enqueue(self, stream, queue):
        """Worker thread function for call with communicate=True and wait=False"""
        
        while True:
            try:
                s = stream.readline()
                if s:
                    queue.put(s)
                else: # s is the empty string or None
                    if self._stop_reader: return
                    sleep(0.1) # avoid too much spinning in this loop
            except: # the stream was closed, so finish
                break
        return None
    
    def __stop_readers(self):
        """Method to cleanly stop the reader threads"""
        
        self.process.stdout.flush()
        self.process.stderr.flush()
        self._stop_reader = True
        if self._stdout_reader.is_alive(): self._stdout_reader.join()
        if self._stderr_reader.is_alive(): self._stderr_reader.join()
        self.process.stdout.close()
        self.process.stderr.close()

        return None
    
    def monitor(self):
        """
        Method to read output of freerflag and check if it has finished.
        
        The user of this method must call it repeatedly to fill stdout_data
        and stderr_data, until call_val is set, or abort is called.
        """
        
        # Use this method only if communicating with a running subprocess
        if not self.communicate: return None
        if not self.isRunning: return None
        
        self.call_val = self.process.poll()

        # read from the queues, allowing time to accumulate if nearly empty
        try: o = self.stdout_queue.get_nowait()
        except Empty: o = None
        try: e = self.stderr_queue.get_nowait()
        except Empty: e = None
        
        if o: self.stdout_data.append(o.rstrip())
        if e: self.stderr_data.append(e.rstrip())
        
        # process finished, clean up
        if self.call_val is not None:
            # first time in this block, end reader threads
            if not self._stop_reader:
                self.__stop_readers()
            # subsequently, check if nothing left in the queues
            elif not o and not e: self.isRunning = False
            else: pass
            
        return o, e
    
    def abort(self):
        """Method to terminate or kill the process"""

        # Use this method only if communicating with a running subprocess
        if not self.communicate: return None
        if not self.isRunning: return None
        
        self.process.terminate()
        if not self.process.poll() is None:
            self.process.kill()
        
        # ask reader threads to close
        self.__stop_readers()
        
        # close the streams
        self.process.stdout.close()
        self.process.stderr.close()
        self.call_val = self.process.poll()
        self.isRunning = False
        
if __name__ == "__main__":
    # This dispatcher is being run as a script, so go ahead and execute freerflag
    
    # Instantiate the dispatcher
    dispatcher = Dispatcher(communicate = False)
    
    # Set the environment
    # dispatcher.set_env()
    
    # Run the program now and exit with its value, or re-raise any error
    # encountered, to exit the interpreter
    dispatcher.set_cmd_args(sys.argv[1:])
    dispatcher.call()
    if dispatcher.call_err:
        raise dispatcher.call_err[1], None, dispatcher.call_err[2]
    else:
        sys.exit(dispatcher.call_val)
