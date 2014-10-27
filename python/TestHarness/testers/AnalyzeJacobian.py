import re, os, sys
from Tester import Tester
from RunParallel import RunParallel # For TIMEOUT value

class AnalyzeJacobian(Tester):

  @staticmethod
  def validParams():
    params = Tester.validParams()
    params.addRequiredParam('input',  "The input file to use for this test.")
    params.addParam('test_name',      "The name of the test - populated automatically")
    params.addParam('expect_out',     "A regular expression that must occur in the input in order for the test to be considered passing.")
    params.addParam('resize_mesh', False, "Resize the input mesh")
    params.addParam('mesh_size',   1, "Resize the input mesh")

    return params

  def __init__(self, name, params):
    Tester.__init__(self, name, params)


  def getCommand(self, options):
    specs = self.specs
    # Create the command line string to run
    command = specs['moose_dir'] + '/python/jacobiandebug/analyzejacobian.py'

    # Check for built application
    if not options.dry_run and not os.path.exists(command):
      print 'Application not found: ' + str(specs['executable'])
      sys.exit(1)

    if specs['resize_mesh'] :
      mesh_options = '-r -s %d' % specs['mesh_size']
    else :
      mesh_options = ''

    command += mesh_options + ' ' + specs['input']

    return command


  def processResults(self, moose_dir, retcode, options, output):
    reason = ''
    specs = self.specs
    if specs.isValid('expect_out'):
      out_ok = self.checkOutputForPattern(output, specs['expect_out'])
      if (out_ok and retcode != 0):
        reason = 'OUT FOUND BUT CRASH'
      elif (not out_ok):
        reason = 'NO EXPECTED OUT'
    if reason == '':
      if retcode == RunParallel.TIMEOUT:
        reason = 'TIMEOUT'
      elif retcode != 0 :
        reason = 'CRASH'

    return (reason, output)

  def checkOutputForPattern(self, output, re_pattern):
    if re.search(re_pattern, output, re.MULTILINE | re.DOTALL) == None:
      return False
    else:
      return True
