## pbs_submit_dijet_imbalance.py

## convenience script to batch submit and monitor
## job status, automatically resubmitting jobs that
## do not produce the expected output files

import os
import argparse
import subprocess
import time
import random
import ROOT

def checkstatus(jobstatus) :
  loop = False
  for i in range(len(jobstatus)) :
    if jobstatus[i] != 2 :
      loop = True
      break
  return loop

def updatestatus(jobstatus, outdir, name) :
  print "Updating job status"
  print "Total: " + str(len(jobstatus))
  
  ## get the qstat job listing
  proccommand = 'qstat | grep dx5412'
  proc = subprocess.Popen( proccommand, stdout=subprocess.PIPE, shell=True)
  qstat_result = proc.stdout.read()
  
  for i in range(len(jobstatus)) :
    
    ## if job is completed, we don't need to check again
    if jobstatus[i] == 2 :
      continue
    
    ## check if the job is still underway
    jobinprocess = qstat_result.find(name + str(i))
    if jobinprocess >= 0 :
      jobstatus[i] = 1
      continue
    
    ## if the job is not still underway,
    ## check to see if the job has completed properly
    ## if not, mark to resubmit
    if outdir.startswith('/') :
      filename = outdir + '/' + name + str(i) + '.root'
    else :
      filename = os.getcwd() + '/' + outdir + '/' + name + str(i) + '.root'
    
    if os.path.isfile(filename) :
      outputfile = ROOT.TFile(filename, "READ")
      if outputfile.IsZombie() :
        print "job " + str(i+1) + " of " + str(len(jobstatus)) + " complete: file is zombie, resubmit"
        jobstatus[i] = 0
        os.remove(filename)
      elif outputfile.IsOpen() :
        print "job " + str(i+1) + " of " + str(len(jobstatus)) + " complete: ROOT file healthy"
        print filename
        jobstatus[i] = 2
        outputfile.Close()
      else :
        print "job " + str(i+1) + " of " + str(len(jobstatus)) + " undefined file status, resubmit"
        jobstatus[i] = 0
    else :
      print "undefined status: job " + str(i+1) + " of " + str(len(jobstatus)) + " marked for submission"
      jobstatus[i] = 0

return jobstatus


def main(args) :
  ## if there are no input files, exit
  files = args.strings
  if files is None :
    return
  
  ## get max number of jobs to be submitted to pbs
  ## at a time
  maxjobs = args.maxjobs
  
  ## read in embedding data file
  ## create a list of embedding data if requested
  embedding_list = []
  if args.embed is not None :
    with open(args.embed, 'r') as fp :
      line = fp.readline()
      while line :
        embedding_list.append(line)
        line = fp.readline()

  ## if no valid embedding data was given, fill w/ empty
  if not embedding_list:
    embedding_list = ['']

  ## some paths
  execpath = os.getcwd()
  executable = './bin/dijet_imbalance/dijet_imbalance'
  ## find the qwrap file
  qwrap = execpath + '/submit/qwrap.sh'

  ## we need to do our own book keeping
  ## for when  job is active and when it
  ## has completed successfully
  
  ## 0 not submitted, 1 for running, 2 for complete,
  ## and -1 for failed - resubmit
  jobstatus = [0 for i in range(len(files))]
  
  ## count the number of qsub submission failures
  qsubfail = 0

  ## select defaults
  embedTrig = ''
  if args.embedTriggers is not None :
    embedTrig = args.embedTriggers

  reader = ''
  if args.readerSetting is not None :
    reader = args.readerSetting

  embedReader = ''
  if args.embedReaderSetting is not None :
    embedReader = args.embedReaderSetting

  triggerEfficiency = ''
  if args.triggerEfficiency is not None :
    triggerEfficiency = args.triggerEfficiency

  embedEfficiency = ''
  if args.embedEfficiency is not None :
    embedEfficiency = args.embedEfficiency

  ## define extra boolean options string
  extraOpts = ''
  if args.triggerEfficiency :
    extraOpts = extraOpts + ' --efficiency'
  if args.embedEfficiency :
    extraOpts = extraOpts + ' --embedEfficiency'
  if args.pp :
    extraOpts = extraOpts + ' --pp'
  if args.forcePPEffCent is not None :
    extraOpts = extraOpts + ' --forcePPCent=' + str(args.forcePPEffCent)
  
  while checkstatus(jobstatus) :
    
    ## update status of all jobs & output
    jobstatus = updatestatus(jobstatus, args.output, args.name)
    
    ## if we have completed all jobs, exit
    if len(jobstatus) == jobstatus.count(2) :
      break
    
    ## find the number of jobs still running via qstat
    ## if its at the maximum set jobsactive or jobsqueue,
    ## then pause
    jobsactive = jobstatus.count(1)
    while jobsactive >= maxjobs :
      print "reached max number of active jobs: pausing"
      time.sleep(60)
      jobstatus = updatestatus(jobstatus, args.output, args.name)
      jobsactive = jobstatus.count(1)
    
    ## now submit jobs up to maxjobs - jobsqueued
    njobs = maxjobs - jobsactive
    
    for i in range(len(jobstatus)) :
      if njobs == 0 :
        break
      if jobstatus[i] == 1 or jobstatus[i] == 2 :
        continue
      ## build log & error locations
      outstream = "log/" + args.name + str(i) + ".log"
      errstream = "log/" + args.name + str(i) + ".err"

      ## select embedding data randomly
      embed_choice = random.choice(embedding_list)

      ## build our qsub execution string
      clargs = '--outDir=' + args.output + ' --input=' + files[i] + ' --id=' + str(i)
      clargs = clargs + ' --name=' + args.name + extraOpts
      clargs = clargs + ' --runList=' + args.badRuns + ' --towList=' + args.badTowers + ' --triggers='
      clargs = clargs + args.triggers + ' --embedTriggers=' + embedTrig + ' --constEta=' + args.constEta
      clargs = clargs + ' --leadConstPt=' + args.leadConstPt + ' --subConstPt=' + args.subConstPt
      clargs = clargs + ' --leadConstPtMatch=' + args.leadConstPtMatch + ' --subConstPtMatch='
      clargs = clargs + args.subConstPtMatch + ' --leadR=' + args.leadR + ' --subR=' + args.subR
      clargs = clargs + ' --leadJetPt=' + args.leadJetPt + ' --subJetPt=' + args.subJetPt
      clargs = clargs + ' --readerSetting=' + reader + ' --embedReaderSetting='
      clargs = clargs + embedReader + ' --reuseTrigger=' + str(args.reuseTrigger)
      clargs = clargs  + ' --embed=' + embed_choice

      
      qsub = 'qsub -V -p ' + str(args.priority) + ' -l mem=' + str(args.mem) + 'GB -l nodes=' + str(args.nodes)
      qsub = qsub + ':ppn=' + str(args.ppn) + ' -q ' + str(args.queue) + ' -o ' + outstream
      qsub = qsub + ' -e ' + errstream + ' -N ' + args.name + str(i) + ' -- '
      qsub = qsub + qwrap + ' ' + execpath + ' ' + executable + ' ' + clargs
      print "submitting job: "
      print qsub
      ret = subprocess.Popen( qsub, shell=True)
      ret.wait()
      if ret.returncode == 0 :
        jobstatus[i] = 1
      else :
        print "warning: qsub submission failed"
        qsubfail = qsubfail + 1
      njobs = njobs - 1
    
    if qsubfail > 100 :
      print "qsub failure too many times - exiting"
      return
    
    ## wait two minutes before rechecking
    print "finished round of submissions: pausing"
    time.sleep(60)
  
  print "all jobs completed: exiting"

if __name__ == "__main__":
  parser = argparse.ArgumentParser(description='Submit jobs via PBS & resubmit if necessary')
  parser.add_argument('strings', metavar='S', nargs='*', help=' input files ')
  parser.add_argument('--name', default='job_', help=' job name (identifier for qsub)')
  parser.add_argument('--mem', type=int, default=2, help=' memory required per job [GB]')
  parser.add_argument('--nodes', type=int, default=1, help=' number of nodes required per job')
  parser.add_argument('--ppn', type=int, default=1, help=' number of processors per node required per job')
  parser.add_argument('--priority', type=int, default=0, help=' queue priority')
  parser.add_argument('--queue', default='erhiq', help=' queue to submit jobs to' )
  parser.add_argument('--maxjobs',type=int, default=100, help=' max number of jobs to have in running or queue states')
  parser.add_argument('--output', default='out/post/tmp', help=' directory for output root files' )
  parser.add_argument('--triggerEfficiency', dest='triggerEfficiency', action='store_true', help='turn on relative efficiency corrections for trigger data (wrt central AuAu)')
  parser.add_argument('--no-triggerEfficiency', dest='triggerEfficiency', action='store_false', help='turn off relative efficiency corrections for trigger data')
  parser.set_defaults(triggerEfficiency=False)
  parser.add_argument('--embedEfficiency', dest='embedEfficiency', action='store_true', help='turn on relative efficiency corrections for the embeding data (wrt central AuAu)')
  parser.add_argument('--no-embedEfficiency', dest='embedEfficiency', action='store_false', help='turn off relative efficiency corrections for embedding data')
  parser.set_defaults(embedEfficiency=False)
  parser.add_argument('--pp', dest='pp', action='store_true', help='specify that the trigger data is p+p data')
  parser.add_argument('--no_pp', dest='pp', action='store_false', help='specify that the trigger data is Au+Au data')
  parser.add_argument('--auau', dest='pp', action='store_false', help='specify that the trigger data is Au+Au data')
  parser.set_defaults(pp=False)
  parser.add_argument('--forcePPEffCent', default=None, help='if efficiency corrections are on, forces the centrality for pp corrections. If not on, the centrality is calculated from embedding, or chosen randomly')
  parser.add_argument('--badRuns', default='submit/y14_bad_run.txt', help=' csv file containing runs to mask')
  parser.add_argument('--badTowers', default='submit/y14_y6_bad_tower.txt', help=' csv file containing towers to mask')
  parser.add_argument('--triggers', default='ALL', help=' event triggers to consider: [y7, y10, y11, y14, y6pp, y9pp, y12pp] + [HT, MB, HT2, HT3, VPDMB30, VPDMB5, MBMON, ALL] (default "ALL": accept all events)')
  parser.add_argument('--embedTriggers', default=None, help=' event triggers to consider for embedding data (see above for options)')
  parser.add_argument('--embed', default=None, help='list file containing root files for pp embedding, not used for AuAu')
  parser.add_argument('--reuseTrigger', type=int, default=1, help='[for pp] how many times a pp event should be embedded')
  parser.add_argument('--readerSetting', default=None, help='can specify non-default reader settings in a text tile')
  parser.add_argument('--embedReaderSetting', default=None, help='can specify non-default reader settings for embedding data in text file')
  parser.add_argument('--constEta', default='1.0', help='list of constituent eta ranges to use during jetfinding')
  parser.add_argument('--leadConstPt', default='2.0', help='list of leading hard jet constituent pt cuts to use during jetfinding')
  parser.add_argument('--leadConstPtMatch', default='0.2', help='list of leading matched jet constituent pt cuts to use during jetfinding')
  parser.add_argument('--subConstPt', default='2.0', help='list of subleading hard jet constituent pt cuts to use during jetfinding')
  parser.add_argument('--subConstPtMatch', default='0.2', help='list of subleading matched jet constituent pt cuts to use during jetfinding')
  parser.add_argument('--leadR', default='0.4', help='list of jet radii to be used for leading jet')
  parser.add_argument('--leadJetPt', default='20.0', help='list of leading jet pt cuts to use during jetfinding')
  parser.add_argument('--subR', default='0.4', help='list of jet radii to be used for subleading jet')
  parser.add_argument('--subJetPt', default='10.0', help='list of subleading jet pt cuts to use during jetfinding')
  
  args = parser.parse_args()
  main( args )
