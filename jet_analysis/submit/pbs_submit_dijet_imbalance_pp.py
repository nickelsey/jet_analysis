## pbs_submit_dijet_imbalance_auau.py

## convenience script to batch submit and monitor
## job status, automatically resubmitting jobs that
## do not produce the expected output files

from __future__ import print_function, division

import os
import argparse
import subprocess
import time
import random
import ROOT

def checkstatus(jobstatus) :
  loop = False
  for key in jobstatus :
    if jobstatus[key][1] != 2 :
      loop = True
      break
  return loop

def updatestatus(jobstatus, outdir, name) :
  print("Updating job status")
  print("Total: " + str(len(jobstatus)))
  
  ## get the qstat job listing
  proccommand = 'qstat | grep dx5412'
  proc = subprocess.Popen( proccommand, stdout=subprocess.PIPE, shell=True)
  qstat_result = proc.stdout.read()
  
  for key in jobstatus :
    
    index = jobstatus[key][0]
    
    ## if job is completed, we don't need to check again
    if jobstatus[key][1] == 2 :
      continue
    
    ## check if the job is still underway
    jobinprocess = qstat_result.find(name + str(index))
    if jobinprocess >= 0 :
      jobstatus[key][1] = 1
      continue
    
    ## if the job is not still underway,
    ## check to see if the job has completed properly
    ## if not, mark to resubmit
    if outdir.startswith('/') :
      filename = outdir + '/' + name + str(index) + '.root'
    else :
      filename = os.getcwd() + '/' + outdir + '/' + name + str(index) + '.root'
    
    if os.path.isfile(filename) :
      outputfile = ROOT.TFile(filename, "READ")
      if outputfile.IsZombie() :
        print("job " + str(index+1) + " of " + str(len(jobstatus)) + " complete: file is zombie, resubmit")
        jobstatus[key][1] = 0
        os.remove(filename)
      elif outputfile.IsOpen() :
        print("job " + str(index+1) + " of " + str(len(jobstatus)) + " complete: ROOT file healthy")
        print(filename)
        jobstatus[key][1] = 2
        outputfile.Close()
      else :
        print("job " + str(index+1) + " of " + str(len(jobstatus)) + " undefined file status, resubmit")
        jobstatus[key][1] = 0
    else :
      print("undefined status: job " + str(index+1) + " of " + str(len(jobstatus)) + " marked for submission")
      jobstatus[key][1] = 0

  return jobstatus


def main(args) :
  ## if there are no input files, exit
  files = args.strings
  if files is None :
    return
  
  ## get max number of jobs to be submitted to pbs
  ## at a time
  maxjobs = args.maxjobs

  ## select defaults
  embedTrig = ''
  if args.embedTriggers is not None :
    embedTrig = args.embedTriggers
  
  ## some paths
  execpath = os.getcwd()
  executable = './bin/dijet_imbalance/dijet_imbalance_pp'
  ## find the qwrap file
  qwrap = execpath + '/submit/qwrap.sh'
  
  ## we need to do our own book keeping
  ## for when  job is active and when it
  ## has completed successfully

  ## for pp we have the option of running systematic uncertainties as well
  ## this is a list of the different variations for each file
  systematic_variations = ['nom', 'ptow', 'mtow', 'ptrack', 'ntrack']

  if args.systematics == False :
    systematic_variations = ['nom']

  ## create the job list - its one job per input file, per systematic variation
  ## 0 not submitted, 1 for running, 2 for complete,
  ## and -1 for failed - resubmit
  counter = 0
  jobstatus = {}
  for file in files :
    for variation in systematic_variations :
      jobstatus[(file, variation)] = [counter, 0]
      counter = counter + 1
  
  ## count the number of qsub submission failures
  qsubfail = 0

  reader = ''
  embedReader = ''
  if args.readerSetting is not None :
    reader = args.readerSetting
  if args.embedReaderSetting is not None :
    embedReader = args.readerSetting
  
  while checkstatus(jobstatus) :
    
    ## update status of all jobs & output
    jobstatus = updatestatus(jobstatus, args.output, args.name)
    
    ## if we have completed all jobs, exit
    completed = 0
    for key in jobstatus :
      if jobstatus[key][1] == 2 :
        completed = completed + 1
    if len(jobstatus) == completed :
      break
    
    ## find the number of jobs still running via qstat
    ## if its at the maximum set jobsactive or jobsqueue,
    ## then pause
    jobsactive = 0
    for key in jobstatus :
      if jobstatus[key][1] == 1 :
        jobsactive = jobsactive + 1
    while jobsactive >= maxjobs :
      print("reached max number of active jobs: pausing")
      time.sleep(60)
      jobstatus = updatestatus(jobstatus, args.output, args.name)
      jobsactive = 0
      for key in jobstatus :
        if jobstatus[key][1] == 1 :
          jobsactive = jobsactive + 1
  
    ## now submit jobs up to maxjobs - jobsqueued
    njobs = maxjobs - jobsactive
    
    for key in jobstatus :
      if njobs == 0 :
        break
      if jobstatus[key][1] == 1 or jobstatus[key][1] == 2 :
        continue

      index = jobstatus[key][0]

      ## build log & error locations
      outstream = "log/" + args.name + str(index) + ".log"
      errstream = "log/" + args.name + str(index) + ".err"

      ## figure out which systematic variation we are running
      tow_sys = 0
      track_sys = 0
      outDir = args.output
      if key[1] == 'mtrack' :
        track_sys = -1
        outDir = outDir + '/tow_0_track_-1'
      if key[1] == 'ptrack' :
        track_sys = 1
        outDir = outDir + '/tow_0_track_1'
      if key[1] == 'mtow' :
        tow_sys = -1
        outDir = outDir + '/tow_-1_track_0'
      if key[1] == 'ptow' :
        tow_sys = 1
        outDir = outDir + '/tow_1_track_0'
      if key[1] == 'nom' :
        outDir = outDir + '/tow_0_track_0'
      
      ## build our qsub execution string
      clargs = '--outDir=' + outDir + ' --input=' + key[0] + ' --id=' + str(index)
      clargs = clargs + ' --name=' + args.name
      clargs = clargs + ' --runList=' + args.badRuns + ' --towList=' + args.badTowers + ' --triggers='
      clargs = clargs + args.triggers + ' --embedTriggers=' + embedTrig + ' --constEta=' + args.constEta
      clargs = clargs + ' --leadConstPt=' + args.leadConstPt + ' --subConstPt=' + args.subConstPt
      clargs = clargs + ' --leadConstPtMatch=' + args.leadConstPtMatch + ' --subConstPtMatch='
      clargs = clargs + args.subConstPtMatch + ' --leadR=' + args.leadR + ' --subR=' + args.subR
      clargs = clargs + ' --leadJetPt=' + args.leadJetPt + ' --subJetPt=' + args.subJetPt
      clargs = clargs + ' --readerSetting=' + reader + ' --embedReaderSetting=' + embedReader
      clargs = clargs + ' --towerUnc=' + str(tow_sys) + ' --trackUnc=' + str(track_sys)
      clargs = clargs + ' --efficFile=' + args.efficiencyFile + ' --embed=' + args.embedFile
      
      qsub = 'qsub -V -p ' + str(args.priority) + ' -l mem=' + str(args.mem) + 'GB -l nodes=' + str(args.nodes)
      qsub = qsub + ':ppn=' + str(args.ppn) + ' -q ' + str(args.queue) + ' -o ' + outstream
      qsub = qsub + ' -e ' + errstream + ' -N ' + args.name + str(index) + ' -- '
      qsub = qsub + qwrap + ' ' + execpath + ' ' + executable + ' ' + clargs
      print("submitting job: ")
      print(qsub)
      ret = subprocess.Popen( qsub, shell=True)
      ret.wait()
      if ret.returncode == 0 :
        jobstatus[key][1] = 1
      else :
        print("warning: qsub submission failed")
        qsubfail = qsubfail + 1
      njobs = njobs - 1
        
    if qsubfail > args.max_failures :
      print("qsub failure too many times - exiting")
      return
  
    ## wait one minute before rechecking
    print("finished round of submissions: pausing")
    time.sleep(60)
  
  print("all jobs completed: exiting")

if __name__ == "__main__":
  parser = argparse.ArgumentParser(description='Submit jobs via PBS & resubmit if necessary')
  parser.add_argument('strings', metavar='S', nargs='*', help=' input files ')
  parser.add_argument('--max_failures', default=200, help=' number of job submissions that can fail before forcing an exit')
  parser.add_argument('--name', default='job_', help=' job name (identifier for qsub)')
  parser.add_argument('--mem', type=int, default=2, help=' memory required per job [GB]')
  parser.add_argument('--nodes', type=int, default=1, help=' number of nodes required per job')
  parser.add_argument('--ppn', type=int, default=1, help=' number of processors per node required per job')
  parser.add_argument('--priority', type=int, default=0, help=' queue priority')
  parser.add_argument('--queue', default='erhiq', help=' queue to submit jobs to' )
  parser.add_argument('--maxjobs',type=int, default=100, help=' max number of jobs to have in running or queue states')
  parser.add_argument('--systematics', type=bool, default=True, help=' will run 4 systematic variations for error estimation')
  parser.add_argument('--output', default='out/post/tmp', help=' directory for output root files' )
  parser.add_argument('--embedFile', default='submit/y14_mb_file_list.txt', help=' file containing list of root files for embedding events')
  parser.add_argument('--efficiencyFile', default='submit/y14_effic_dca1.root', help=' root file containing run 14 efficiency curves')
  parser.add_argument('--badRuns', default='submit/y14_bad_run.txt', help=' csv file containing runs to mask')
  parser.add_argument('--badTowers', default='submit/y14_y6_bad_tower.txt', help=' csv file containing towers to mask')
  parser.add_argument('--triggers', default='ALL', help=' event triggers to consider: [y7, y10, y11, y14, y6pp, y9pp, y12pp] + [HT, MB, HT2, HT3, VPDMB30, VPDMB5, MBMON, ALL] (default "ALL": accept all events)')
  parser.add_argument('--embedTriggers', default='y14vpdmb30', help=' event triggers to consider for embedding data (see above for options)')
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

