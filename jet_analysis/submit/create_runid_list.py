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
import numpy as np

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

  ## some paths
  execpath = os.getcwd()
  executable = './bin/data_qa/generate_runid_list'
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

      ## build our qsub execution string
      clargs = '--outDir=' + args.output + ' --input=' + files[i] + ' --id=' + str(i)
      clargs = clargs + ' --name=' + args.name + ' --towList=' + args.badTowers

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
  
  print "all jobs completed: adding root files"
  
  ## now we add all root files into a temporary file
  ## first get list
  output_files = [os.path.join(args.output, args.name + str(i) + ".root") for i in range(len(files))]
  root_tmp_file_name = "runid_tmp.root"
  proccommand = "hadd -f " + root_tmp_file_name + " " + " ".join(output_files)
  ret = subprocess.Popen(proccommand, shell=True)
  ret.wait()
  if ret.returncode != 0 :
    print "warning: hadd command failed, exiting"
    return

  infile = ROOT.TFile.Open(root_tmp_file_name, "READ")

  ## now I create a set to hold only unique runIDs
  runid_set = set()

  ## loop over the tree
  for event in infile.runid :
    runid_set.add(event.runid)

  ## fill a new TTree into the final output file
  root_output_file_name = "runid.root"
  outfile = ROOT.TFile.Open(root_output_file_name, "RECREATE")
  tree = ROOT.TTree("runid", "tree of runids")

  ## we will write an np array so that it gives us the address...
  runid = np.zeros(1, dtype=np.uint32)
  tree.Branch('runid', runid, 'runid/i')
  for id in runid_set :
    runid[0] = id
    tree.Fill()

  outfile.Write()
  outfile.Close()

  print "completed processing runids: cleaning up and exiting"
  infile.Close()
  os.remove(root_tmp_file_name)

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
  parser.add_argument('--output', default='out/tmp', help=' directory for output root files' )
  parser.add_argument('--badTowers', default='submit/empty_list.txt', help=' csv file containing towers to mask (needed to keep TStarJetPicoReader happy)')
  args = parser.parse_args()
  main( args )
