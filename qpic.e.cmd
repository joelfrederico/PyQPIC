#!/bin/csh -f
#  qpic.e.cmd
#
#  UGE job for qpic.e built Tue Aug 12 12:23:48 PDT 2014
#
#  The following items pertain to this script
#  Use current working directory
#$ -cwd
#  input           = /dev/null
#  output          = /u/scratch/l/litos/work/test2/qpic.e.joblog
#$ -o /u/scratch/l/litos/work/test2/qpic.e.joblog.$JOB_ID
#  error           = Merged with joblog
#$ -j y
#  The following items pertain to the user program
#  user program    = /u/scratch/l/litos/work/test2/qpic.e
#  arguments       = 
#  program input   = Specified by user program
#  program output  = Specified by user program
#  Parallelism:  8-way parallel
#  Resources requested
#$ -pe shared 8
#$ -l h_data=4096M,h_rt=8:00:00,highp
# # #$ -q mori_msa.q
#
#  Name of application for log
#$ -v QQAPP=openmpi
#  Email address to notify
#$ -M litos@mail
#  Notify at beginning and end of job
#$ -m bea
#  Job is not rerunable
#$ -r n
#  Uncomment the next line to have your environment variables used by SGE
#$ -V
#
# Initialization for mpi parallel execution
#
  unalias *
  set qqversion = 
  set qqapp     = "openmpi parallel"
  set qqptasks  = 8
  set qqidir    = /u/scratch/l/litos/work/test2
  set qqjob     = qpic.e
  set qqodir    = /u/scratch/l/litos/work/test2
  cd     /u/scratch/l/litos/work/test2
  source /u/local/bin/qq.sge/qr.runtime
  if ($status != 0) exit (1)
#
  echo "UGE job for qpic.e built Tue Aug 12 12:23:48 PDT 2014"
  echo ""
  echo "  qpic.e directory:"
  echo "    "/u/scratch/l/litos/work/test2
  echo "  Submitted to UGE:"
  echo "    "$qqsubmit
  echo "  'scratch' directory (on each node):"
  echo "    $qqscratch"
  echo "  qpic.e 8-way parallel job configuration:"
  echo "    $qqconfig" | tr "\\" "\n"
#
  echo ""
  echo "qpic.e started on:   "` hostname -s `
  echo "qpic.e started at:   "` date `
  echo ""
#
# Run the user program
#

  source /u/local/Modules/default/init/modules.csh
  module load intel/11.1
  module load openmpi/1.4

  echo qpic.e -nt 8 "" \>\& qpic.e.output.$JOB_ID
  echo ""

  setenv PATH /u/local/bin:$PATH

  time $MPI_BIN/mpiexec -n 8 -machinefile $QQ_NODEFILE  \
         /u/scratch/l/litos/work/test2/qpic.e  >& /u/scratch/l/litos/work/test2/qpic.e.output.$JOB_ID 

  echo ""
  echo "qpic.e finished at:  "` date `

#
# Cleanup after mpi parallel execution
#
  source /u/local/bin/qq.sge/qr.runtime
#
  echo "-------- /u/scratch/l/litos/work/test2/qpic.e.joblog.$JOB_ID --------" >> /u/local/apps/queue.logs/openmpi.log.parallel
 if (`wc -l /u/scratch/l/litos/work/test2/qpic.e.joblog.$JOB_ID  | awk '{print $1}'` >= 1000) then
        head -50 /u/scratch/l/litos/work/test2/qpic.e.joblog.$JOB_ID >> /u/local/apps/queue.logs/openmpi.log.parallel
        echo " "  >> /u/local/apps/queue.logs/openmpi.log.parallel
        tail -10 /u/scratch/l/litos/work/test2/qpic.e.joblog.$JOB_ID >> /u/local/apps/queue.logs/openmpi.log.parallel
  else
        cat /u/scratch/l/litos/work/test2/qpic.e.joblog.$JOB_ID >> /u/local/apps/queue.logs/openmpi.log.parallel
  endif
#  cat            /u/scratch/l/litos/work/test2/qpic.e.joblog.$JOB_ID           >> /u/local/apps/queue.logs/openmpi.log.parallel
  exit (0)
