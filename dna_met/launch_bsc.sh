#!/usr/bin/env bash
#SBATCH --no-requeue
#SBATCH --mem 12G
#SBATCH -c 2
#SBATCH -t 168:00:00
#SBATCH --account=cli39  ## Our BSC account
 
set -e
set -u
 
_term() {
        echo "Caught SIGTERM signal!"
        kill -s SIGTERM $pid
        wait $pid
}
 
trap _term TERM
 
"$@" & pid=$!
 
echo "Waiting for ${pid}"
wait $pid
 
exit 0
