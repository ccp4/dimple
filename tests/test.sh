#!/bin/sh -e

test_dir=$(dirname $0)
export PYTHONPATH=$test_dir/..${PYTHONPATH:+:${PYTHONPATH}}

echo "Run dimple on artificial thaumatin example"
python -m dimple $test_dir/thaumatin.* $test_dir/th_out/

echo "output directory inspection"
python -m c4.workflow $test_dir/th_out
echo "job details directory inspection"
python -m c4.workflow $test_dir/th_out 1 5
echo "rerun 4th job"
python -c "import pickle; wf=pickle.load('$test_dir/th_out'); wf.jobs[3].run()"
