#!/bin/sh -e

test_dir=$(dirname $0)
export PYTHONPATH=$test_dir/..${PYTHONPATH:+:${PYTHONPATH}}

echo "Run dimple on artificial thaumatin example"
python -m dimple $test_dir/thaumatin.* $test_dir/th_out/

# the rest in not essential
# workflow is pickled as workflow.pickle in output dir, let's check it
echo "output directory inspection"
python -m c4.workflow $test_dir/th_out
echo "job details directory inspection"
python -m c4.workflow $test_dir/th_out 1 5
echo "rerun 4th job"
python -c "import pickle; pickle.load(open('$test_dir/th_out/workflow.pickle')).jobs[3].run()"
