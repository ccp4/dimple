#!/bin/sh -e

test_dir=$(dirname $0)

# test local version rather than installed one
dimple() { $test_dir/../dimple.py "$@"; }
export PYTHONPATH=$test_dir/..${PYTHONPATH:+:${PYTHONPATH}}

echo "Run dimple on artificial thaumatin example"
dimple $test_dir/thaumatin.* $test_dir/th_out/

# the rest in not essential
# workflow is pickled as workflow.pickle in output dir, let's check it
echo
echo "-> output directory inspection"
dimple info $test_dir/th_out
echo
echo "-> job details directory inspection"
dimple info $test_dir/th_out 1 5
echo
echo "-> rerun 4th job"
dimple repeat $test_dir/th_out 4
echo
echo "-> the same in another way"
python -c "import pickle; pickle.load(open('$test_dir/th_out/workflow.pickle')).jobs[3].run()"
