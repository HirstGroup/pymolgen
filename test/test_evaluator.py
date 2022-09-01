import sys,os
import filecmp

path = os.path.dirname(os.path.realpath(__file__))

def test_evaluator():

	os.system(f'python ../evaluator.py -a ../datasets/database1000/fragments1000.sdf -f ../datasets/database1000/fragments1000.txt -d ../datasets/database1000/frequencies1000.txt -p ../datasets/database1000/phenylisoxazole.sdf -x ../datasets/database1000/benzene.sdf ../datasets/database1000/benzene.sdf -r 20 21 -R 11 11 -o outputs/evaluator1000.sdf -n 10 --parent_mapping 15 0 16 0 --batch_size 2 --log outputs/evaluator1000.csv --filters --rules --rules_file rules1000.smi --mpo -5.0 --unique --mw_check -s 100 --no_time')

	compare = filecmp.cmp('outputs/evaluator1000.sdf', '../datasets/database1000/evaluator1000.sdf', shallow=False)
	assert compare is True

	compare = filecmp.cmp('outputs/evaluator1000.csv', '../datasets/database1000/evaluator1000.csv', shallow=False)
	assert compare is True

