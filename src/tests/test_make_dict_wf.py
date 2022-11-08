import unittest
import sys
sys.path.insert(1, '/home/marinedjaffardjy/Documents/wf_features/src/parsing')
import parsing_snkmk as ps


file = "data/1.snakefile" #bioxfu
file2 = "data/3.snakefile" #davemcg

class MyTestCase(unittest.TestCase):

    def test_load(self):
        wf = ps.parse_wf(file)
        wf = ps.make_dict_wf(wf)
        wf_keys = ['rulecount', 'compilation', 'filename', 'rules_list', 'rules_names', 'input_names', 'output_names', 'extracted_toolnames']
        #all the keys are here
        for k in wf_keys :
            self.assertTrue(k in wf.keys(),msg="Cannot find key "+str(k))

        #the values associated with the keys are not empty
        #self.assert




if __name__ == '__main__':
    unittest.main()
