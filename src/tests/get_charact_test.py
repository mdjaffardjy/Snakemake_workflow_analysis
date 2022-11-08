import unittest
import sys
sys.path.insert(1, '/home/marinedjaffardjy/Documents/wf_features/src/parsing')
import parsing_snkmk as ps
import get_characteristics as t

class MyTestCase(unittest.TestCase):

    def test_make_df(self):
        path_800 = "/home/marinedjaffardjy/Documents/wf_features/data/inputs/data/"
        # workflow parsing
        wf_800 = ps.parse_wf_folder(path_800)
        dict_tool, total_extr = ps.make_dict_wf_with_biotools_dump(wf_800)
        df_800 = t.get_df_attributes(wf_800)
        print(df_800)
        assert(wf_800.len() == df_800.size)





if __name__ == '__main__':
    unittest.main()
