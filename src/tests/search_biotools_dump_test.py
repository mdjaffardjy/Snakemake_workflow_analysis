import unittest
from rdflib import ConjunctiveGraph
from rdflib.namespace import Namespace, RDF, RDFS

import sys
sys.path.insert(1, '/home/marinedjaffardjy/Documents/wf_features/src/parsing')
import search_biotools_dump as sb

tool_name = 'bedtool'
tool_dict = {'name': 'BEDTools', 'uri': 'https://bio.tools/bedtools' , 'topic': [{'uri': 'http://edamontology.org/topic_0622', 'term': 'Genomics'}], 'function': [{'operation': [{'uri': 'http://edamontology.org/operation_2429', 'term': 'Mapping'}, {'uri': 'http://edamontology.org/operation_2429', 'term': 'Cartography'}], 'input': [], 'output': []}]}

class MyTestCase(unittest.TestCase):

    def test_navigating(self):
        sorted_matches = sb.get_sorted_matches(tool_name)
        liste_tool_info = sb.extract_match_info(sorted_matches)
        try:
            assert (liste_tool_info[0] == tool_dict)
            print('Annotations trouvées !')

        except AssertionError:
            print('Erreur !')
            print('réf = '+str(tool_dict))
            print('réponse = '+str(liste_tool_info[0]))



if __name__ == '__main__':

    unittest.main()


