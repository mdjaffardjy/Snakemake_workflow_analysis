import unittest
import sys
sys.path.insert(1, '/home/marinedjaffardjy/Documents/wf_features/src/parsing')
import parsing_snkmk as ps
import get_characteristics as chr
import pandas as pd
from statistics import *

# file = "/home/marinedjaffardjy/Documents/wf_features/data/inputs/bon/select/davemcg/scEiaD/3.snakefile"
# file2 = "/home/marinedjaffardjy/Documents/wf_features/data/inputs/bon/select/bioxfu/RNA-Seq/1.snakefile"
#
# path_15 = "/home/marinedjaffardjy/Documents/wf_features/data/inputs/bon/select/"
# path_800 = "/home/marinedjaffardjy/Documents/wf_features/data/inputs/data/"

def quality_score (found, ref):
    #takes a list of tools and
    # returns a list of similar, different tools as well as a similarity score and a divergent score
    ref = ref.replace("{",'').replace("}",'').replace(' ','').split(",")
    ref = set(ref)
    found = set(found)
    if ref == "set()":
        ref = set()
    similar = []
    different = []

    if found:
        prec = len(ref.intersection(found)) / len(found)
    else :
        if ref :
            prec = 0
        else :
            prec = 1

    if ref:
        rapp = len(ref.intersection(found)) / len(ref)
    else :
        if found:
            rapp = 0
        else:
            rapp = 1

    #we add each element from the reference that is in found
    for element in ref:
        #print(element)
        if element in found:
            similar.append(element)
    #we add each element from the found elements not in the reference
    for element2 in found:
        #print(element2)
        if(element2 not in ref):
            different.append(element2)

    return similar, different, prec, rapp




class MyTestCase(unittest.TestCase):

    # def test_make_dict_wf(self):
    #     file = "/home/marinedjaffardjy/Documents/wf_features/data/inputs/bon/select/davemcg/scEiaD/3.snakefile"
    #
    #     wf = ps.parse_wf(file)
    #     wf = ps.make_dict_wf(wf)
    #     try:
    #         assert (wf['rulecount'] == 21)
    #         print('Bonne réponse !')
    #
    #     except AssertionError:
    #         print('Erreur !')

    # def test_tool_extraction(self):
    #     file = "/home/marinedjaffardjy/Documents/wf_features/data/inputs/bon/select/davemcg/scEiaD/3.snakefile"
    #     file2 = "/home/marinedjaffardjy/Documents/wf_features/data/inputs/bon/select/bioxfu/RNA-Seq/1.snakefile"
    #
    #     wf = ps.parse_wf(file)
    #     wf = ps.make_dict_wf(wf)
    #     wf2 = ps.parse_wf(file2)
    #     wf2 = ps.make_dict_wf(wf2)
    #     try:
    #         assert (wf['extracted_toolnames'] == [['kallisto', 'index'], ['bustools', 'sort'], ['bustools', 'whitelist'], ['bustools', 'correct'], ['bustools', 'count']])
    #         assert (wf2['extracted_toolnames'] == [['fastqc'], ['trimmomatic', 'PE'], ['fastqc'], ['hisat2'], ['samtools', 'view'], ['samtools', 'sort'], ['samtools', 'index'], ['qualimap', 'bamqc'], ['samtools', 'sort'], ['htseq-count'], ['bedtools', 'bamtobed'], ['bedtools', 'bed12tobed6'], ['samtools', 'view'], ['bedtools', 'genomecov'], ['igvtools', 'toTDF']])
    #         print('Bons outils extraits !')
    #     except AssertionError:
    #         print('Mauvais outils extraits !')
    #         print(wf['extracted_toolnames'] )
    #         print(wf2['extracted_toolnames'] )

    def test_tool_annotations(self):
        # file = "/home/marinedjaffardjy/Documents/wf_features/data/inputs/bon/select/davemcg/scEiaD/3.snakefile"
        # file2 = "/home/marinedjaffardjy/Documents/wf_features/data/inputs/bon/select/bioxfu/RNA-Seq/1.snakefile"
        # file3 = "/home/marinedjaffardjy/Documents/wf_features/data/inputs/data/bhattlab/kraken2_classification/1.snakefile"
        path = "/home/marinedjaffardjy/Documents/wf_features/data/inputs/bon/select/"

        # file = "/home/marinedjaffardjy/Documents/wf_features/data/inputs/bon/select/anndoan/class-18/1.snakefile"

        wf = ps.parse_wf_folder(path)
        # wf = [parse_wf(file)]
        dict_tool, total_extr = ps.make_dict_wf_with_biotools_dump(wf)
        #for el in wf :
         #   ps.write_json_in_txt(el, "/home/marinedjaffardjy/Documents/wf_features/data/outputs", "/home/marinedjaffardjy/Documents/wf_features/data/inputs/bon/select/")
        # print("extracted_toolnames")
        # print(wf[13]["extracted_toolnames"])
        # print(wf[13]["rules_list"])
        print("dict tools")
        print(len(dict_tool))
        #print(dict_tool)
        print("total extr")
        print(len(total_extr))
        #print(total_extr)
        print("dataframe of the wf attributes")
        df = chr.get_df_attributes(wf)

        #df.to_csv("/home/marinedjaffardjy/Documents/wf_features/data/outputs/out15.csv")
        # print(df)
        ref_file = "/home/marinedjaffardjy/Documents/wf_features/data/outputs/out15_annote.csv"
        df_ref = pd.read_csv(ref_file)
        #IMPORT THE CORRECTLY ANNOTATED ONE
        i = 0
        scorelist = []
        divlist = []
        difflist = []
        for liste_outils in df["tool names"] :
            similaire, different, score, div = quality_score(liste_outils, df_ref["tool names"][i] )
            scorelist.append(score)
            divlist.append(div)
            difflist.append(different)
            #print(df["name"][i])
            #print("score = "+str(score)+", div = "+str(div))
            try :
                assert(score == 1 and div == 0)
                print("Les outils trouvés correspondent pour tout le workflows.")
                print("réf : " + str(df_ref["tool names"][i]))
            except AssertionError :
                print("Les outils trouvés ne correspondent pas tous :")
                print("réf : "+str(df_ref["tool names"][i]))
                print("similaires : "+str(similaire))
                print("différents : "+str(different))
            i+=1
        print("precision = "+ str(scorelist) + "; mean : "+str(mean(scorelist)))
        print("rappel : "+str(divlist)+ "; mean : "+str(mean(divlist)))
        print("differences : "+ str(difflist))


if __name__ == '__main__':

    unittest.main()


