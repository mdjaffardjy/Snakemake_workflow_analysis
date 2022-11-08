import unittest
import logging
import sys

from src.r2p2.Snakemake_WF import Snakemake_WF
from src.r2p2.Nextflow_WF import Nextflow_WF
from src.r2p2.Factory import Factory, Implem



logging.basicConfig(
    level=logging.DEBUG,
    format="[%(asctime)s] {%(filename)s:%(lineno)d} %(levelname)-8s %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
)
LOGGER = logging.getLogger()
if not LOGGER.handlers:
    LOGGER.addHandler(logging.StreamHandler(sys.stdout))

class MyTestCase(unittest.TestCase):

    def test_several_instances(self):
        wf1 = Snakemake_WF("wf1", "/tmp/snakefile_1.snakefile")
        logging.info(str(wf1))

        wf2 = Nextflow_WF("wf2", "/tmp/2.nf")
        logging.info(str(wf2))

        wf_db = [wf1, wf2]
        for wf in wf_db:
            wf.parse()

    def test_factory(self):
        factory = Factory()
        wfs = []
        for i in range (0,6):
            wfs.append(factory.create_workflow(impl=Implem.SNAKEMAKE, name=f"N{i}", file_name=f"/tmp/F{i}"))
        for i in range (0,3):
            wfs.append(factory.create_workflow(impl=Implem.NEXTFLOW, name=f"N{i}", file_name=f"/tmp/F{i}"))

        for wf in wfs:
            print(wf)
            wf.get_steps()



if __name__ == '__main__':
    unittest.main()
