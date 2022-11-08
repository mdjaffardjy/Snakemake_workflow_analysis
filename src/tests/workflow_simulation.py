import os
import unittest
from pathlib import Path
from io import StringIO
import sys
import subprocess

from snakemake import snakemake
from snakemake.exceptions import MissingInputException


class MyTestCase(unittest.TestCase):

    def test_launch_snakemake(self):

        #self.assertEqual(True, True)  # add assertion here
        base_path = Path(__file__).parent.parent.parent  # current directory
        #snakefile = str(base_path) + "/data/exome_pipeline/Snakefile"
        snakefile = str(base_path) + "/data/exome_pipeline/Snakefile"
        snakefile = str(base_path) + "/data/Example_Workflows/simple.snakefile"

        # static_file_path = str((base_path / "static/data/jsonldcontext.json").resolve())
        snakemake(snakefile, dryrun=True, printdag=True)

    def test_launch_snakemake(self):

        #self.assertEqual(True, True)  # add assertion here
        base_path = Path(__file__).parent.parent.parent  # current directory
        #snakefile = str(base_path) + "/data/exome_pipeline/Snakefile"
        snakefile = str(base_path) + "/data/exome_pipeline/Snakefile"
        snakefile = str(base_path) + "/data/Example_Workflows/simple.snakefile"

        # static_file_path = str((base_path / "static/data/jsonldcontext.json").resolve())
        snakemake(snakefile, dryrun=True, printdag=True)

    def dry_run_and_touch(sef, snakefile):
        base_path = Path(__file__).parent  # current directory

        command = f"snakemake -s {snakefile} --dryrun --rulegraph".split()
        out = subprocess.run(command, capture_output=True)

        retry = False
        if "MissingInputException" in str(out.stderr):
            lines = out.stderr.decode().split("\n")
            for i in range(3, len(lines)):
                missing_file = lines[i]
                if len(missing_file) > 0:
                    missing_file_path = Path(str(base_path) + os.sep + str(missing_file))
                    dir = Path(str(base_path) + os.sep + str(missing_file)).parent
                    if not dir.exists():
                        dir.mkdir(parents=True, exist_ok=True)
                    missing_file_path.touch()
                    print(f"CREATED missing input file {missing_file}")
                    retry = True
        else :
            print(out.stderr)
        return retry, out.stdout, out.stderr, out.returncode

    def test_missing_inputs(self):

        #self.assertEqual(True, True)  # add assertion here
        base_path = Path(__file__).parent.parent.parent  # data directory
        # snakefile = str(base_path) + "/Example_Workflows/more_complex.snakefile"
        # snakefile = str(base_path) + "/Example_Workflows/simple.snakefile"
        snakefile = str(base_path) + "/data/Example_Workflows/exome.snakefile"

        retry = True
        while (retry == True):
            retry, stdout, stderr, returncode = self.dry_run_and_touch(snakefile)

        print(stdout.decode())



if __name__ == '__main__':
    unittest.main()
