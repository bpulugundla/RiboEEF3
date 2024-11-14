import os
import subprocess
import gzip
import argparse
from pathlib import Path

def parse_args():
    parser = argparse.ArgumentParser(
        description="A script for quality control of FASTQ sequences"
    )

    parser.add_argument("--rawdir", type=str, help="Raw FASTQ file directory")
    parser.add_argument("--names", type=str, help="Names of samples")
    parser.add_argument("--threshold", type=float, help="Output file name (default: 0.80)")
    parser.add_argument("--read_len_min", type=int, help="Minimum length of reads (default: 25)")
    parser.add_argument("--read_len_max", type=int, help="Maximum length of reads (default: 35)")

    return parser.parse_args()


def qualityFilter(args):
    PHREDDict = {
        "!": 9.999999e-01, "\"": 7.943282e-01, "#": 6.309573e-01, "$": 5.011872e-01, "%": 3.981072e-01,
        "&": 3.162278e-01, "\'": 2.511886e-01, "(": 1.995262e-01, ")": 1.584893e-01, "*": 1.258925e-01,
        "+": 1.000000e-01, ",": 7.943282e-02, "-": 6.309573e-02, ".": 5.011872e-02, "/": 3.981072e-02,
        "0": 3.162278e-02, "1": 2.511886e-02, "2": 1.995262e-02, "3": 1.584893e-02, "4": 1.258925e-02,
        "5": 1.000000e-02, "6": 7.943282e-03, "7": 6.309573e-03, "8": 5.011872e-03, "9": 3.981072e-03,
        ":": 3.162278e-03, ";": 2.511886e-03, "<": 1.995262e-03, "=": 1.584893e-03, ">": 1.258925e-03,
        "?": 1.000000e-03, "@": 7.943282e-04, "A": 6.309573e-04, "B": 5.011872e-04, "C": 3.981072e-04,
        "D": 3.162278e-04, "E": 2.511886e-04, "F": 1.995262e-04, "G": 1.584893e-04, "H": 1.258925e-04,
        "I": 1.000000e-04, "J": 7.943282e-05
    }

    basedir = str(Path(args.rawdir).resolve().parent)
    fdir = basedir + "/2-Filtered"
    frdir = basedir + "/2-Filtered/Reports"

    cmd = "mkdir -p " + fdir + " " + frdir
    os.system(cmd)
    LogFileName = frdir + "/Quality_filtering" + "_iv_log.txt"
    LOG_FILE = open(LogFileName, "wt")

    files = args.names.split()
    #files = ["ERR2660266"]
    for name in files:
        low_qual      = 0
        short_reads   = 0
        long_reads    = 0
        included_reads= 0

        sample = basedir + "/1-Trimmed/" + name + "_trimmed.fastq.gz"
        cmd = "echo $(zcat " + sample + " | wc -l)/4 | bc"
        length = int(subprocess.check_output(cmd, shell=True, text=True).strip())
        
        report = " {:16}: {:>12,}".format(name, length); print(report)
        LOG_FILE.write(report + "\n")

        File    = gzip.open(basedir + "/1-Trimmed/" + name + "_trimmed.fastq.gz", "rt")
        FileOut = open(basedir + "/2-Filtered/" + name + "_filtered.fastq", "w")

        for i in range(0, length):
            Identifier  = File.readline().rstrip("\n")
            Sequence    = File.readline().rstrip("\n")
            QIdentifier = File.readline().rstrip("\n")
            PHRED       = File.readline().rstrip("\n")
            Score       = 1.0
            Len         = len(PHRED)
           
            if Len < args.read_len_min:
                short_reads += 1

            elif Len > args.read_len_max:
                long_reads  += 1

            else:
                for IdxL in range(0,Len):
                    Score = Score * (1 - PHREDDict[PHRED[IdxL]])

                if (Score > args.threshold):
                    included_reads += 1
                    FileOut.write(Identifier + "\n" + Sequence + "\n" + QIdentifier + "\n" + PHRED + "\n")
                else:
                    low_qual += 1

        report  = " Reads len < {:>3} : {:>12,}\n".format(args.read_len_min, short_reads)
        report += " Reads len > {:>3} : {:>12,}\n".format(args.read_len_max, long_reads)
        report += " Quality   < {:>3} : {:>12,}\n".format(args.threshold, low_qual)
        report += " Reads left      : {:>12,}".format(included_reads)

        print(report, "\n"); LOG_FILE.write(report + "\n\n")

        File.close()
        FileOut.close()

    LOG_FILE.close()


def main():
    args = parse_args()
    qualityFilter(args)

if __name__ == "__main__":
    main()
