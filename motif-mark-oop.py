#!/usr/bin/env python
import argparse
import re
import math

def get_args():
  parser = argparse.ArgumentParser(description="Finds binding motifs around an exon and displays them in a picture.")
  parser.add_argument("-f", "--fasta", help="input fasta file", required=True)
  parser.add_argument("-m", "--motifFile", help="motif file", required=True)
  return parser.parse_args()

class MotifSearcher:
  def __init__(self, motifSeq):
    self.motifSeq = motifSeq.upper()
    self.motifLen = len(motifSeq)
    self.motifSearch = ""
    self._generalizeMotif()

  def _generalizeMotif(self):
    tempString = self.motifSeq.replace("U", "T")
    tempString = tempString.replace("W", "[AT]")
    tempString = tempString.replace("S", "[CG]")
    tempString = tempString.replace("M", "[AC]")
    tempString = tempString.replace("K", "[GT]")
    tempString = tempString.replace("R", "[AG]")
    tempString = tempString.replace("Y", "[CT]")
    tempString = tempString.replace("N", "[ACGT]")
    
    #a lookahead allows for overlapping motifs I think we want this
    #https://stackoverflow.com/questions/5616822/how-to-use-regex-to-find-all-overlapping-matches
    tempString = "(?=(" + tempString
    tempString += "))"
    self.motifSearch = tempString

  def findMotifs(self, searchSeq):
    matches = re.finditer(self.motifSearch, searchSeq.seq, re.IGNORECASE)
    print([m.start() for m in matches])

class Seq:
  def __init__(self, header, seq):
    self.header = header
    self.seq = seq

  def getBoundaries(self):
    exonStart = math.inf
    if "A" in self.seq:
      exonStart = min(self.seq.index("A"), exonStart)
    if "C" in self.seq:
      exonStart = min(self.seq.index("C"), exonStart)
    if "G" in self.seq:
      exonStart = min(self.seq.index("G"), exonStart)
    if "T" in self.seq:
      exonStart = min(self.seq.index("T"), exonStart)

    exonEnd = math.inf
    if "a" in self.seq[exonStart:]:
      exonEnd = min(self.seq[exonStart:].index("a"), exonEnd)
    if "c" in self.seq[exonStart:]:
      exonEnd = min(self.seq[exonStart:].index("c"), exonEnd)
    if "g" in self.seq[exonStart:]:
      exonEnd = min(self.seq[exonStart:].index("g"), exonEnd)
    if "t" in self.seq[exonStart:]:
      exonEnd = min(self.seq[exonStart:].index("t"), exonEnd)
    exonEnd += exonStart

    return (0, exonStart, exonEnd, len(self.seq))

def getHeadersAndSeqs(fileHandle):
  allSeqs = []
  header = ""
  seq = ""
  for line in fileHandle:
    line = line.strip()
    if line.startswith(">"):
      if header != "" and seq != "":
        allSeqs.append(Seq(header, seq))
        header = ""
        seq = ""
      header = line[1:]
    else:
      seq += line
  else:
    if header != "" and seq != "":
        allSeqs.append(Seq(header, seq))
  return allSeqs

def main():
  args=get_args()

  searchers = []
  with open(args.motifFile) as inFileMotif:
    for line in inFileMotif:
      searchers.append(MotifSearcher(line.strip()))

  with open(args.fasta) as inFileFasta:
    seqs = getHeadersAndSeqs(inFileFasta)

  for seq in seqs:
    print(seq.header)
    for searcher in searchers:
      print(searcher.motifSeq)
      searcher.findMotifs(seq)
    print("")

if __name__ == "__main__":
  main()