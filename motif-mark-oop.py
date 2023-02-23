#!/usr/bin/env python
import argparse
import re
import math
import cairo

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
    return [m.start() for m in matches]

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

    return [0, exonStart, exonEnd, len(self.seq)]

def getHeadersAndSeqs(fileHandle):
  allSeqs = [] #[Seq, Seq, ...]
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

def scalePositions(positionDict, boundaries, width, margin, bpWidth):
  #TODO: change scaling to all bps the same
  
  scaledBoundaries = [(entry * bpWidth) + margin for entry in boundaries]

  # scaledPostions = {} #motif:[(start,end), ...]
  # for motif, positions in positionDict.items():
  #   tempPositions = [((width - margin) - margin) * ((entry - boundaries[0]) / (boundaries[-1] - boundaries[0])) + margin for entry in positions]
  #   scaledMotifLen = ((width - margin) - margin) * ((len(motif) - boundaries[0]) / (boundaries[-1] - boundaries[0])) + margin
  #   scaledPostions[motif] = [(tempPositions[i], tempPositions[i] + scaledMotifLen) for i in range(len(tempPositions))]

  scaledPostions = {} #motif:[(start,end), ...]
  for motif, positions in positionDict.items():
    tempPositions = [(entry * bpWidth) + margin for entry in positions]
    scaledPostions[motif] = [(tempPositions[i], tempPositions[i] + (len(motif) * bpWidth)) for i in range(len(tempPositions))]
  
  print(scaledBoundaries)
  return scaledBoundaries, scaledPostions

def draw(seqs, searchers, fastaName):  
  colorPalette = [(213/255, 94/255, 0),
                  (0, 114/255, 178/255),
                  (240/255, 228/255, 66/255),
                  (0, 158/255, 115/255),
                  (204/255, 121/255, 167/255)]

  margin = 20
  width = 1000
  heightPerGroup = 100
  legendHeight = heightPerGroup
  height = (heightPerGroup * len(seqs)) + legendHeight
  bpWidth = (width - (margin * 2)) / max([len(seq.seq) for seq in seqs])
  
  fileBaseName = ".".join(fastaName.split(".")[:-1])
  surface = cairo.SVGSurface(f"{fileBaseName}.svg", width, height)
  context = cairo.Context(surface)
  context.select_font_face(
        "Arial", cairo.FONT_SLANT_NORMAL, cairo.FONT_WEIGHT_NORMAL)
  
  context.set_source_rgb(1,1,1)
  context.rectangle(0,0,width,height)
  context.fill()

  for i, searcher in enumerate(searchers):
    yCenter = legendHeight / 2
    xStart = margin + ((width - (margin * 2)) * (i / len(searchers)))

    print(xStart)
    
    context.move_to(xStart, yCenter)
    context.set_line_width(heightPerGroup / 2)
    color = colorPalette[i]
    context.set_source_rgba(color[0], color[1], color[2], 1)
    context.line_to(xStart + heightPerGroup / 2, yCenter)
    context.stroke()

    context.set_font_size(heightPerGroup / 4)
    context.move_to(xStart + heightPerGroup/2 + margin / 2, legendHeight / 1.5)
    context.show_text(searcher.motifSeq)

  for i, seq in enumerate(seqs):
    boundaries = seq.getBoundaries()
    print(boundaries)

    positionDict = {}
    for searcher in searchers:
      positionDict[searcher.motifSeq] = searcher.findMotifs(seq)
    boundaries, positionDict = scalePositions(positionDict, boundaries, width, margin, bpWidth)  #positionDict values are [(start, end),...], but scaled for printing
    print((boundaries[2] - boundaries[1]) / boundaries[3])

    yCenter = ((i + 1) * heightPerGroup) - (heightPerGroup / 2) + legendHeight
    context.set_source_rgba(0, 0, 0, 1)
    context.set_font_size(heightPerGroup / 3)
    
    context.move_to(margin, yCenter - heightPerGroup / 6)
    context.show_text(seq.header.split(" ")[0])
    
    #draw line for entire sequence
    context.move_to(boundaries[0], yCenter)
    context.set_line_width(1)
    context.line_to(boundaries[3], yCenter)
    context.stroke()
    
    #draw really thick line for exon (rectangles are boring)
    context.move_to(boundaries[1], yCenter)
    context.set_line_width(heightPerGroup / 2)
    context.line_to(boundaries[2], yCenter)
    context.stroke()

    #draw motifs
    context.set_line_width(heightPerGroup/4)
    for j, positions in enumerate(positionDict.values()):
      color = colorPalette[j]
      context.set_source_rgba(color[0], color[1], color[2], .75)
      for positionGroup in positions:
        context.move_to(positionGroup[0], yCenter)
        context.line_to(positionGroup[1], yCenter)
        context.stroke()
  
  surface.write_to_png(f"{fileBaseName}.png")

def main():
  args = get_args()

  searchers = []
  with open(args.motifFile) as inFileMotif:
    for line in inFileMotif:
      searchers.append(MotifSearcher(line.strip()))

  with open(args.fasta) as inFileFasta:
    seqs = getHeadersAndSeqs(inFileFasta)

  draw(seqs, searchers, args.fasta)
  

if __name__ == "__main__":
  main()