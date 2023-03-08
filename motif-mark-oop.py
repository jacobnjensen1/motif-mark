#!/usr/bin/env python
import argparse
import re
import cairo

def get_args():
  """Gets user defined inputs."""
  parser = argparse.ArgumentParser(description="Finds binding motifs around an exon and displays them in a picture.")
  parser.add_argument("-f", "--fasta", help="input fasta file", required=True)
  parser.add_argument("-m", "--motifFile", help="motif file", required=True)
  return parser.parse_args()

class MotifSearcher:  
  """Holds a motif and finds that motif in sequences."""
  def __init__(self, motifSeq, motifIndex):
    """Makes a MotifSearcher object. Requires the motif sequence of course,
    and requires a motifIndex so that Drawer can choose a color and position the motif in the legend."""
    self.motifSeq = motifSeq.upper()
    self.motifIndex = motifIndex
    self.motifLen = len(motifSeq)
    self.motifSearch = ""
    self._generalizeMotif()

  def _generalizeMotif(self):
    """Replaces ambiguous bases and generates a regex from the input motif."""
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
    """Finds object's motif in a seq object's sequence"""
    matches = re.finditer(self.motifSearch, searchSeq.seq, re.IGNORECASE)
    return [m.start() for m in matches]

class Seq:
  """Stores all necessary information for a sequence and finds exon boundaries"""
  bpWidth = 0
  def __init__(self, header, seq, seqIndex):
    """Requires header(will be printed in entrirety) and seq for obvious reasons.
    Also requires seqIndex so Drawer knows where to put that sequence on the plot."""
    self.header = header
    self.seq = seq
    self.seqIndex = seqIndex

  def getBoundaries(self):
    """finds the exon boundaries and returns a list of [0,exonStart, exonEnd, <sequence length>]"""
    try:
      exonSpan = re.search(r"[agct]+([AGTC]+)[agct]+", self.seq).span(1)
      return [0, exonSpan[0], exonSpan[1] - 1, len(self.seq)]
    except:
      return [0,1,2,len(self.seq)] #return ridiculous numbers

class Drawer:
  """Drawer objects deal with the pycairo surface and context, and draw all shapes."""
  colorPalette = [(245/255, 84/255, 66/255),   # class variable. I want all objects to refer to this palette.
                  (194/255, 110/255, 183/255),
                  (45/255, 150/255, 44/255),
                  (67/255, 78/255, 209/255),
                  (210/255, 245/255, 66/255)]
  def __init__(self, heightPerGroup, width, margin, seqs, fastaName):
    """Builds a Drawer object, requires information about the desired size of the png output
    as well as all sequences to calculate scaling factors. Creates the pycairo surface and context."""
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

    self.fileBaseName = fileBaseName
    self.surface = surface
    self.context = context
    self.bpWidth = bpWidth
    self.legendHeight = legendHeight
    self.height = height
    self.width = width
    self.margin = margin
    self.heightPerGroup = heightPerGroup

  def saveOutput(self):
    """Saves the surface we've been building as a png"""
    self.surface.write_to_png(f"{self.fileBaseName}.png")

  def drawLegend(self, searchers):
    """Draws the legend, requires all motifSearcher objects"""
    for searcher in searchers:
      i = searcher.motifIndex
      yCenter = self.legendHeight / 2
      xStart = self.margin + ((self.width - (self.margin * 2)) * (i / len(searchers)))

      self.context.move_to(xStart, yCenter)
      self.context.set_line_width(self.heightPerGroup / 2)
      color = self.colorPalette[i]
      self.context.set_source_rgba(color[0], color[1], color[2], 1)
      self.context.line_to(xStart + self.heightPerGroup / 2, yCenter)
      self.context.stroke()

      self.context.set_font_size(self.heightPerGroup / 4)
      self.context.move_to(xStart + self.heightPerGroup/2 + self.margin / 2, self.legendHeight / 1.5)
      self.context.show_text(searcher.motifSeq)

  def _scalePosition(self, position):
    """Scales positions according to the already calculated width per base pair"""
    return (position * self.bpWidth) + self.margin

  def drawSequence(self, seq):    #first - introns and exons
    """Draws only the introns and exon of a specific sequence"""
    i = seq.seqIndex
    boundaries = seq.getBoundaries()
    boundaries = [self._scalePosition(entry) for entry in boundaries]

    yCenter = ((i + 1) * self.heightPerGroup) - (self.heightPerGroup / 2) + self.legendHeight
    self.context.set_source_rgba(0, 0, 0, 1)
    self.context.set_font_size(self.heightPerGroup / 3)
    
    #print first portion of header
    self.context.move_to(self.margin, yCenter - self.heightPerGroup / 6)
    self.context.show_text(seq.header.split(" ")[0])
    
    #draw line for entire sequence
    self.context.move_to(boundaries[0], yCenter)
    self.context.set_line_width(1)
    self.context.line_to(boundaries[3], yCenter)
    self.context.stroke()
    
    #draw really thick line for exon (rectangles are boring)
    self.context.move_to(boundaries[1], yCenter)
    self.context.set_line_width(self.heightPerGroup / 1.75)
    self.context.line_to(boundaries[2], yCenter)
    self.context.stroke()

  def drawMotifsOnSequence(self, seq, searchers):   #call this second, just motifs
    """Draws all motifs on a specific sequence"""
    testPositionsAll = []
    for searcher in searchers:
      testPositions = searcher.findMotifs(seq)
      for testPosition in testPositions:
        maxRequiredRow = 1
        #compare sets of ranges, check if there is overlap
        #if there is overlap, the new motif needs to be shown below the other one
        testPositionRange = {*range(testPosition, testPosition + searcher.motifLen + 1)}
        for testPositionAll in testPositionsAll:
          testPositionAllRange = {*range(testPositionAll[0], testPositionAll[1] + 1)}
          if len(testPositionRange.intersection(testPositionAllRange)) > 1:
            maxRequiredRow = max(maxRequiredRow, testPositionAll[3] + 1)
        testPositionsAll.append([testPosition, testPosition + searcher.motifLen, searcher.motifIndex, maxRequiredRow])
    testPositionsAll = [[self._scalePosition(entry[0]), self._scalePosition(entry[1]), entry[2], entry[3]] for entry in testPositionsAll]
    #holds [scaled start, scaled end, motif index, row index]
    numRows = max([entry[3] for entry in testPositionsAll])
    print(numRows)
    
    lineWidth = self.heightPerGroup / (numRows + 1)
    self.context.set_line_width(lineWidth)
    for testMotif in testPositionsAll:
      yCenter = ((testMotif[3] / (numRows + 1)) * self.heightPerGroup) + self.legendHeight + (self.heightPerGroup * seq.seqIndex)
      self.context.set_source_rgba(self.colorPalette[testMotif[2]][0], self.colorPalette[testMotif[2]][1],self.colorPalette[testMotif[2]][2], .85)
      self.context.move_to(testMotif[0], yCenter)
      self.context.line_to(testMotif[1], yCenter)
      self.context.stroke()

def getHeadersAndSeqs(fileHandle):
  """Reads through the fasta and returns a list of Seq objects"""
  allSeqs = [] #[Seq, Seq, ...]
  header = ""
  seq = ""
  counter = 0
  for line in fileHandle:
    line = line.strip()
    if line.startswith(">"):
      if header != "" and seq != "":
        allSeqs.append(Seq(header, seq, counter))
        counter += 1
        header = ""
        seq = ""
      header = line[1:]
    else:
      seq += line
  else:
    if header != "" and seq != "":
        allSeqs.append(Seq(header, seq, counter))
  return allSeqs

def main():
  """Runs primary functionality of motif-mark"""
  args = get_args()

  searchers = []
  with open(args.motifFile) as inFileMotif:
    for i, line in enumerate(inFileMotif):
      searchers.append(MotifSearcher(line.strip(), i))

  with open(args.fasta) as inFileFasta:
    seqs = getHeadersAndSeqs(inFileFasta)

  margin = 20
  width = 1000
  heightPerGroup = 100

  drawer = Drawer(heightPerGroup, width, margin, seqs, args.fasta)
  drawer.drawLegend(searchers)
  for seq in seqs:
    drawer.drawSequence(seq)
    drawer.drawMotifsOnSequence(seq, searchers)
  drawer.saveOutput()

if __name__ == "__main__":
  main()