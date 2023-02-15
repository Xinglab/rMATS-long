#!/usr/bin/env python3

'''
Author: Robert Wang (Xing Lab)
Date: 2022.03.01

This is a script to enumerate all transcript structure differences
between any given pair of transcript isoforms. Alternative transcript
structure events can be classified into 7 simple categories:
* Exon skipping (SE)
* Mutually exclusive exons (MXE)
* Alternative 5'-splice site usage (A5SS)
* Alternative 3'-splice site usage (A3SS)
* Intron retention (RI)
* Alternative first exon usage (AFE)
* Alternative last exon usage (ALE)

Any alternative transcript structure event that fails to fall into any
of the seven categories listed above will, by default, get characterized
as complex (COMPLEX). It is possible to have combinations of alternative
transcript structure events for any given pair of transcript isoforms.
'''

# Load required libraries
import sys, argparse
from numpy import sort
from pandas import DataFrame, read_csv, concat
from networkx import all_simple_paths, get_node_attributes, DiGraph

def ClassifyBubble(bubble, attributes):
    # Extract paths within bubble object
    path1, path2 = bubble

    # Assign bubble to an event label based on path1, path2, and attributes
    if min(len(path1), len(path2)) == 2 and max(len(path1), len(path2)) == 4:
        '''
        Check if bubble involves exon skipping (SE):
            * bubble start node has attribute 'end'
            * bubble end node has attribute 'start'
        Check if bubble involves intron retention (RI):
            * bubble start node has attribute 'start'
            * bubble end node has attribute 'end'
        '''
        if attributes[path1[0]] == 'start' and attributes[path1[-1]] == 'end':
            # Bubble corresponds to intron retention (RI); return coordinates of retained intron
            return 'RI', ':'.join(list(map(str,sort(list(set(path1).symmetric_difference(set(path2))))+[1,-1])))
        elif attributes[path1[0]] == 'end' and attributes[path1[-1]] == 'start':
            # Bubble corresponds to exon skipping (SE); return coordinates of skipped exon
            return 'SE', ':'.join(list(map(str,sort(list(set(path1).symmetric_difference(set(path2)))))))
        else:
            # Bubble cannot be classified as either SE or RI; default to COMPLEX
            return 'COMPLEX', 'NA'

    elif min(len(path1), len(path2)) == 3 and max(len(path1), len(path2)) == 3:
        '''
        Check if bubble involves an alternative 5'-splice site (A5SS):
            * bubble start and end nodes have the same attribute 'start'
        Check if bubble involves an alternative 3'-splice site (A3SS):
            * bubble start and end nodes have the same attribute 'end'
        '''
        if attributes[path1[0]] == 'start' and attributes[path1[-1]] == 'start':
            # Bubble corresponds to an alternative 5'-splice site (A5SS); return coordinates of short and long exons
            return 'A5SS', ':'.join(list(map(str,sort(path1[:2])))) + ';' + ':'.join(list(map(str,sort(path2[:2]))))
        elif attributes[path1[0]] == 'end' and attributes[path1[-1]] == 'end':
            # Bubble corresponds to an alternative 3'-splice site (A3SS); return coordinates of short and long exons
            return 'A3SS', ':'.join(list(map(str,sort(path1[-2:])))) + ';' + ':'.join(list(map(str,sort(path2[-2:]))))
        else:
            # Bubble cannot be classified as either A5SS or A3SS; default to COMPLEX
            return 'COMPLEX', 'NA'

    elif min(len(path1), len(path2)) == 4 and max(len(path1), len(path2)) == 4:
        '''
        Check if bubble involves mutually exclusive exons:
            * bubble start node has attribute 'end'
            * bubble end node has attribute 'start'
            * internal exons are not overlapping
        Check if bubble involves an alternative first exon:
            * bubble start node is the root
            * bubble end node has attribute 'start'
            * first exons are not overlapping
        Check if bubble involves an alternative last exon:
            * bubble start node has attribute 'end'
            * bubble end node is the sink
            * last exons are not overlapping
        '''
        if attributes[path1[0]] == 'end' and attributes[path1[-1]] == 'start' and not CheckOverlap(sort(path1[1:3]), sort(path2[1:3])):
            # Bubble involves mutually exclusive exons; return coordinates of these mutually exclusive exons
            return 'MXE', ':'.join(list(map(str,sort(path1[1:3])))) + ';' + ':'.join(list(map(str,sort(path2[1:3]))))
        elif attributes[path1[0]] == 'root' and attributes[path1[-1]] == 'start' and not CheckOverlap(sort(path1[1:3]), sort(path2[1:3])):
            # Bubble involves an alternative first exon; return coordinates of these first exons
            return 'AFE', ':'.join(list(map(str,sort(path1[1:3])))) + ';' + ':'.join(list(map(str,sort(path2[1:3]))))
        elif attributes[path1[0]] == 'end' and attributes[path1[-1]] == 'sink' and not CheckOverlap(sort(path1[1:3]), sort(path2[1:3])):
            # Bubble involves an alternative last exon; return coordinates of these last exons
            return 'ALE', ':'.join(list(map(str,sort(path1[1:3])))) + ';' + ':'.join(list(map(str,sort(path2[1:3]))))
        else:
            # Bubble cannot be classified as MXE, AFE, or ALE; default to COMPLEX
            return 'COMPLEX', 'NA'
    else:
        # Bubble cannot be classified into a simple event; default to COMPLEX
        return 'COMPLEX', 'NA'

def FindBubbles(spliceGraph):
    '''
    Algorithm:
        * Set 'root' node as current node
        * If the current node has an indegree of two:
            * We have reached the end of a bubble; enumerate all simple paths between the bubble start and bubble end
        * If the current node has outdegree of two:
            * Mark current node as the start of a bubble
            * Pick a child at random and set it to current node
        * Else if current node has outdegree of one:
            * Set child to current node
        * Else:
            * We have reached the root node; terminate procedure
    '''

    # Create an empty list to store all bubbles found in the splice graph
    bubbles = list()

    # Instantiate currentNode and bubbleStart to 'root'
    currentNode, bubbleStart = 'root', 'root'

    # Repeat procedure until bubbleStart is set to 'sink'
    while not bubbleStart == 'sink':
        # Check if currentNode has an indegree of two
        if spliceGraph.in_degree(currentNode) == 2:
            # We have reached the end of a bubble; enumerate all simple paths between bubbleStart and currentNode
            bubbles.append(list(all_simple_paths(spliceGraph, source=bubbleStart, target=currentNode)))

        if spliceGraph.out_degree(currentNode) == 2:
            # We have reached the beginning of a bubble
            bubbleStart = currentNode

            # Select first child of bubbleStart and set it as currentNode
            currentNode = list(spliceGraph.successors(currentNode))[0]

        elif spliceGraph.out_degree(currentNode) == 1:
            # Set child to currentNode
            currentNode = list(spliceGraph.successors(currentNode))[0]

        else:
            # We have reached the root node; set bubbleStart to currentNode, which is 'root'
            bubbleStart = currentNode

    return bubbles

def BuildSpliceGraph(inputDF):
    # Establish an empty directed splice graph
    spliceGraph = DiGraph()

    # Establish a list of nodes for spliceGraph
    startNodes = [(V, {'attribute': 'start'}) for V in set(inputDF['exonStart'].sum())]
    endNodes = [(V, {'attribute': 'end'}) for V in set(inputDF['exonEnd'].sum())]
    auxNodes = [('sink', {'attribute': 'sink'}), ('root', {'attribute': 'root'})]

    # Establish a list of edges for spliceGraph
    exonEdges = list(set(sum([list(zip(inputDF.iloc[idx]['exonStart'], inputDF.iloc[idx]['exonEnd'])) for idx in range(inputDF.shape[0])], [])))
    intronEdges = list(set(sum([list(zip(inputDF.iloc[idx]['exonEnd'][:-1], inputDF.iloc[idx]['exonStart'][1:])) for idx in range(inputDF.shape[0])], [])))
    auxEdges = list(set(sum([[('root', inputDF.iloc[idx]['exonStart'][0]), (inputDF.iloc[idx]['exonEnd'][-1], 'sink')] for idx in range(inputDF.shape[0])], [])))

    # Add nodes and edges to empty spliceGraph
    spliceGraph.add_nodes_from(startNodes + endNodes + auxNodes)
    spliceGraph.add_edges_from(exonEdges + intronEdges + auxEdges)

    return spliceGraph

def CheckOverlap(tuple1, tuple2):
    # Checks if tuple1 overlaps with tuple2
    return max(tuple1[0], tuple2[0]) <= min(tuple1[1], tuple2[1])

def SmoothEnds(inputDF):
    '''
    Transcript ends cannot be accurately determined using RNA sequencing.
    Differences in transcript start and end points will, however, manifest as
    bubbles in our bubble-search algorithm. If the terminal exons of the two
    transcript isoforms in inputDF are overlapping, then we will 'smooth out'
    the transcript ends by selecting a common coordinate value. By consequence,
    this obviates the capacity to detect all bona fide alternative promoter
    or alternative polyadenylation events.
    '''

    # Create a copy of inputDF
    copyDF = inputDF.copy()

    # Extract coordinates of terminal exons for both transcript isoforms
    firstExon1, firstExon2 = (inputDF['exonStart'].iloc[0][0], inputDF['exonEnd'].iloc[0][0]), (inputDF['exonStart'].iloc[1][0], inputDF['exonEnd'].iloc[1][0])
    lastExon1, lastExon2 = (inputDF['exonStart'].iloc[0][-1], inputDF['exonEnd'].iloc[0][-1]), (inputDF['exonStart'].iloc[1][-1], inputDF['exonEnd'].iloc[1][-1])

    # Adjust transcript start positions; if firstExon1 is identical to firstExon2, do nothing
    if firstExon1 != firstExon2:
        # Check if firstExon1 overlaps with firstExon2 (first sort the tuple coordinates)
        # If there's no overlap, then do nothing
        if CheckOverlap(sort(firstExon1), sort(firstExon2)):
            # Pick a common start coordinate for both first exons (default to start for firstExon1)
            copyDF['exonStart'].iloc[1][0] = copyDF['exonStart'].iloc[0][0]

    # Adjust transcript end positions; if lastExon1 is identical to lastExon2, do nothing
    if lastExon1 != lastExon2:
        # Check if lastExon1 overlaps with lastExon2 (first sort the tuple coordinates)
        # If there's no overlap, then do nothing
        if CheckOverlap(sort(lastExon1), sort(lastExon2)):
            # Pick a common end coordinate for both last exons (default to end for lastExon1)
            copyDF['exonEnd'].iloc[1][-1] = copyDF['exonEnd'].iloc[0][-1]

    return copyDF

def ParseGTF(infile):
    # Read infile as a pandas dataframe
    gtfDF = read_csv(infile, sep='\t', header=None)

    # Pull out transcript annotations in gtfDF
    transcriptDF = gtfDF[gtfDF[2] == 'transcript']

    # Double check that the GTF file only has annotations for two transcript isoforms
    if transcriptDF.shape[0] == 2:
        # Isolate transcript IDs for each transcript
        transcriptID = [[x for x in anno.split(';') if 'transcript_id' in x][0].split('\"')[1] for anno in transcriptDF[8]]

        # Establish an empty pandas dataframe for output
        inputDF = DataFrame(columns = ['transcript_ID', 'chrom', 'strand', 'exonStart', 'exonEnd'])

        # For each transcript in transcriptID, record exon start and end coordinates
        for transcript in transcriptID:
            # Determine strand associated with transcript
            strand = transcriptDF[transcriptDF[8].str.contains(transcript)][6].item()

            # Determine chromosome associated with transcript
            chrom = transcriptDF[transcriptDF[8].str.contains(transcript)][0].item()

            # Pull out rows of gtfDF corresponding to exons of given transcript
            exonDF = gtfDF[(gtfDF[2] == 'exon') & (gtfDF[8].str.contains(transcript))]

            # Establish sorted lists for exon start and end coordinates based on strand
            if strand == '+':
                exonStart, exonEnd = sort(exonDF[3]), sort(exonDF[4])
            else:
                exonStart, exonEnd = -sort(-exonDF[4]), -sort(-exonDF[3])

            # Update inputDF with information about current transcript
            inputDF = concat([
                inputDF,
                DataFrame([[transcript, chrom, strand, exonStart.tolist(), exonEnd.tolist()]],
                          columns = ['transcript_ID', 'chrom', 'strand', 'exonStart', 'exonEnd'])
            ])

        # Smooth out transcript ends if terminal exons of both transcript isoforms are overlapping
        inputDF = SmoothEnds(inputDF)

        return inputDF

    else:
        sys.exit('ERROR: Input GTF file does not contain two transcript isoforms!')

def main():
    moduleSummary = 'This is a script to enumerate all transcript structure differences between a pair of transcript isoforms'
    parser = argparse.ArgumentParser(description=moduleSummary)

    # Add arguments
    parser.add_argument('-i', metavar='/path/to/input/GTF', required=True,
        help='path to GTF file describing structures of two transcript isoforms')
    parser.add_argument('-o', metavar='/path/to/output/file', required=True,
        help='path to output file')

    # Parse command-line arguments
    args = parser.parse_args()
    infile, outfile = args.i, args.o

    # Parse input GTF of transcript isoforms
    inputDF = ParseGTF(infile)

    # ASSUME: Both transcript isoforms are on the same chromosome and strand
    chrom, strand = inputDF['chrom'].unique()[0], inputDF['strand'].unique()[0]
    chrom = str(chrom)  # ensure chrom is parsed as a string

    # Build a splice graph
    spliceGraph = BuildSpliceGraph(inputDF)

    # Generate dictionary of attributes associated with each node in spliceGraph
    attributes = get_node_attributes(spliceGraph, 'attribute')

    # Find all bubbles in spliceGraph
    bubbles = FindBubbles(spliceGraph)

    # Establish an empty pandas dataframe to describe each of the bubbles
    outputDF = DataFrame(columns = ['transcript1', 'transcript2', 'event', 'coordinates'])

    # Iterate over bubbles identified in spliceGraph
    for bubble in bubbles:
        # For each bubble, identify the underlying event and report coordinates
        # Note: COMPLEX events will not have reported coordinates
        event, coord = ClassifyBubble(bubble, attributes)

        # Update coord with chrom and strand information (only if coord is not 'NA')
        if coord != 'NA':
            coord = ';'.join([chrom + ':' + x + ':' + strand for x in coord.split(';')])

        # Update outputDF with information about each bubble
        outputDF = concat([
            outputDF,
            DataFrame([[inputDF['transcript_ID'].iloc[0], inputDF['transcript_ID'].iloc[1], event, coord]],
                      columns = ['transcript1', 'transcript2', 'event', 'coordinates'])
        ])

    # Print contents of outputDF to outfile
    outputDF.to_csv(outfile, sep='\t', index=False)

if __name__ == '__main__':
    main()
