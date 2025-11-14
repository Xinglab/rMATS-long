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
import sys
import argparse
from numpy import sort
from pandas import DataFrame, read_csv, concat
from networkx import all_simple_paths, DiGraph

def GetCoordStringForPath(path):
    pairs = list()
    start = None
    for node in path:
        attribute = node.attribute
        if attribute in ['root', 'sink']:
            if start is not None:
                pairs.append((start, start))
                start = None

            pairs.append((attribute, attribute))
            continue

        if attribute == 'start':
            if start is not None:
                pairs.append((start, start))

            start = node.value
            continue

        if attribute == 'end':
            if start is None:
                pairs.append((node.value, node.value))
                continue

            pairs.append((start, node.value))
            start = None

    if start is not None:
        pairs.append((start, start))

    parts = list()
    for pair in pairs:
        part = '{}:{}'.format(pair[0], pair[1])
        parts.append(part)

    return ';'.join(parts)

def GetComplexCoordString(path1, path2):
    path1_len = len(path1)
    path2_len = len(path2)
    if path1_len > path2_len:
        longer = path1
        shorter = path2
    else:
        longer = path2
        shorter = path1

    shorter_string = GetCoordStringForPath(shorter)
    longer_string = GetCoordStringForPath(longer)
    return ';'.join([shorter_string, longer_string])

def get_ri_coord_string(path1, path2):
    set_diff = set(path1).symmetric_difference(set(path2))
    sorted_coords = sort([node.value for node in set_diff])
    # adjust coords from the exon into the intron
    sorted_coords[0] += 1
    sorted_coords[1] -= 1
    return ':'.join([str(x) for x in sorted_coords])

def get_se_coord_string(path1, path2):
    set_diff = set(path1).symmetric_difference(set(path2))
    sorted_coords = sort([node.value for node in set_diff])
    return ':'.join([str(x) for x in sorted_coords])

def get_a5ss_coord_string(path1, path2):
    sorted_coords_1 = sort([node.value for node in path1[:2]])
    sorted_coords_2 = sort([node.value for node in path2[:2]])
    first_string = ':'.join([str(x) for x in sorted_coords_1])
    second_string = ':'.join([str(x) for x in sorted_coords_2])
    return '{};{}'.format(first_string, second_string)

def get_a3ss_coord_string(path1, path2):
    sorted_coords_1 = sort([node.value for node in path1[-2:]])
    sorted_coords_2 = sort([node.value for node in path2[-2:]])
    first_string = ':'.join([str(x) for x in sorted_coords_1])
    second_string = ':'.join([str(x) for x in sorted_coords_2])
    return '{};{}'.format(first_string, second_string)

def get_mxe_coord_string(path1, path2):
    sorted_coords_1 = sort([node.value for node in path1[1:3]])
    sorted_coords_2 = sort([node.value for node in path2[1:3]])
    first_string = ':'.join([str(x) for x in sorted_coords_1])
    second_string = ':'.join([str(x) for x in sorted_coords_2])
    return '{};{}'.format(first_string, second_string)

def get_afe_coord_string(path1, path2):
    return get_mxe_coord_string(path1, path2)

def get_ale_coord_string(path1, path2):
    return get_mxe_coord_string(path1, path2)

def ClassifyBubble(bubble):
    # Extract paths within bubble object
    path1, path2 = bubble

    path1_len = len(path1)
    path2_len = len(path2)
    min_path_len = min(path1_len, path2_len)
    max_path_len = max(path1_len, path2_len)
    # Assign bubble to an event label based on path1, path2, and attributes
    if min_path_len == 2 and max_path_len == 4:
        '''
        Check if bubble involves exon skipping (SE):
            * bubble start node has attribute 'end'
            * bubble end node has attribute 'start'
        Check if bubble involves intron retention (RI):
            * bubble start node has attribute 'start'
            * bubble end node has attribute 'end'
        '''
        if path1[0].attribute == 'start' and path1[-1].attribute == 'end':
            # Bubble corresponds to intron retention (RI); return coordinates of retained intron
            return 'RI', get_ri_coord_string(path1, path2)
        elif path1[0].attribute == 'end' and path1[-1].attribute == 'start':
            # Bubble corresponds to exon skipping (SE); return coordinates of skipped exon
            return 'SE', get_se_coord_string(path1, path2)
        else:
            # Bubble cannot be classified as either SE or RI; default to COMPLEX
            return 'COMPLEX', GetComplexCoordString(path1, path2)

    elif min_path_len == 3 and max_path_len == 3:
        '''
        Check if bubble involves an alternative 5'-splice site (A5SS):
            * bubble start and end nodes have the same attribute 'start'
        Check if bubble involves an alternative 3'-splice site (A3SS):
            * bubble start and end nodes have the same attribute 'end'
        '''
        if path1[0].attribute == 'start' and path1[-1].attribute == 'start':
            # Bubble corresponds to an alternative 5'-splice site (A5SS); return coordinates of short and long exons
            return 'A5SS', get_a5ss_coord_string(path1, path2)
        elif path1[0].attribute == 'end' and path1[-1].attribute == 'end':
            # Bubble corresponds to an alternative 3'-splice site (A3SS); return coordinates of short and long exons
            return 'A3SS', get_a3ss_coord_string(path1, path2)
        else:
            # Bubble cannot be classified as either A5SS or A3SS; default to COMPLEX
            return 'COMPLEX', GetComplexCoordString(path1, path2)

    elif min_path_len == 4 and max_path_len == 4:
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
        path_1_coords_to_check = sort([node.value for node in path1[1:3]])
        path_2_coords_to_check = sort([node.value for node in path2[1:3]])
        has_overlap = CheckOverlap(path_1_coords_to_check, path_2_coords_to_check)
        if path1[0].attribute == 'end' and path1[-1].attribute == 'start' and not has_overlap:
            # Bubble involves mutually exclusive exons; return coordinates of these mutually exclusive exons
            return 'MXE', get_mxe_coord_string(path1, path2)
        elif path1[0].attribute == 'root' and path1[-1].attribute == 'start' and not has_overlap:
            # Bubble involves an alternative first exon; return coordinates of these first exons
            return 'AFE', get_afe_coord_string(path1, path2)
        elif path1[0].attribute == 'end' and path1[-1].attribute == 'sink' and not has_overlap:
            # Bubble involves an alternative last exon; return coordinates of these last exons
            return 'ALE', get_ale_coord_string(path1, path2)
        else:
            # Bubble cannot be classified as MXE, AFE, or ALE; default to COMPLEX
            return 'COMPLEX', GetComplexCoordString(path1, path2)
    else:
        # Bubble cannot be classified into a simple event; default to COMPLEX
        return 'COMPLEX', GetComplexCoordString(path1, path2)

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
            * We have reached the sink node; terminate procedure
    '''

    # Create an empty list to store all bubbles found in the splice graph
    bubbles = list()

    # Instantiate currentNode and bubbleStart to 'root'
    root_node = Node('root', 'root')
    currentNode = root_node
    bubbleStart = root_node

    # Repeat procedure until bubbleStart is set to 'sink'
    while bubbleStart.value != 'sink':
        if spliceGraph.in_degree(currentNode) == 2:
            # We have reached the end of a bubble;
            # enumerate all simple paths between bubbleStart and currentNode.
            all_paths = list(all_simple_paths(spliceGraph, source=bubbleStart, target=currentNode))
            # The nodes can be an integer coordinate or a string like root or sink.
            # Consider everything as a string for sorting here.
            # The purpose of the sort is to ensure that each bubble has a deterministic
            # order. Then if the same bubble is found when running for different pairs of
            # isoforms, the output event coordinates will be the same.
            all_paths.sort(key=lambda path: [str(node.value) for node in path])
            bubbles.append(all_paths)

        if spliceGraph.out_degree(currentNode) == 2:
            # We have reached the beginning of a bubble
            bubbleStart = currentNode

            # Select first child of bubbleStart and set it as currentNode
            currentNode = list(spliceGraph.successors(currentNode))[0]

        elif spliceGraph.out_degree(currentNode) == 1:
            # Set child to currentNode
            currentNode = list(spliceGraph.successors(currentNode))[0]

        else:
            # We have reached the sink node; set bubbleStart to currentNode, which is 'sink'
            bubbleStart = currentNode

    return bubbles


class Node:
    def __init__(self, value, attribute):
        self.value = value
        self.attribute = attribute

    def __hash__(self):
        return hash((self.value, self.attribute))

    def __eq__(self, other):
        return (self.value, self.attribute) == (other.value, other.attribute)


def BuildSpliceGraph(inputDF):
    # Establish an empty directed splice graph
    spliceGraph = DiGraph()

    nodes = list()
    root_node = Node('root', 'root')
    nodes.append(root_node)
    sink_node = Node('sink', 'sink')
    nodes.append(sink_node)
    # The exonStart column in each row is a list of coordinates.
    # .sum() concats the lists
    start_coords = inputDF['exonStart'].sum()
    for start_coord in set(start_coords):
        nodes.append(Node(start_coord, 'start'))

    end_coords = inputDF['exonEnd'].sum()
    for end_coord in set(end_coords):
        nodes.append(Node(end_coord, 'end'))

    edges = set()
    for transcript_i in [0, 1]:
        exon_starts = inputDF.iloc[transcript_i]['exonStart']
        exon_ends = inputDF.iloc[transcript_i]['exonEnd']
        num_exons = len(exon_starts)
        root_edge = (root_node, Node(exon_starts[0], 'start'))
        edges.add(root_edge)
        sink_edge = (Node(exon_ends[-1], 'end'), sink_node)
        edges.add(sink_edge)
        for exon_i, exon_start in enumerate(exon_starts):
            exon_end = exon_ends[exon_i]
            exon_start_node = Node(exon_start, 'start')
            exon_end_node = Node(exon_end, 'end')
            exon_edge = (exon_start_node, exon_end_node)
            edges.add(exon_edge)
            next_exon_i = exon_i + 1
            if next_exon_i < num_exons:
                next_exon_start = exon_starts[next_exon_i]
                next_exon_start_node = Node(next_exon_start, 'start')
                intron_edge = (exon_end_node, next_exon_start_node)
                edges.add(intron_edge)

    # Add nodes and edges to empty spliceGraph
    spliceGraph.add_nodes_from(nodes)
    spliceGraph.add_edges_from(edges)

    return spliceGraph

def CheckOverlap(tuple1, tuple2):
    # Checks if tuple1 overlaps with tuple2
    return max(tuple1[0], tuple2[0]) <= min(tuple1[1], tuple2[1])

def SmoothEnds(inputDF, is_minus_strand):
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
    tran_1_first_exon_start = inputDF['exonStart'].iloc[0][0]
    tran_1_first_exon_end = inputDF['exonEnd'].iloc[0][0]
    tran_2_first_exon_start = inputDF['exonStart'].iloc[1][0]
    tran_2_first_exon_end = inputDF['exonEnd'].iloc[1][0]
    tran_1_last_exon_start = inputDF['exonStart'].iloc[0][-1]
    tran_1_last_exon_end = inputDF['exonEnd'].iloc[0][-1]
    tran_2_last_exon_start = inputDF['exonStart'].iloc[1][-1]
    tran_2_last_exon_end = inputDF['exonEnd'].iloc[1][-1]

    tran_1_first_exon = (tran_1_first_exon_start, tran_1_first_exon_end)
    tran_2_first_exon = (tran_2_first_exon_start, tran_2_first_exon_end)
    tran_1_last_exon = (tran_1_last_exon_start, tran_1_last_exon_end)
    tran_2_last_exon = (tran_2_last_exon_start, tran_2_last_exon_end)
    # Tuples are sorted for the overlap check to handle minus strand
    # coordinates being reversed.
    # Adjust transcript start positions if there's overlap.
    if CheckOverlap(sort(tran_1_first_exon), sort(tran_2_first_exon)):
        # Pick a common start coordinate for both first exons.
        # Use the coordinate that increases the exon length.
        if is_minus_strand:
            higher_start = max(tran_1_first_exon_start, tran_2_first_exon_start)
            copyDF['exonStart'].iloc[0][0] = higher_start
            copyDF['exonStart'].iloc[1][0] = higher_start
        else:
            lower_start = min(tran_1_first_exon_start, tran_2_first_exon_start)
            copyDF['exonStart'].iloc[0][0] = lower_start
            copyDF['exonStart'].iloc[1][0] = lower_start

    if CheckOverlap(sort(tran_1_last_exon), sort(tran_2_last_exon)):
        # Pick a common end coordinate for both last exons.
        # Use the coordinate that increases the exon length.
        if is_minus_strand:
            lower_end = min(tran_1_last_exon_end, tran_2_last_exon_end)
            copyDF['exonEnd'].iloc[0][-1] = lower_end
            copyDF['exonEnd'].iloc[1][-1] = lower_end
        else:
            higher_end = max(tran_1_last_exon_end, tran_2_last_exon_end)
            copyDF['exonEnd'].iloc[0][-1] = higher_end
            copyDF['exonEnd'].iloc[1][-1] = higher_end

    return copyDF

def ParseAttributeString(attributes_str):
    attributes = dict()
    attribute_pairs = attributes_str.split(';')
    for attribute_pair in attribute_pairs:
        attribute_pair = attribute_pair.strip()
        first_space = attribute_pair.find(' ')
        if first_space <= 0:
            continue

        key = attribute_pair[:first_space]
        value = attribute_pair[first_space + 1:]
        if value[0] == '"' and value[-1] == '"':
            # remove quotes
            value = value[1:-1]

        existing_attributes = attributes.get(key)
        if existing_attributes:
            if isinstance(existing_attributes, list):
                existing_attributes.append(value)
            else:
                attributes[key] = [existing_attributes, value]
        else:
            attributes[key] = value

    return attributes

def getTranscriptIdFromRow(row):
    attribute_string = row[8]
    attributes = ParseAttributeString(attribute_string)
    transcript_id = attributes.get('transcript_id')
    return transcript_id

def ParseGTF(infile):
    # Read infile as a pandas dataframe
    gtfDF = read_csv(infile, sep='\t', header=None)
    allTranscriptID = gtfDF.apply(getTranscriptIdFromRow, axis=1)
    gtfDF.insert(9, 9, allTranscriptID)

    # Pull out transcript annotations in gtfDF
    transcriptDF = gtfDF[gtfDF[2] == 'transcript']
    exonDF = gtfDF[gtfDF[2] == 'exon']

    col_names = ['transcript_ID', 'chrom', 'strand', 'exonStart', 'exonEnd']
    # Double check that the GTF file only has annotations for two transcript isoforms
    if transcriptDF.shape[0] == 2:
        transcriptID = transcriptDF[9]

        # Establish an empty pandas dataframe for output
        inputDF = DataFrame(columns=col_names)
        is_minus_strand = False

        orig_ends_by_transcript = dict()
        # For each transcript in transcriptID, record exon start and end coordinates
        for transcript in transcriptID:
            # Determine strand associated with transcript
            strand = transcriptDF[transcriptDF[9] == transcript][6].item()

            # Determine chromosome associated with transcript
            chrom = transcriptDF[transcriptDF[9] == transcript][0].item()

            # Pull out rows of gtfDF corresponding to exons of given transcript
            transcriptExonDF = exonDF[exonDF[9] == transcript]

            # Establish sorted lists for exon start and end coordinates based on strand
            if strand == '+':
                exonStart, exonEnd = sort(transcriptExonDF[3]), sort(transcriptExonDF[4])
            else:
                is_minus_strand = True
                exonStart, exonEnd = -sort(-transcriptExonDF[4]), -sort(-transcriptExonDF[3])

            orig_ends_by_transcript[transcript] = (exonStart[0], exonEnd[-1])
            # Update inputDF with information about current transcript
            inputDF = concat([
                inputDF,
                DataFrame([[transcript, chrom, strand, exonStart.tolist(), exonEnd.tolist()]],
                          columns=col_names)
            ])

        # Smooth out transcript ends if terminal exons of both transcript isoforms are overlapping
        inputDF = SmoothEnds(inputDF, is_minus_strand)

        return {'inputDF': inputDF, 'orig_ends': orig_ends_by_transcript}

    else:
        sys.exit('ERROR: Input GTF file does not contain two transcript isoforms!')

def get_endpoint_event_and_coords(transcript_1, transcript_2,
                                  orig_ends_by_transcript, chrom, strand):
    event = 'ENDS'
    coord_parts = list()
    for transcript in [transcript_1, transcript_2]:
        start, end = orig_ends_by_transcript[transcript]
        coord_part = '{}:{}:{}:{}'.format(chrom, start, end, strand)
        coord_parts.append(coord_part)

    if len(set(coord_parts)) == 1:
        # ends are the same
        return None, None

    coord = ';'.join(coord_parts)
    return event, coord

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
    parse_gtf_result = ParseGTF(infile)
    inputDF = parse_gtf_result['inputDF']
    orig_ends_by_transcript = parse_gtf_result['orig_ends']

    # ASSUME: Both transcript isoforms are on the same chromosome and strand
    chrom, strand = inputDF['chrom'].unique()[0], inputDF['strand'].unique()[0]
    chrom = str(chrom)  # ensure chrom is parsed as a string

    # Build a splice graph
    spliceGraph = BuildSpliceGraph(inputDF)

    # Find all bubbles in spliceGraph
    bubbles = FindBubbles(spliceGraph)

    # Establish an empty pandas dataframe to describe each of the bubbles
    output_columns = ['transcript1', 'transcript2', 'event', 'coordinates']
    outputDF = DataFrame(columns=output_columns)

    transcript_1 = inputDF['transcript_ID'].iloc[0]
    transcript_2 = inputDF['transcript_ID'].iloc[1]
    # Iterate over bubbles identified in spliceGraph
    for bubble in bubbles:
        # For each bubble, identify the underlying event and report coordinates
        event, coord = ClassifyBubble(bubble)
        coord = ';'.join([chrom + ':' + x + ':' + strand for x in coord.split(';')])

        # Update outputDF with information about each bubble
        outputDF = concat([
            outputDF,
            DataFrame([[transcript_1, transcript_2, event, coord]],
                      columns=output_columns)
        ])

    if not bubbles:
        event, coord = get_endpoint_event_and_coords(
            transcript_1, transcript_2, orig_ends_by_transcript, chrom, strand)
        if event:
            outputDF = concat([
                outputDF,
                DataFrame([[transcript_1, transcript_2, event, coord]],
                          columns=output_columns)
            ])

    # Print contents of outputDF to outfile
    outputDF.to_csv(outfile, sep='\t', index=False)


if __name__ == '__main__':
    main()
