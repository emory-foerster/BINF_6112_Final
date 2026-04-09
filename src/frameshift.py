#!/usr/bin/env python3
#Mekhi Lucas
import orf
import sys
import fasta_io

file = '../datasets/Covid_GCF_009858895.2/sequence.fasta'
'''
Input:
Using the longest sequence from the ORF function \

Do:
Analyzing codons in 3 frames, each of the following offset by 1 codon \
Looking from a stop codon in frame 0 \
If frames 1 or 2 continue for 30+ amino acids = possible frameshift \

Output:
Return that there is a possible frameshift at that site (True or False)\
Return the site and type  (+1, +2) of the possible frameshift \
'''


class FrameshiftDetector:
    '''
    Purpose:
        All frameshift detection logic for a set of ORFs.
        Stores all_orfs and window as instance attributes so helper
        methods do not need them passed on every call.

    Input:
        all_orfs (list[dict]): All ORFs found across all reading frames
        window   (int):        Nucleotide window to search around an ORF's
                               end position when looking for neighbors —
                               default is 200

    Usage:
        detector = FrameshiftDetector(all_orfs, window=200)
        long_orf = detector.longest_orf()
        result   = detector.analyze(long_orf)
    '''

    def __init__(self, all_orfs: list[dict], window: int = 200):
        self.all_orfs = all_orfs
        self.window   = window

    def longest_orf(self) -> dict:
        '''
        Purpose:
            Find and return the longest ORF in self.all_orfs by length.

        Input:
            self

        Output:
            dict: The ORF dictionary with the greatest "length" value

        High-level steps:
            1. Initialize longest to the first ORF in self.all_orfs
            2. Iterate through all ORFs
            3. Update longest whenever a greater length is found
            4. Return longest
        '''
        longest = self.all_orfs[0]
        for o in self.all_orfs:
            if o['length'] > longest['length']:
                longest = o
        return longest

    def orf_coverage_proportion(self, long_orf: dict) -> float:
        '''
        Purpose:
            Calculate how dominant the longest ORF is relative to all ORFs
            found, as a proportion of total nucleotides across all ORFs.

        Input:
            long_orf (dict): The longest ORF dictionary
            self

        Output:
            float: long_orf length divided by total nucleotides across all ORFs

        High-level steps:
            1. Sum the "length" value across every dict in self.all_orfs
            2. Divide long_orf["length"] by that total
            3. Return the rounded proportion
        '''
        total = 0
        for o in self.all_orfs:
            total += o['length']
        return round(long_orf['length'] / total, 3)

    def find_neighboring_orfs(self, long_orf: dict) -> list:
        '''
        Purpose:
            Search self.all_orfs for ORFs in a different reading frame whose
            start position falls within self.window of the longest ORF's end
            position, indicating possible frameshift continuations.

        Input:
            long_orf (dict): The longest ORF dictionary
            self

        Output:
            A list of dicts, each being a neighboring ORF in a different frame
            whose start falls within self.window nucleotides of long_orf end —
            empty list if none found

        High-level steps:
            1. Get long_orf["end"] and long_orf["frame"] as reference points
            2. Iterate through self.all_orfs
                a. Skip any ORF that shares the same frame as long_orf
                b. Check if the candidate ORF's start falls within
                   long_orf["end"] ± self.window
                c. If yes, append that ORF dict to the neighbors list
            3. Return the neighbors list (empty if none found)
        '''
        end_pos   = long_orf['end']
        neighbors = []

        for o in self.all_orfs:
            if o['frame'] == long_orf['frame']:
                continue
            if end_pos - self.window <= o['start'] <= end_pos + self.window:
                neighbors.append(o)

        return neighbors

    def shift_type(self, long_orf: dict) -> list:
        '''
        Purpose:
            Determine the type and magnitude of the frameshift between the
            longest ORF and each of its neighboring ORFs.

        Input:
            long_orf (dict): The longest ORF dictionary
            self

        Output:
            A list of tuples, one per neighboring ORF, each containing:
                - shift_type      (str): Frameshift type as a string e.g. "+1" or "+2"
                - shift_magnitude (int): Raw integer shift value (1 or 2)
            Returns None if no neighboring ORFs are found

        High-level steps:
            1. Call find_neighboring_orfs to get all neighboring ORFs
            2. If no neighbors found, return None
            3. For each neighbor, subtract long_orf["frame"] from neighbor["frame"]
            4. Apply modulo 3 to handle frame wraparound
            5. Format as a "+N" string for shift_type
            6. Append (shift_str, shift_magnitude) tuple to shifts list
            7. Return shifts list
        '''
        neighbors = self.find_neighboring_orfs(long_orf)

        if not neighbors:
            return None

        shifts = []
        for neighbor in neighbors:
            shift_magnitude = (neighbor['frame'] - long_orf['frame']) % 3
            shift_str       = '+' + str(shift_magnitude)
            shifts.append((shift_str, shift_magnitude))

        return shifts

    def build_frameshift_details(self, long_orf: dict, neighboring_orfs: list[dict]) -> list[dict]:
        """
        Purpose:
            Assemble a list of frameshift detail dicts, one per neighboring ORF,
            combining the shift position, type, and neighboring ORF data into
            structured output.

        Input:
            long_orf         (dict):       The longest ORF dictionary
            neighboring_orfs (list[dict]): Neighboring ORFs from find_neighboring_orfs
            self

        Output:
            A list of dicts, one per neighboring ORF, each containing:
                - "shift_position"   (int):  The neighboring ORF's start position
                - "neighboring_frame"(int):  The frame of the neighboring ORF
                - "shift_type"       (str):  "+1" or "+2"
                - "shift_magnitude"  (int):  1 or 2
                - "neighboring_orf"  (dict): The full neighboring ORF dict

        High-level steps:
            1. Call shift_type to get a list of (shift_str, shift_magnitude) tuples
            2. Zip neighboring_orfs with shift_data
            3. For each (neighbor, shift) pair, build a details dict
            4. Append each details dict to details_list
            5. Return details_list
        """
        shift_data   = self.shift_type(long_orf)
        details_list = []

        for neighbor, (shift_str, shift_magnitude) in zip(neighboring_orfs, shift_data):
            details = {
                'shift_position'   : neighbor['start'],
                'neighboring_frame': neighbor['frame'],
                'shift_type'       : shift_str,
                'shift_magnitude'  : shift_magnitude,
                'neighboring_orf'  : neighbor
            }
            details_list.append(details)

        return details_list

    def analyze(self, long_orf: dict) -> dict:
        '''
        Purpose:
            Orchestrate the full frameshift analysis of the longest ORF by
            calling each helper method and assembling the final result.

        Input:
            long_orf (dict): The longest ORF dictionary from longest_orf()
            self

        Output:
            The long_orf dictionary with three keys added:
                - "dominance_ratio"    (float):      Proportion of total ORF
                                                     nucleotides in this ORF
                - "frameshift_boolean" (bool):       True if neighboring ORFs found
                - "frameshift_details" (list[dict]): One detail dict per frameshift
                                                     candidate, or None if none found

        High-level steps:
            1. Validate long_orf is not empty — raise ValueError if it is
            2. Call orf_coverage_proportion and store as dominance_ratio
            3. Call find_neighboring_orfs to get all frameshift candidates
            4. If neighboring ORFs are found:
                a. Set frameshift_flag to True
                b. Call build_frameshift_details to assemble list of frameshift_details
            5. If no neighboring ORFs are found:
                a. Set frameshift_boolean to False
                b. Set frameshift_details to None
            6. Add dominance_ratio, frameshift_boolean, and frameshift_details
               to long_orf and return it
        '''
        if not long_orf:
            raise ValueError("long_orf is empty")

        neighboring_orfs = self.find_neighboring_orfs(long_orf)

        if neighboring_orfs:
            is_frameshift      = True
            frameshift_details = self.build_frameshift_details(long_orf, neighboring_orfs)
        else:
            is_frameshift      = False
            frameshift_details = None

        long_orf['dominance_ratio']    = self.orf_coverage_proportion(long_orf)
        long_orf['frameshift_boolean'] = is_frameshift
        long_orf['frameshift_details'] = frameshift_details

        return long_orf


if __name__ == '__main__':
    fasta_dict = fasta_io.read_fasta(file)
    seq        = fasta_dict[0]['Sequence']
    all_orfs   = orf.detect_all_frames(seq)

    detector = FrameshiftDetector(all_orfs, window=200)
    long_orf = detector.longest_orf()
    result   = detector.analyze(long_orf)
    print(result)
