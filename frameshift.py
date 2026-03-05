#!/usr/bin/env python3
#Mekhi Lucas
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
all_orfs = [
    {'ORF_1': 'seq',
        'start': 10,
        'end': 150,
        'length': int(140),
        'frame': 0
    },
    {'ORF_2': 2,
        'start': 25,
        'end': 380,
        'length': int(355),
        'frame': 1
    },
    {'ORF_3': 3,
        'start': 5,
        'end': 200,
        'length': int(195),
        'frame': 2
    },
    {'ORF_4': 4,
        'start': 60,
        'end': 500,
        'length': int(440),
        'frame': 0
    },
    {'ORF_5': 5,
        'start': 120,
        'end': 280,
        'length': int(160),
        'frame': 1
    }
]


import sys
print(type(all_orfs))
def longest_orf(all_orfs: list):
    top = 0 
    longest = all_orfs[0]
    for orf in all_orfs:
            if orf['length'] > longest['length']:
                longest = orf
    return longest
print(longest_orf(all_orfs))


# def ORF_coverage_proportion(long_orf:dict, all_orfs:list[dict]):
# 	'''
# 	Purpose:
#         Calculate how dominant the longest ORF is relative to all ORFs found.

#     Input:
#         long_orf  (dict):       The longest ORF dictionary from detect_ORF
#         all_orfs  (list[dict]): All ORFs found across all three reading frames

#     Output:
#         float: The longest ORF length divided by total nucleotides across all ORFs

#     High-level steps:
#         1. Sum the "length" value across every dict in all_orfs
#         2. Divide long_orf["length"] by that total
# 	'''
# 	len_long_orf = len(long_orf[seq])
# 	frame = long_orf[frame]

# 	for orf in all_orfs:
# 		if orf[frame] == frame:
# 			total += orf[length]

# 	coverage_proportion = round(len_long_orf/ total, 3)

# 	return long_orf[proportion] = coverage_proportion #I think this works


# def find_neighboring_orf(long_orf: dict, all_orfs: list[dict], window: int = 15):
# 	'''
# 	Purpose:
#         Search all_orfs for an ORF in a different reading frame whose start
#         position falls within a window of the longest ORF's end position,
#         indicating a possible frameshift continuation.

#     Input:
#         long_orf  (dict):       The longest ORF dictionary
#         all_orfs  (list[dict]): All ORFs found across all three reading frames
#         window    (int):        Nucleotide window to search around long_orf end
#                                 position — default is 15

#     Output:
#         A dictionary of the neighboring ORF if one is found, or None if not

#     High-level steps:
#         1. Get long_orf["end"] and long_orf["frame"] as reference points
#         2. Iterate through all_orfs
#             a. Skip any ORF that shares the same frame as long_orf
#             b. Check if the candidate ORF's start falls within
#                long_orf["end"] ± window
#             c. If yes, return that ORF dict immediately
#         3. Return None if no neighboring ORF is found
# 	'''

	
# def shift_type(long_orf:dict, all_orfs:list[dict]):
#  '''
#  Purpose:
#         Determine the type and magnitude of the frameshift between the longest
#         ORF and its neighboring ORF.

#     Input:
#         long_orf       (dict): The longest ORF dictionary
#         neighboring_orf(dict): The neighboring ORF found by find_neighboring_orf

#     Output:
#         A tuple containing:
#             - shift_type      (str): The frameshift type as a string e.g. "+1" or "+2"
#             - shift_magnitude (int): The raw integer shift value (1 or 2)

#     High-level steps:
#         1. Subtract long_orf["frame"] from neighboring_orf["frame"]
#         2. Apply modulo 3 to handle frame wraparound
#         3. Format as a "+N" string for the shift_type
#         4. Return both values as a tuple
#  '''


# def build_frameshift_details(long_orf: dict, neighboring_orf: dict):
#   """
#     Purpose:
#         Assemble the frameshift_details dictionary once a neighboring ORF
#         has been confirmed, by combining the shift position, type, and
#         neighboring ORF data into one structured output.

#     Input:
#         long_orf        (dict): The longest ORF dictionary
#         neighboring_orf (dict): The neighboring ORF found by find_neighboring_orf

#     Output:
#         A dictionary containing:
#             - "shift_position"   (int): long_orf end position as the approximate
#                                         frameshift location
#             - "neighboring_frame"(int): The frame of the neighboring ORF
#             - "shift_type"       (str): "+1" or "+2"
#             - "shift_magnitude"  (int): 1 or 2
#             - "neighboring_orf"  (dict): The full neighboring ORF dict

#     High-level steps:
#         1. Call calculate_shift_type to get shift_type and shift_magnitude
#         2. Build and return the details dictionary using those values
#            plus long_orf["end"] and neighboring_orf["frame"]
#     """
    



# def frameshift_detector_longest_ORF(long_orf: dict, all_orfs: list[dict]):
# 	'''
# 	 Purpose:
#         Orchestrate the full frameshift analysis of the longest ORF by calling
#         each modular helper function and assembling the final result.

#     Input:
#         long_orf  (dict):       The longest ORF dictionary from detect_ORF
#         all_orfs  (list[dict]): All ORFs found across all three reading frames

#     Output:
#         The long_orf dictionary gets:
#             - dominance_ratio
#             - frameshift_boolean
#             - frameshift_details (type)

#     High-level steps:
#         1. Validate long_orf is not empty — raise ValueError if it is
#         2. Call calculate_dominance_ratio and store the result
#         3. Call find_neighboring_orf to search for a frameshift candidate
#         4. If a neighboring ORF is found:
#             a. Set frameshift_flag to True
#             b. Call build_frameshift_details to assemble frameshift_details
#         5. If no neighboring ORF is found:
#             a. Set frameshift_boolean to False
#             b. Set frameshift_details to None
#         6. Add dominance_ratio, frameshift_flag, and frameshift_details
#            to long_orf and return it

# 	'''

# def Frame_detect(ORF, threshold):
#   #detection across all open reading frames
	


# if __name__ == '__main__':
# 	main()



