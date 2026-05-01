02/26/2026 \
Function: frameshift.py \
Group member: Mekhi Lucas \



# Example Input

List of dictionaries: \
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


#Code for the Function

def longest_orf(all_orfs: list):
    top = 0 
    longest = all_orfs[0]
    for orf in all_orfs:
            if orf['length'] > longest['length']:
                longest = orf
    return longest
print(longest_orf(all_orfs))

# Output as:

This will return the dictionary of the longest ORF. To be used for frameshift detection